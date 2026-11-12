#!/usr/bin/env python3
"""Align pairs of selector axes and write combined molecular structures."""

import argparse
import csv
import json
import math
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import yaml
from rdkit import Chem
from rdkit.Geometry import Point3D

from alignment.polar_scan import PositionMerger, spherical_position
from alignment.reference_frame import ReferenceSelector
from core.logging import configure_logging, get_logger
from selection.axis_config import UniqueKeyLoader
from selection.matching import MoleculeSelector

logger = get_logger("alignment.aligner")
FILTER_FIELDS = ("molecule_name", "selector", "axis_type")
REQUIRED_FIELDS = (*FILTER_FIELDS, "sdf_file", "record_index", "center_atom_index", "center", "axis")


@dataclass(frozen=True)
class AxisRow:
    line: int
    fields: dict
    path: Path
    record_index: int
    center_index: int | None
    center: np.ndarray
    axis: np.ndarray


class AxisInput:
    """Load selected axes and their referenced source molecules."""

    def __init__(self):
        self.molecules = {}

    @staticmethod
    def vector(value, field):
        try:
            data = json.loads(value)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{field} must be a JSON array of three numbers") from exc
        if (not isinstance(data, list) or len(data) != 3
                or any(type(x) not in (int, float) or not math.isfinite(x) for x in data)):
            raise ValueError(f"{field} must be a JSON array of three finite numbers")
        return np.asarray(data, dtype=float)

    def read(self, path, rigid, flexible):
        left, right = [], []
        with Path(path).open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, strict=True)
            if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise ValueError("CSV must have unique headers")
            missing = set(REQUIRED_FIELDS) - set(reader.fieldnames)
            if missing:
                raise ValueError(f"CSV missing required columns: {', '.join(sorted(missing))}")
            for row in reader:
                line = reader.line_num
                if None in row or any(row[key] is None for key in reader.fieldnames):
                    raise ValueError(f"CSV line {line}: incorrect number of columns")
                is_left = rigid.matches(row)
                is_right = flexible.matches(row)
                if not (is_left or is_right):
                    continue
                try:
                    sdf = Path(row["sdf_file"])
                    if not row["sdf_file"] or not sdf.is_file():
                        raise ValueError(f"SDF file does not exist: {sdf}")
                    record_index = int(row["record_index"])
                    if record_index < 0:
                        raise ValueError("record_index must be non-negative")
                    center_index = int(row["center_atom_index"]) if row["center_atom_index"] else None
                    center = self.vector(row["center"], "center")
                    axis = self.vector(row["axis"], "axis")
                    length = np.linalg.norm(axis)
                    if not math.isfinite(length) or length < 1e-12:
                        raise ValueError("axis must have nonzero finite length")
                    sdf = sdf.resolve()
                    if sdf not in self.molecules:
                        self.molecules[sdf] = list(MoleculeSelector._read_sdf(sdf))
                    molecules = self.molecules[sdf]
                    if record_index >= len(molecules):
                        raise ValueError("record_index is outside the SDF")
                    mol = molecules[record_index]
                    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                        raise ValueError("SDF record must have 3D coordinates")
                    if center_index is not None:
                        if center_index < 0 or center_index >= mol.GetNumAtoms():
                            raise ValueError("center_atom_index is outside the SDF")
                        point = mol.GetConformer().GetAtomPosition(center_index)
                        if not np.allclose(center, (point.x, point.y, point.z), rtol=0, atol=1e-3):
                            raise ValueError("center does not match the SDF atom coordinates")
                    elif row["axis_type"] != "aromatic_normal":
                        raise ValueError("center_atom_index is required except for aromatic ring axes")
                    entry = AxisRow(line, dict(row), sdf, record_index, center_index, center, axis / length)
                except (ValueError, OSError, RuntimeError) as exc:
                    raise ValueError(f"CSV line {line}: {exc}") from exc
                if is_left:
                    left.append(entry)
                if is_right:
                    right.append(entry)
        if not left or not right:
            raise ValueError(f"filters matched {len(left)} rigid and {len(right)} flexible axes")
        return left, right

    def molecule(self, row):
        return self.molecules[row.path][row.record_index]


class AxisFilter:
    def __init__(self, expressions):
        if isinstance(expressions, dict):
            expressions = [expressions]
        if not isinstance(expressions, list) or not expressions:
            raise ValueError("each side requires at least one filter triplet")
        self.filters = []
        for item in expressions:
            if not isinstance(item, dict) or set(item) != set(FILTER_FIELDS):
                raise ValueError("each filter must contain molecule_name, selector and axis_type")
            try:
                self.filters.append({key: re.compile(item[key]) for key in FILTER_FIELDS})
            except (re.error, TypeError) as exc:
                raise ValueError(f"invalid filter regex: {exc}") from exc

    def matches(self, row):
        return any(all(pattern.search(row[key]) is not None for key, pattern in group.items())
                   for group in self.filters)


class PairAligner:
    @staticmethod
    def rotation(source, target):
        cross = np.cross(source, target)
        sine = np.linalg.norm(cross)
        cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
        if sine < 1e-12:
            if cosine > 0:
                return np.eye(3)
            basis = np.eye(3)[np.argmin(np.abs(source))]
            axis = np.cross(source, basis)
            axis /= np.linalg.norm(axis)
            return 2 * np.outer(axis, axis) - np.eye(3)
        axis = cross / sine
        skew = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + skew * sine + skew @ skew * (1 - cosine)

    @staticmethod
    def radius(row, mol):
        if row.center_index is None:
            return 0.0
        atom = mol.GetAtomWithIdx(row.center_index)
        if atom.GetAtomicNum() == 0:
            raise ValueError(f"CSV line {row.line}: a pseudo-atom has no van der Waals radius")
        radius = Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError(f"CSV line {row.line}: no valid van der Waals radius for {atom.GetSymbol()}")
        return float(radius)

    def align(self, rigid, flexible, rigid_mol, flexible_mol, distance, offset=None,
              orientation="stable"):
        if distance is None:
            distance = self.radius(rigid, rigid_mol) + self.radius(flexible, flexible_mol)
        matrix = self.rotation(flexible.axis, -rigid.axis)
        target = rigid.center + distance * rigid.axis
        if offset is not None:
            target = target + offset
        if orientation == "toward-center":
            direction = rigid.center - target
            length = np.linalg.norm(direction)
            if not math.isfinite(length):
                raise ValueError("non-finite center separation")
            if length > 1e-12:
                matrix = self.rotation(-rigid.axis, direction / length) @ matrix
        translation = target - matrix @ flexible.center
        moved = Chem.Mol(flexible_mol)
        conformer = moved.GetConformer()
        for index in range(moved.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            xyz = matrix @ np.array((point.x, point.y, point.z)) + translation
            if not np.isfinite(xyz).all():
                raise ValueError("non-finite coordinates in flexible molecule")
            conformer.SetAtomPosition(index, Point3D(*xyz))
        rigid_conf = rigid_mol.GetConformer()
        for index in range(rigid_mol.GetNumAtoms()):
            point = rigid_conf.GetAtomPosition(index)
            if not all(math.isfinite(value) for value in (point.x, point.y, point.z)):
                raise ValueError("non-finite coordinates in rigid molecule")
        pair = Chem.CombineMols(rigid_mol, moved)
        return pair, matrix, translation, target, distance


@dataclass(frozen=True)
class DistancePlan:
    mode: str
    start: Fraction | None = None
    step: Fraction | None = None
    first: int = 0
    count: int = 1
    distance: float | None = None

    @staticmethod
    def number(value, name):
        if type(value) not in (int, float, str):
            raise ValueError(f"{name} must be a finite number")
        try:
            number = Decimal(str(value))
        except InvalidOperation as exc:
            raise ValueError(f"{name} must be a finite number") from exc
        if not number.is_finite() or not math.isfinite(float(number)):
            raise ValueError(f"{name} must be a finite number")
        return Fraction(number)

    @classmethod
    def from_options(cls, settings, args, allow_negative=False, default_distance=None):
        def option(name, default=None):
            value = getattr(args, name)
            return value if value is not None else settings.get(name, default)

        mode = option("alignment_mode", "default")
        distance = option("distance")
        scan_names = ("scan_start", "scan_end", "include_start", "include_end", "scan_step", "scan_count")
        scan = {name: option(name) for name in scan_names}
        if mode not in ("default", "distance-scan"):
            raise ValueError("alignment_mode must be 'default' or 'distance-scan'")
        if mode == "default":
            if any(value is not None for value in scan.values()):
                raise ValueError("scan options require alignment_mode distance-scan")
            if distance is None:
                return cls(mode, distance=default_distance)
            value = cls.number(distance, "distance")
            if value < 0 and not allow_negative:
                raise ValueError("distance must be non-negative")
            return cls(mode, distance=float(value))
        if distance is not None:
            raise ValueError("distance cannot be used with distance-scan")
        if scan["scan_start"] is None or scan["scan_end"] is None:
            raise ValueError("distance-scan requires scan_start and scan_end")
        start = cls.number(scan["scan_start"], "scan_start")
        end = cls.number(scan["scan_end"], "scan_end")
        if (start < 0 and not allow_negative) or end <= start:
            raise ValueError("scan range must have scan_start < scan_end" +
                             ("" if allow_negative else " and 0 <= scan_start"))
        for name in ("include_start", "include_end"):
            if scan[name] is None:
                scan[name] = True
            elif type(scan[name]) is not bool:
                raise ValueError(f"{name} must be a boolean")
        # A CLI sampling choice overrides the alternative from the config.
        if args.scan_step is not None:
            scan["scan_count"] = None
        elif args.scan_count is not None:
            scan["scan_step"] = None
        if (scan["scan_step"] is None) == (scan["scan_count"] is None):
            raise ValueError("distance-scan requires exactly one of scan_step and scan_count")
        include_start, include_end = scan["include_start"], scan["include_end"]
        first = 0 if include_start else 1
        if scan["scan_step"] is not None:
            step = cls.number(scan["scan_step"], "scan_step")
            if step <= 0:
                raise ValueError("scan_step must be positive")
            intervals = (end - start) / step
            if intervals.denominator != 1:
                raise ValueError("scan range must be evenly divisible by scan_step")
            count = intervals.numerator + 1 - int(not include_start) - int(not include_end)
            if count <= 0:
                raise ValueError("scan produces no poses")
        else:
            interior_count = scan["scan_count"]
            if type(interior_count) is not int or interior_count < 1:
                raise ValueError("scan_count must be a positive integer")
            intervals = interior_count + 1
            count = interior_count + int(include_start) + int(include_end)
            step = (end - start) / intervals
        return cls(mode, start, step, first, count)

    def distances(self):
        if self.mode == "default":
            yield self.distance
            return
        for index in range(self.first, self.first + self.count):
            yield float(self.start + index * self.step)


def load_config(path):
    if path is None:
        return {}
    try:
        with Path(path).open(encoding="utf-8") as source:
            data = yaml.load(source, Loader=UniqueKeyLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"invalid configuration {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("configuration must be a YAML mapping")
    unknown = set(data) - {"rigid", "flexible", "output", "output_format", "overwrite_output",
                           "grid_based_scan", "polar_coordinate_based_scan",
                           "reference_axis", "flexible_orientation"}
    if unknown:
        raise ValueError(f"unknown configuration keys: {', '.join(sorted(map(str, unknown)))}")
    return data


def triplets(values):
    return [dict(zip(FILTER_FIELDS, group)) for group in values]


def axis_plan(axis, scans, args, *, polar=False):
    scan_keys = ("scan_start", "scan_end", "include_start", "include_end", "scan_step", "scan_count")
    keys = ("distance", *scan_keys)
    values = scans.get(axis, {})
    section = "polar_coordinate_based_scan" if polar else "grid_based_scan"
    if not isinstance(values, dict) or set(values) - set(keys):
        raise ValueError(f"{section}.{axis} must be a mapping of distance and scan options")
    prefix = f"polar_{axis}_" if polar else ("" if axis == "y" else f"{axis}_")
    options = argparse.Namespace(**{key: getattr(args, prefix + key) for key in ("alignment_mode", *keys)})
    requested = bool(values) or any(getattr(options, key) is not None for key in ("alignment_mode", *keys))
    values = dict(values)
    cli_scan = any(getattr(options, key) is not None for key in scan_keys)
    if options.distance is not None and not cli_scan and options.alignment_mode != "distance-scan":
        values = {key: value for key, value in values.items() if key not in scan_keys}
    elif options.alignment_mode == "distance-scan" or cli_scan:
        values.pop("distance", None)
    if options.alignment_mode == "default":
        values = {key: value for key, value in values.items() if key not in scan_keys}
    values["alignment_mode"] = "distance-scan" if any(key in values for key in scan_keys) or cli_scan else "default"
    try:
        plan = DistancePlan.from_options(values, options, allow_negative=not polar and axis != "y",
                                         default_distance=None if axis in ("y", "radius") else 0.0)
    except ValueError as exc:
        raise ValueError(f"{section}.{axis}: {exc}") from exc
    return plan, requested


def polar_options(scans, args):
    if not isinstance(scans, dict) or set(scans) - {"radius", "phi", "theta", "phi_unit", "theta_unit", "merge_tolerance"}:
        raise ValueError("polar_coordinate_based_scan must contain only radius, phi, theta, angle units and merge_tolerance")
    units = []
    plans = []
    for angle in ("phi", "theta"):
        unit = getattr(args, f"polar_{angle}_unit") or scans.get(f"{angle}_unit", "degrees")
        if unit not in ("degrees", "radians"):
            raise ValueError(f"polar {angle}_unit must be degrees or radians")
        plan, _ = axis_plan(angle, scans, args, polar=True)
        maximum = (360 if angle == "phi" else 180) if unit == "degrees" else (2 * math.pi if angle == "phi" else math.pi)
        bounds = [plan.distance] if plan.mode == "default" else [float(plan.start), float(plan.start + plan.step * (plan.count - 1 + plan.first))]
        # Validate declared boundaries too, including excluded endpoints.
        if plan.mode == "distance-scan":
            raw_end = getattr(args, f"polar_{angle}_scan_end")
            bounds.append(float(DistancePlan.number(raw_end if raw_end is not None else scans[angle]["scan_end"], "scan_end")))
        if any(value is None or value < 0 or value > maximum + 1e-12 for value in bounds):
            raise ValueError(f"polar {angle} must be between 0 and {maximum} {unit}")
        units.append(unit)
        plans.append(plan)
    tolerance = args.polar_merge_tolerance if args.polar_merge_tolerance is not None else scans.get("merge_tolerance", 0.05)
    tolerance = float(DistancePlan.number(tolerance, "merge_tolerance"))
    if tolerance <= 0:
        raise ValueError("polar merge_tolerance must be positive")
    radius, _ = axis_plan("radius", scans, args, polar=True)
    return radius, plans[0], plans[1], units[0], units[1], tolerance


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path, help="alignment CSV produced by selector.py")
    parser.add_argument("--config", type=Path, help="YAML parameters (input CSV is always positional)")
    parser.add_argument("--rigid", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--flexible", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--output", "-o", type=Path, help="new pair directory or new multi-record SDF")
    parser.add_argument("--output-format", choices=("pair-files", "single"), help="SDF layout (default: pair-files)")
    parser.add_argument("--overwrite-output", action=argparse.BooleanOptionalAction, default=None,
                        help="replace matching output files (default: false)")
    parser.add_argument("--flexible-orientation", choices=("stable", "toward-center"),
                        help="orientation after placement (default: stable)")
    parser.add_argument("--coordinate-system", choices=("grid", "polar"))
    parser.add_argument("--distance", type=float, help="fixed center separation in Å; default: sum of van der Waals radii")
    parser.add_argument("--alignment-mode", choices=("default", "distance-scan"))
    parser.add_argument("--scan-start", help="scan lower boundary in Å")
    parser.add_argument("--scan-end", help="scan upper boundary in Å")
    parser.add_argument("--include-start", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--include-end", action=argparse.BooleanOptionalAction, default=None)
    sampling = parser.add_mutually_exclusive_group()
    sampling.add_argument("--scan-step", help="fixed scan step in Å")
    sampling.add_argument("--scan-count", type=int, help="number of interior poses per axis pair; boundaries are additional")
    for axis in ("x", "z"):
        parser.add_argument(f"--{axis}-alignment-mode", choices=("default", "distance-scan"))
        parser.add_argument(f"--{axis}-distance", help=f"fixed signed {axis} offset in Å")
        parser.add_argument(f"--{axis}-scan-start", help=f"{axis} scan lower boundary in Å")
        parser.add_argument(f"--{axis}-scan-end", help=f"{axis} scan upper boundary in Å")
        parser.add_argument(f"--{axis}-include-start", action=argparse.BooleanOptionalAction, default=None)
        parser.add_argument(f"--{axis}-include-end", action=argparse.BooleanOptionalAction, default=None)
        lateral_sampling = parser.add_mutually_exclusive_group()
        lateral_sampling.add_argument(f"--{axis}-scan-step", help=f"{axis} scan step in Å")
        lateral_sampling.add_argument(f"--{axis}-scan-count", type=int, help=f"number of interior {axis} poses")
    for axis in ("radius", "phi", "theta"):
        parser.add_argument(f"--polar-{axis}-alignment-mode", choices=("default", "distance-scan"))
        parser.add_argument(f"--polar-{axis}", dest=f"polar_{axis}_distance", help=f"fixed polar {axis}")
        parser.add_argument(f"--polar-{axis}-scan-start")
        parser.add_argument(f"--polar-{axis}-scan-end")
        parser.add_argument(f"--polar-{axis}-include-start", action=argparse.BooleanOptionalAction, default=None)
        parser.add_argument(f"--polar-{axis}-include-end", action=argparse.BooleanOptionalAction, default=None)
        sampling_group = parser.add_mutually_exclusive_group()
        sampling_group.add_argument(f"--polar-{axis}-scan-step")
        sampling_group.add_argument(f"--polar-{axis}-scan-count", type=int)
    parser.add_argument("--polar-phi-unit", choices=("degrees", "radians"))
    parser.add_argument("--polar-theta-unit", choices=("degrees", "radians"))
    parser.add_argument("--polar-merge-tolerance", help="duplicate position tolerance in Å (default: 0.05)")
    args = parser.parse_args(argv)
    configure_logging()
    try:
        settings = load_config(args.config)
        orientation = args.flexible_orientation if args.flexible_orientation is not None else settings.get("flexible_orientation", "stable")
        if orientation not in ("stable", "toward-center"):
            raise ValueError("flexible_orientation must be 'stable' or 'toward-center'")
        rigid = AxisFilter(triplets(args.rigid) if args.rigid is not None else settings.get("rigid"))
        flexible = AxisFilter(triplets(args.flexible) if args.flexible is not None else settings.get("flexible"))
        output_format = args.output_format or settings.get("output_format", "pair-files")
        if output_format not in ("pair-files", "single"):
            raise ValueError("output_format must be 'pair-files' or 'single'")
        output_value = args.output or settings.get("output", "aligned_pairs" if output_format == "pair-files" else "aligned.sdf")
        if not isinstance(output_value, (str, Path)) or not str(output_value):
            raise ValueError("output must be a nonempty path")
        output = Path(output_value)
        if output.name in ("", ".", ".."):
            raise ValueError("output must name a new directory or file")
        csv_output = output.with_suffix(".csv") if output_format == "single" else output.parent / (output.name + ".csv")
        if output == csv_output:
            raise ValueError("SDF output and CSV output must be different paths")
        grid_present = "grid_based_scan" in settings
        polar_present = "polar_coordinate_based_scan" in settings
        grid_cli_names = ("alignment_mode", "distance", "scan_start", "scan_end", "include_start",
                          "include_end", "scan_step", "scan_count")
        grid_cli = any(getattr(args, prefix + name) is not None
                       for prefix in ("", "x_", "z_") for name in grid_cli_names)
        polar_cli = any(value is not None for name, value in vars(args).items() if name.startswith("polar_"))
        if (grid_present or grid_cli) and (polar_present or polar_cli):
            raise ValueError("grid and polar scan settings cannot be combined")
        coordinate_system = args.coordinate_system or ("polar" if polar_present or polar_cli else "grid")
        if (coordinate_system == "polar" and (grid_present or grid_cli)) or (
                coordinate_system == "grid" and (polar_present or polar_cli)):
            raise ValueError("coordinate-system conflicts with configured scan settings")
        if coordinate_system == "grid":
            scans = settings.get("grid_based_scan", {})
            if not isinstance(scans, dict) or set(scans) - {"x", "y", "z"}:
                raise ValueError("grid_based_scan must be a mapping containing only x, y and z")
            plan, _ = axis_plan("y", scans, args)
            x_plan, x_requested = axis_plan("x", scans, args)
            z_plan, z_requested = axis_plan("z", scans, args)
        else:
            polar_config = settings.get("polar_coordinate_based_scan", {})
            radius_plan, phi_plan, theta_plan, phi_unit, theta_unit, tolerance = polar_options(polar_config, args)
        reference_selector = ReferenceSelector(settings.get("reference_axis"))
        overwrite = args.overwrite_output if args.overwrite_output is not None else settings.get("overwrite_output", False)
        if type(overwrite) is not bool:
            raise ValueError("overwrite_output must be a boolean")
        if not overwrite and (os.path.lexists(output) or os.path.lexists(csv_output)):
            raise ValueError(f"output already exists: {output} or {csv_output}")
        reader = AxisInput()
        rigid_rows, flexible_rows = reader.read(args.input_csv, rigid, flexible)
        needs_frame = coordinate_system == "polar" or x_requested or z_requested
        frames = {row.line: reference_selector.select(row, reader.molecule(row))
                  for row in rigid_rows} if needs_frame else {}
        if output_format == "pair-files":
            destinations = [output / f"axis-pair-{index:06d}.sdf"
                            for index in range(1, len(rigid_rows) * len(flexible_rows) + 1)]
            if os.path.lexists(output) and (output.is_symlink() or not output.is_dir()):
                raise ValueError("pair-files output must be a directory, not a file or symlink")
        else:
            destinations = [output]
        destinations.append(csv_output)
        inputs = [args.input_csv, *reader.molecules]
        if args.config is not None:
            inputs.append(args.config)
        for destination in destinations:
            if destination.is_dir():
                raise ValueError(f"output file is a directory: {destination}")
            for source in inputs:
                if destination.resolve() == source.resolve() or (
                        destination.exists() and source.exists() and os.path.samefile(destination, source)):
                    raise ValueError(f"output would overwrite an input file: {destination}")
        aligner = PairAligner()
        def generate_pairs():
            axis_pair = 0
            for r in rigid_rows:
                frame = frames.get(r.line)
                for f in flexible_rows:
                    axis_pair += 1
                    pose_index = 0
                    if coordinate_system == "polar":
                        merger = PositionMerger(tolerance)
                        duplicates = 0
                        def placements():
                            nonlocal duplicates
                            for radius in radius_plan.distances():
                                if radius is None:
                                    radius = (aligner.radius(r, reader.molecule(r))
                                              + aligner.radius(f, reader.molecule(f)))
                                for phi in phi_plan.distances():
                                    for theta in theta_plan.distances():
                                        point = spherical_position(radius, phi, theta, phi_unit, theta_unit)
                                        if not merger.add((point.x, point.y, point.z)):
                                            duplicates += 1
                                            continue
                                        yield point.y, point.x, point.z, point
                        positions = placements()
                    else:
                        def grid_placements():
                            for distance in plan.distances():
                                for x in x_plan.distances():
                                    for z in z_plan.distances():
                                        yield distance, x, z, None
                        positions = grid_placements()
                    for distance, x, z, point in positions:
                        pose_index += 1
                        offset = None if frame is None else x * frame.x + z * frame.z
                        mol, matrix, translation, target, separation = aligner.align(
                            r, f, reader.molecule(r), reader.molecule(f), distance, offset, orientation
                        )
                        yield (axis_pair, pose_index, r, f, mol, matrix, translation, target,
                               separation, x, z, frame, orientation, point,
                               phi_unit if point is not None else "", theta_unit if point is not None else "")
                    if coordinate_system == "polar" and duplicates:
                        logger.warning("axis pair %d: merged %d polar positions within %g Å",
                                       axis_pair, duplicates, tolerance)

        pairs = generate_pairs()
        count = write_output(output, output_format, pairs, destinations, overwrite)
        logger.info("wrote %d aligned poses and %s", count, csv_output)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


def write_output(output, output_format, pairs, destinations, overwrite):
    fields = ("pair_id", "axis_pair_id", "pose_index", "sdf_file", "sdf_record_index", "rigid_csv_line", "flexible_csv_line",
              "rigid_source", "flexible_source", "distance_angstrom", "rigid_center",
              "flexible_original_center", "flexible_aligned_center", "rigid_axis",
              "flexible_original_axis", "flexible_aligned_axis", "rotation", "translation",
              "x_offset_angstrom", "z_offset_angstrom", "x_axis", "z_axis",
              "reference_method", "reference_atom_index", "reference_point", "flexible_orientation",
              "coordinate_system", "polar_radius_angstrom", "polar_phi", "polar_phi_unit", "polar_theta", "polar_theta_unit")
    staging = Path(tempfile.mkdtemp(prefix=".aligner-", dir=output.parent))
    preserve_staging = False
    writer = None
    try:
        with (staging / "metadata.csv").open("x", encoding="utf-8", newline="") as metadata:
            if output_format == "single":
                writer = Chem.SDWriter(str(staging / "single.sdf"))
            current_axis_pair = None
            table = csv.DictWriter(metadata, fieldnames=fields, lineterminator="\n")
            table.writeheader()
            index = 0
            for index, (axis_pair, pose_index, rigid, flexible, molecule, matrix, translation, target,
                        distance, x, z, frame, orientation, point, phi_unit, theta_unit) in enumerate(pairs, 1):
                pair_id = f"pair-{index:06d}"
                axis_pair_id = f"axis-pair-{axis_pair:06d}"
                sdf = output / f"{axis_pair_id}.sdf" if output_format == "pair-files" else output
                molecule.SetProp("_Name", pair_id)
                molecule.SetProp("quick_qm_pair_id", pair_id)
                molecule.SetProp("quick_qm_axis_pair_id", axis_pair_id)
                molecule.SetProp("quick_qm_pose_index", str(pose_index))
                molecule.SetProp("quick_qm_distance_angstrom", str(distance))
                molecule.SetProp("quick_qm_flexible_orientation", orientation)
                molecule.SetProp("quick_qm_coordinate_system", "polar" if point is not None else "grid")
                if point is not None:
                    molecule.SetProp("quick_qm_polar_radius_angstrom", str(point.radius))
                    molecule.SetProp("quick_qm_polar_phi", str(point.phi))
                    molecule.SetProp("quick_qm_polar_phi_unit", phi_unit)
                    molecule.SetProp("quick_qm_polar_theta", str(point.theta))
                    molecule.SetProp("quick_qm_polar_theta_unit", theta_unit)
                molecule.SetProp("quick_qm_x_offset_angstrom", str(x))
                molecule.SetProp("quick_qm_z_offset_angstrom", str(z))
                if frame is not None:
                    molecule.SetProp("quick_qm_reference_method", frame.method)
                    molecule.SetProp("quick_qm_x_axis", json.dumps(frame.x.tolist()))
                    molecule.SetProp("quick_qm_z_axis", json.dumps(frame.z.tolist()))
                molecule.SetProp("quick_qm_rigid_source", json.dumps(rigid.fields, separators=(",", ":")))
                molecule.SetProp("quick_qm_flexible_source", json.dumps(flexible.fields, separators=(",", ":")))
                if output_format == "pair-files" and current_axis_pair != axis_pair:
                    if writer is not None:
                        writer.flush()
                        writer.close()
                    writer = Chem.SDWriter(str(staging / f"{axis_pair_id}.sdf"))
                    current_axis_pair = axis_pair
                writer.write(molecule)
                encode = lambda value: json.dumps(value, separators=(",", ":"), allow_nan=False)
                table.writerow({
                    "pair_id": pair_id, "axis_pair_id": axis_pair_id, "pose_index": pose_index,
                    "sdf_file": os.path.relpath(sdf),
                    "sdf_record_index": pose_index - 1 if output_format == "pair-files" else index - 1,
                    "rigid_csv_line": rigid.line, "flexible_csv_line": flexible.line,
                    "rigid_source": encode(rigid.fields), "flexible_source": encode(flexible.fields),
                    "distance_angstrom": distance, "rigid_center": encode(rigid.center.tolist()),
                    "flexible_original_center": encode(flexible.center.tolist()),
                    "flexible_aligned_center": encode(target.tolist()), "rigid_axis": encode(rigid.axis.tolist()),
                    "flexible_original_axis": encode(flexible.axis.tolist()),
                    "flexible_aligned_axis": encode((matrix @ flexible.axis).tolist()),
                    "rotation": encode(matrix.tolist()), "translation": encode(translation.tolist()),
                    "x_offset_angstrom": x, "z_offset_angstrom": z,
                    "x_axis": encode(frame.x.tolist()) if frame is not None else "",
                    "z_axis": encode(frame.z.tolist()) if frame is not None else "",
                    "reference_method": frame.method if frame is not None else "",
                    "reference_atom_index": frame.atom_index if frame is not None and frame.atom_index is not None else "",
                    "reference_point": encode(frame.point.tolist()) if frame is not None else "",
                    "flexible_orientation": orientation,
                    "coordinate_system": "polar" if point is not None else "grid",
                    "polar_radius_angstrom": point.radius if point is not None else "",
                    "polar_phi": point.phi if point is not None else "",
                    "polar_phi_unit": phi_unit,
                    "polar_theta": point.theta if point is not None else "",
                    "polar_theta_unit": theta_unit,
                })
            if writer is not None:
                writer.flush()
                writer.close()
                writer = None
        staged_files = ([staging / path.name for path in destinations[:-1]]
                        if output_format == "pair-files" else [staging / "single.sdf"])
        staged_files.append(staging / "metadata.csv")
        installed, backups = [], []
        created_directory = False
        try:
            if output_format == "pair-files":
                if not overwrite or not output.exists():
                    output.mkdir()
                    created_directory = True
                elif output.is_symlink() or not output.is_dir():
                    raise ValueError("pair-files output must be a directory, not a file or symlink")
            for file_index, (source, destination) in enumerate(zip(staged_files, destinations)):
                if destination.is_dir():
                    raise ValueError(f"output file is a directory: {destination}")
                if overwrite and os.path.lexists(destination):
                    backup = staging / f"backup-{file_index}"
                    os.replace(destination, backup)
                    backups.append((backup, destination))
                if overwrite:
                    os.replace(source, destination)
                else:
                    os.link(source, destination)
                installed.append(destination)
        except (OSError, ValueError, RuntimeError):
            for path in reversed(installed):
                try:
                    path.unlink()
                except OSError as cleanup_error:
                    logger.error("could not remove incomplete output %s: %s", path, cleanup_error)
            for backup, destination in reversed(backups):
                try:
                    os.replace(backup, destination)
                except OSError as cleanup_error:
                    preserve_staging = True
                    logger.error("could not restore %s: %s; backup retained at %s",
                                 destination, cleanup_error, backup)
            if created_directory:
                try:
                    output.rmdir()
                except OSError as cleanup_error:
                    logger.error("could not remove output directory %s: %s", output, cleanup_error)
            raise
        return index
    finally:
        if writer is not None:
            writer.close()
        if not preserve_staging:
            try:
                shutil.rmtree(staging)
            except OSError as cleanup_error:
                logger.error("could not clean up temporary output %s: %s", staging, cleanup_error)


if __name__ == "__main__":
    sys.exit(main())
