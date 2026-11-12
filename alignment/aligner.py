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

    def align(self, rigid, flexible, rigid_mol, flexible_mol, distance):
        if distance is None:
            distance = self.radius(rigid, rigid_mol) + self.radius(flexible, flexible_mol)
        matrix = self.rotation(flexible.axis, -rigid.axis)
        target = rigid.center + distance * rigid.axis
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
    def from_options(cls, settings, args):
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
                return cls(mode)
            value = cls.number(distance, "distance")
            if value < 0:
                raise ValueError("distance must be non-negative")
            return cls(mode, distance=float(value))
        if distance is not None:
            raise ValueError("distance cannot be used with distance-scan")
        if scan["scan_start"] is None or scan["scan_end"] is None:
            raise ValueError("distance-scan requires scan_start and scan_end")
        start = cls.number(scan["scan_start"], "scan_start")
        end = cls.number(scan["scan_end"], "scan_end")
        if start < 0 or end <= start:
            raise ValueError("scan range must have 0 <= scan_start < scan_end")
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
    unknown = set(data) - {"rigid", "flexible", "output", "output_format", "distance",
                           "alignment_mode", "scan_start", "scan_end", "include_start",
                           "include_end", "scan_step", "scan_count", "overwrite_output"}
    if unknown:
        raise ValueError(f"unknown configuration keys: {', '.join(sorted(map(str, unknown)))}")
    return data


def triplets(values):
    return [dict(zip(FILTER_FIELDS, group)) for group in values]


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
    parser.add_argument("--distance", type=float, help="fixed center separation in Å; default: sum of van der Waals radii")
    parser.add_argument("--alignment-mode", choices=("default", "distance-scan"))
    parser.add_argument("--scan-start", help="scan lower boundary in Å")
    parser.add_argument("--scan-end", help="scan upper boundary in Å")
    parser.add_argument("--include-start", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--include-end", action=argparse.BooleanOptionalAction, default=None)
    sampling = parser.add_mutually_exclusive_group()
    sampling.add_argument("--scan-step", help="fixed scan step in Å")
    sampling.add_argument("--scan-count", type=int, help="number of interior poses per axis pair; boundaries are additional")
    args = parser.parse_args(argv)
    configure_logging()
    try:
        settings = load_config(args.config)
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
        plan = DistancePlan.from_options(settings, args)
        overwrite = args.overwrite_output if args.overwrite_output is not None else settings.get("overwrite_output", False)
        if type(overwrite) is not bool:
            raise ValueError("overwrite_output must be a boolean")
        if not overwrite and (os.path.lexists(output) or os.path.lexists(csv_output)):
            raise ValueError(f"output already exists: {output} or {csv_output}")
        reader = AxisInput()
        rigid_rows, flexible_rows = reader.read(args.input_csv, rigid, flexible)
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
                for f in flexible_rows:
                    axis_pair += 1
                    for pose_index, distance in enumerate(plan.distances(), 1):
                        mol, matrix, translation, target, separation = aligner.align(
                            r, f, reader.molecule(r), reader.molecule(f), distance
                        )
                        yield (axis_pair, pose_index, r, f, mol, matrix, translation, target, separation)

        pairs = generate_pairs()
        write_output(output, output_format, pairs, destinations, overwrite)
        logger.info("wrote %d aligned poses and %s", len(rigid_rows) * len(flexible_rows) * plan.count, csv_output)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


def write_output(output, output_format, pairs, destinations, overwrite):
    fields = ("pair_id", "axis_pair_id", "pose_index", "sdf_file", "sdf_record_index", "rigid_csv_line", "flexible_csv_line",
              "rigid_source", "flexible_source", "distance_angstrom", "rigid_center",
              "flexible_original_center", "flexible_aligned_center", "rigid_axis",
              "flexible_original_axis", "flexible_aligned_axis", "rotation", "translation")
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
            for index, (axis_pair, pose_index, rigid, flexible, molecule, matrix, translation, target, distance) in enumerate(pairs, 1):
                pair_id = f"pair-{index:06d}"
                axis_pair_id = f"axis-pair-{axis_pair:06d}"
                sdf = output / f"{axis_pair_id}.sdf" if output_format == "pair-files" else output
                molecule.SetProp("_Name", pair_id)
                molecule.SetProp("quick_qm_pair_id", pair_id)
                molecule.SetProp("quick_qm_axis_pair_id", axis_pair_id)
                molecule.SetProp("quick_qm_pose_index", str(pose_index))
                molecule.SetProp("quick_qm_distance_angstrom", str(distance))
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
            for index, (source, destination) in enumerate(zip(staged_files, destinations)):
                if destination.is_dir():
                    raise ValueError(f"output file is a directory: {destination}")
                if overwrite and os.path.lexists(destination):
                    backup = staging / f"backup-{index}"
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
