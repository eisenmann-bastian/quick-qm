#!/usr/bin/env python3
"""Align pairs of selector axes and write combined molecular structures."""

import argparse
import csv
import os
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from alignment import LOGGER_NAMESPACE
from alignment.axis_input import AxisFilter, AxisInput, AxisRow
from alignment.geometry import PairAligner
from alignment.output import write_output
from alignment.pose_generation import PoseGenerator
from alignment.reference_frame import ReferenceSelector
from alignment.scan_planning import DistancePlan, load_config, scan_plans, triplets
from alignment.torsion_optimization import TorsionOptimizer
from molecule_utilities.logging import configure_logging, get_logger

logger = get_logger("alignment.aligner", namespace=LOGGER_NAMESPACE)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path, help="alignment CSV produced by selector.py")
    parser.add_argument("--config", type=Path, help="YAML parameters (input CSV is always positional)")
    parser.add_argument("--rigid", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--flexible", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--rigid-exclude", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--flexible-exclude", nargs=3, action="append", metavar=("MOLECULE", "SELECTOR", "AXIS"))
    parser.add_argument("--output", "-o", type=Path, help="new pair directory or new multi-record SDF")
    parser.add_argument("--output-format", choices=("pair-files", "single"), help="SDF layout (default: pair-files)")
    parser.add_argument("--overwrite-output", action=argparse.BooleanOptionalAction, default=None,
                        help="replace matching output files (default: false)")
    parser.add_argument("--flexible-orientation", choices=("stable", "toward-center"),
                        help="orientation after placement (default: stable)")
    parser.add_argument("--optimize-torsion", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--torsion-score", choices=("minimum_distance", "average_distance", "mmff94_interaction"))
    parser.add_argument("--torsion-objective", choices=("maximize", "minimize"))
    parser.add_argument("--torsion-steps", type=int)
    parser.add_argument("--torsion-iterations", type=int)
    parser.add_argument("--torsion-include-vdw", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--torsion-exclude-flexible-site", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--torsion-exclude-rigid-site", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--torsion-exclude-flexible-atoms", nargs="+", type=int)
    parser.add_argument("--torsion-exclude-rigid-atoms", nargs="+", type=int)
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
    parser.add_argument("--scan-merge-tolerance", help="cross-scan position tolerance in Å (default: 0.05)")
    args = parser.parse_args(argv)
    configure_logging(namespace=LOGGER_NAMESPACE)
    try:
        settings = load_config(args.config)
        optimization = settings.get("torsion_optimization", {})
        if not isinstance(optimization, dict):
            raise ValueError("torsion_optimization must be a mapping")
        optimization = dict(optimization)
        for key, arg in (("enabled", "optimize_torsion"), ("score", "torsion_score"),
                         ("objective", "torsion_objective"), ("steps_per_iteration", "torsion_steps"),
                         ("iterations", "torsion_iterations"), ("include_vdw", "torsion_include_vdw"),
                         ("exclude_flexible_site", "torsion_exclude_flexible_site"),
                         ("exclude_rigid_site", "torsion_exclude_rigid_site"),
                         ("exclude_flexible_atoms", "torsion_exclude_flexible_atoms"),
                         ("exclude_rigid_atoms", "torsion_exclude_rigid_atoms")):
            if getattr(args, arg) is not None:
                optimization[key] = getattr(args, arg)
        optimizer = TorsionOptimizer(optimization)
        orientation = args.flexible_orientation if args.flexible_orientation is not None else settings.get("flexible_orientation", "stable")
        if orientation not in ("stable", "toward-center"):
            raise ValueError("flexible_orientation must be 'stable' or 'toward-center'")
        rigid = AxisFilter(settings.get("rigid"), triplets(args.rigid) if args.rigid is not None else None,
                           triplets(args.rigid_exclude) if args.rigid_exclude is not None else None)
        flexible = AxisFilter(settings.get("flexible"), triplets(args.flexible) if args.flexible is not None else None,
                              triplets(args.flexible_exclude) if args.flexible_exclude is not None else None)
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
        scans = scan_plans(settings, args)
        merge_tolerance = args.scan_merge_tolerance if args.scan_merge_tolerance is not None else settings.get("scan_merge_tolerance", 0.05)
        merge_tolerance = float(DistancePlan.number(merge_tolerance, "scan_merge_tolerance"))
        if merge_tolerance <= 0:
            raise ValueError("scan_merge_tolerance must be positive")
        reference_selector = ReferenceSelector(settings.get("reference_axis"))
        overwrite = args.overwrite_output if args.overwrite_output is not None else settings.get("overwrite_output", False)
        if type(overwrite) is not bool:
            raise ValueError("overwrite_output must be a boolean")
        if not overwrite and (os.path.lexists(output) or os.path.lexists(csv_output)):
            raise ValueError(f"output already exists: {output} or {csv_output}")
        reader = AxisInput()
        rigid_rows, flexible_rows = reader.read(args.input_csv, rigid, flexible)
        needs_frame = any(scan.needs_frame for scan in scans)
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
        pairs = PoseGenerator(reader, rigid_rows, flexible_rows, frames, scans, merge_tolerance,
                              orientation, optimizer).generate()
        count = write_output(output, output_format, pairs, destinations, overwrite)
        logger.info("wrote %d aligned poses and %s", count, csv_output)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
