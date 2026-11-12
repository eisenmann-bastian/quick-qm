#!/usr/bin/env python3
"""Derive candidate alignment axes from selector CSV results and 3D SDFs."""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Geometry import Point3D

from core.logging import configure_logging, get_logger
from selection.axis_config import AXIS_NAMES, DEFAULT_CONFIG, AxisConfig
from selection.axis_calculator import AlignmentCalculator
from selection.definitions import carboxylate_sites
from selection.matching import MoleculeSelector

logger = get_logger("selection.alignment")

INPUT_COLUMNS = (
    "csv_line", "sdf_file", "record_index", "molecule_name", "selector",
    "smarts", "atom_index_base", "match_index", "atom_indices", "site_atom_indices",
)
GROUP_COLUMNS = ("group_type", "group_atom_indices", "group_charge")
OUTPUT_COLUMNS = INPUT_COLUMNS + ("axis_type", "axis_name", "center_atom_index", "center", "axis") + GROUP_COLUMNS
SDF_AXES_PROPERTY = "quick_qm_alignment_axes"
SDF_SOURCE_PROPERTY = "quick_qm_source"


class AlignmentCsv:
    """Read selector CSV rows and preserve their provenance for each axis."""

    def __init__(self, axis_config=None):
        self.calculator = AlignmentCalculator()
        self.axis_config = axis_config if axis_config is not None else AxisConfig.load()
        self.molecules = {}
        self.records = {}

    def read(self, path, *, include_unmatched=False):
        with Path(path).open(encoding="utf-8-sig", newline="") as source:
            return self.read_stream(source, include_unmatched=include_unmatched)

    def read_stream(self, source, *, include_unmatched=False):
        """Process selection CSV from a stream without an intermediate file."""
        results = []
        self.molecules = {}
        self.records = {}
        reader = csv.DictReader(source)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("CSV must have unique column headers")
        missing = set(INPUT_COLUMNS) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            try:
                if None in row or any(row[key] is None for key in INPUT_COLUMNS):
                    raise ValueError("incorrect number of CSV columns")
                if row["atom_index_base"] != "0":
                    raise ValueError("only zero-based atom indices are supported")
                if not row["match_index"]:
                    if row["atom_indices"] or row["site_atom_indices"]:
                        raise ValueError("empty match_index with nonempty atom indices")
                    if not include_unmatched:
                        continue
                elif int(row["match_index"]) < 0:
                    raise ValueError("match_index must be non-negative")
                record_index = int(row["record_index"])
                if record_index < 0:
                    raise ValueError("record_index must be non-negative")
                sdf_path = Path(row["sdf_file"])
                if not row["sdf_file"]:
                    raise ValueError("sdf_file must not be empty")
                if sdf_path not in self.molecules:
                    self.molecules[sdf_path] = list(MoleculeSelector._read_sdf(sdf_path))
                molecules = self.molecules[sdf_path]
                if record_index >= len(molecules):
                    raise ValueError("record_index is outside the SDF")
                key = (sdf_path.resolve(), record_index)
                self.records.setdefault(key, molecules[record_index])
                if not row["match_index"]:
                    continue
                groups = carboxylate_sites(molecules[record_index])
                direction_counts = {}
                for kind, center_index, center, axis in self.calculator.calculate(molecules[record_index], row):
                    name = AXIS_NAMES[kind]
                    key = (name, center_index)
                    direction = direction_counts.get(key, 0) + 1
                    direction_counts[key] = direction
                    if not self.axis_config.enabled(name, direction):
                        continue
                    output = {key: row[key] for key in INPUT_COLUMNS}
                    output["sdf_file"] = os.path.relpath(sdf_path)
                    output.update(axis_type=kind, axis_name=f"{name}_{direction}", center_atom_index=center_index,
                                  center=json.dumps(center, separators=(",", ":"), allow_nan=False),
                                  axis=json.dumps(axis, separators=(",", ":"), allow_nan=False))
                    group = groups.get(center_index)
                    output.update(group_type="carboxylate" if group else "",
                                  group_atom_indices=json.dumps(group, separators=(",", ":")) if group else "",
                                  group_charge=-1 if group else "")
                    results.append(output)
            except (ValueError, OSError, RuntimeError) as exc:
                raise ValueError(f"CSV line {reader.line_num}: {exc}") from exc
        return results

    @staticmethod
    def write(results, destination):
        writer = csv.DictWriter(destination, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)


class AlignmentSdf:
    """Create one annotated SDF record for each referenced source record."""

    @staticmethod
    def build(records, results):
        grouped = {key: [] for key in records}
        for row in results:
            key = (Path(row["sdf_file"]).resolve(), int(row["record_index"]))
            center = json.loads(row["center"])
            axis = json.loads(row["axis"])
            existing = next((entry for entry in grouped[key]
                             if all(abs(a - b) <= 1e-8 for a, b in zip(entry["center"], center))
                             and all(abs(a - b) <= 1e-8 for a, b in zip(entry["axis"], axis))), None)
            provenance = {name: row[name] for name in INPUT_COLUMNS}
            provenance["axis_type"] = row["axis_type"]
            provenance["axis_name"] = row["axis_name"]
            provenance["center_atom_index"] = row["center_atom_index"]
            if row["group_type"]:
                provenance.update(group_type=row["group_type"],
                                  group_atom_indices=json.loads(row["group_atom_indices"]),
                                  group_charge=row["group_charge"])
            if existing is None:
                grouped[key].append({"center": center, "axis": axis, "sources": [provenance]})
            else:
                existing["sources"].append(provenance)

        output = []
        for (path, record_index), source in records.items():
            for name in (SDF_AXES_PROPERTY, SDF_SOURCE_PROPERTY):
                if source.HasProp(name):
                    raise ValueError(f"SDF {path}, record {record_index} already has property {name}")
            mol = Chem.RWMol(source)
            entries = grouped[(path, record_index)]
            for entry in entries:
                atom_index = mol.AddAtom(Chem.Atom(0))
                coordinates = [center + direction for center, direction in zip(entry["center"], entry["axis"])]
                mol.GetConformer().SetAtomPosition(atom_index, Point3D(*coordinates))
                entry["pseudo_atom_index"] = atom_index
            if not entries:
                logger.warning("no supported alignment axes or all axes disabled: SDF %s, record %s; writing unchanged molecule",
                               os.path.relpath(path), record_index)
            mol.SetProp(SDF_SOURCE_PROPERTY, json.dumps(
                {"sdf_file": os.path.relpath(path), "record_index": record_index}, separators=(",", ":")
            ))
            mol.SetProp(SDF_AXES_PROPERTY, json.dumps(entries, separators=(",", ":"), allow_nan=False))
            output.append(mol.GetMol())
        return output

    @staticmethod
    def write(molecules, destination):
        writer = Chem.SDWriter(destination)
        try:
            for molecule in molecules:
                writer.write(molecule)
            writer.flush()
        finally:
            writer.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path, help="CSV output from the selector")
    parser.add_argument("--format", choices=("csv", "sdf"), default="csv", help="output format (default: csv)")
    parser.add_argument("--axes-config", type=Path, default=DEFAULT_CONFIG,
                        help="YAML axis configuration (default: bundled axes.config)")
    parser.add_argument("-o", "--output", type=Path, help="write to a new file instead of stdout")
    args = parser.parse_args(argv)
    if args.format == "sdf" and args.output is None:
        parser.error("--format sdf requires --output")
    configure_logging()
    try:
        processor = AlignmentCsv(AxisConfig.load(args.axes_config))
        results = processor.read(args.input_csv, include_unmatched=args.format == "sdf")
        if args.format == "sdf":
            molecules = AlignmentSdf.build(processor.records, results)
            with args.output.open("x", encoding="utf-8", newline="") as destination:
                AlignmentSdf.write(molecules, destination)
        elif args.output is None:
            processor.write(results, sys.stdout)
        else:
            with args.output.open("x", encoding="utf-8", newline="") as destination:
                processor.write(results, destination)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
