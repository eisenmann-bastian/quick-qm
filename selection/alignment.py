#!/usr/bin/env python3
"""Derive candidate alignment axes from selector CSV results and 3D SDFs."""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Geometry import Point3D

from core.logging import configure_logging, get_logger
from selection.axis_config import AXIS_NAMES, DEFAULT_CONFIG, AxisConfig
from selection.definitions import carboxylate_sites
from selection.matching import MoleculeSelector

logger = get_logger("selection.alignment")

HYDROGEN = 1
CARBON = 6
NITROGEN = 7
OXYGEN = 8
FLUORINE = 9
SULFUR = 16
CHLORINE = 17
SELENIUM = 34
BROMINE = 35
IODINE = 53
TELLURIUM = 52

HALOGENS = (FLUORINE, CHLORINE, BROMINE, IODINE)
CHALCOGENS = (OXYGEN, SULFUR, SELENIUM)
SINGLE_BOND_CHALCOGENS = (*CHALCOGENS, TELLURIUM)
HYDROGEN_DONOR_ELEMENTS = (NITROGEN, *CHALCOGENS)

INPUT_COLUMNS = (
    "csv_line", "sdf_file", "record_index", "molecule_name", "selector",
    "smarts", "atom_index_base", "match_index", "atom_indices", "site_atom_indices",
)
GROUP_COLUMNS = ("group_type", "group_atom_indices", "group_charge")
OUTPUT_COLUMNS = INPUT_COLUMNS + ("axis_type", "axis_name", "center_atom_index", "center", "axis") + GROUP_COLUMNS
SDF_AXES_PROPERTY = "quick_qm_alignment_axes"
SDF_SOURCE_PROPERTY = "quick_qm_source"


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _unit(vector):
    length = math.sqrt(sum(x * x for x in vector))
    if not math.isfinite(length) or length < 1e-8:
        raise ValueError("degenerate or non-finite geometry cannot define an axis")
    return tuple(x / length for x in vector)


def _position(conformer, index):
    p = conformer.GetAtomPosition(index)
    if not all(math.isfinite(x) for x in (p.x, p.y, p.z)):
        raise ValueError("non-finite atomic coordinates")
    return (p.x, p.y, p.z)


class AlignmentCalculator:
    """Derive outward unit vectors without modifying the input conformer."""

    def calculate(self, mol, row):
        atom_indices = self._indices(row["atom_indices"], mol, "atom_indices")
        sites = self._indices(row["site_atom_indices"], mol, "site_atom_indices")
        if not set(sites).issubset(atom_indices):
            raise ValueError("site_atom_indices must be contained in atom_indices")
        if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
            raise ValueError("SDF record must have a 3D conformer")
        conformer = mol.GetConformer()
        if row["selector"] == "aromatic ring":
            return self._ring(mol, conformer, sites)

        axes = []
        groups = carboxylate_sites(mol)
        for index in sites:
            found = self._site(mol, conformer, index, groups.get(index))
            if not found:
                logger.warning(
                    "not yet supported: selector %r, SDF %s, record %s, atom %s",
                    row["selector"] or row["smarts"], row["sdf_file"],
                    row["record_index"], index,
                )
            axes.extend(found)
        return axes

    @staticmethod
    def _indices(value, mol, label):
        try:
            indices = json.loads(value)
        except (json.JSONDecodeError, TypeError) as exc:
            raise ValueError(f"{label} must be a JSON array of atom indices") from exc
        if (not isinstance(indices, list) or not indices
                or any(type(i) is not int or i < 0 or i >= mol.GetNumAtoms() for i in indices)
                or len(set(indices)) != len(indices)):
            raise ValueError(f"{label} must contain distinct, valid zero-based atom indices")
        return indices

    @staticmethod
    def _ring(mol, conformer, sites):
        rings = mol.GetRingInfo().AtomRings()
        if not any(set(ring) == set(sites) and all(
            mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)]).GetIsAromatic()
            for i in range(len(ring))
        ) for ring in rings):
            raise ValueError("aromatic ring site indices do not identify an aromatic ring")
        # Use the SDF's ring ordering, not the order of the CSV cell.
        ring = next(r for r in rings if set(r) == set(sites) and all(
            mol.GetBondBetweenAtoms(r[i], r[(i + 1) % len(r)]).GetIsAromatic()
            for i in range(len(r))
        ))
        points = [_position(conformer, i) for i in ring]
        center = tuple(sum(p[j] for p in points) / len(points) for j in range(3))
        normal = [0.0, 0.0, 0.0]
        for i, point in enumerate(points):
            area = _cross(_sub(point, center), _sub(points[(i + 1) % len(points)], center))
            normal = [a + b for a, b in zip(normal, area)]
        normal = _unit(normal)
        return [
            ("aromatic_normal", None, center, normal),
            ("aromatic_normal", None, center, tuple(-x for x in normal)),
        ]

    def _site(self, mol, conformer, index, carboxylate=None):
        atom = mol.GetAtomWithIdx(index)
        element = atom.GetAtomicNum()
        center = _position(conformer, index)
        axes = []

        def extension(neighbor, axis_type, center_index=index):
            origin = _position(conformer, center_index)
            direction = _unit(_sub(origin, _position(conformer, neighbor.GetIdx())))
            axes.append((axis_type, center_index, origin, direction))

        if carboxylate is not None:
            carbon_index, oxygen1, oxygen2 = carboxylate
            other_oxygen = oxygen2 if index == oxygen1 else oxygen1
            carbon = _position(conformer, carbon_index)
            direction = _unit(_sub(center, carbon))
            relative = _sub(_position(conformer, other_oxygen), carbon)
            perpendicular = _unit(_cross(_cross(direction, relative), direction))
            axes.append(("carboxylate_sigma", index, center, direction))
            for sign in (-1, 1):
                vector = tuple(0.5 * d + sign * math.sqrt(3) / 2 * p
                               for d, p in zip(direction, perpendicular))
                axes.append(("carboxylate_lone_pair", index, center, _unit(vector)))
            return axes

        bonds = list(atom.GetBonds())
        if element == HYDROGEN and len(bonds) == 1 and bonds[0].GetBondType() == Chem.BondType.SINGLE:
            extension(bonds[0].GetOtherAtom(atom), "hydrogen_sigma")
        elif element in HALOGENS and len(bonds) == 1 and bonds[0].GetBondType() == Chem.BondType.SINGLE:
            extension(bonds[0].GetOtherAtom(atom), "halogen_sigma")
        elif element in SINGLE_BOND_CHALCOGENS:
            for bond in bonds:
                neighbor = bond.GetOtherAtom(atom)
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() > HYDROGEN:
                    extension(neighbor, "chalcogen_sigma")
                    continue
                if element not in CHALCOGENS or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                extension(neighbor, "chalcogen_sigma")
                if neighbor.GetAtomicNum() == CARBON:
                    substituents = [n for n in neighbor.GetNeighbors() if n.GetIdx() != index]
                    if len(substituents) == 2:
                        direction = axes[-1][3]
                        carbon = _position(conformer, neighbor.GetIdx())
                        sides = []
                        for substituent in substituents:
                            relative = _sub(_position(conformer, substituent.GetIdx()), carbon)
                            sides.append(_unit(_cross(_cross(direction, relative), direction)))
                        perpendicular = _unit(_sub(sides[0], sides[1]))
                        for sign in (-1, 1):
                            vector = tuple(0.5 * d + sign * math.sqrt(3) / 2 * p
                                           for d, p in zip(direction, perpendicular))
                            axes.append(("chalcogen_lone_pair", index, center, _unit(vector)))
                    elif neighbor.GetTotalNumHs() > 0:
                        logger.warning(
                            "not yet supported: carbonyl lone-pair plane at atom %s needs explicit substituent coordinates",
                            index,
                        )
            neighbor_elements = sorted(bond.GetOtherAtom(atom).GetAtomicNum() for bond in bonds)
            sulfur_geometry = None
            if (element == SULFUR and len(bonds) == 2
                    and all(bond.GetBondType() == Chem.BondType.SINGLE for bond in bonds)):
                if neighbor_elements == [CARBON, CARBON] and atom.GetTotalNumHs(includeNeighbors=True) == 0:
                    sulfur_geometry = "thioether"
                elif (atom.GetFormalCharge() == 0 and neighbor_elements == [HYDROGEN, CARBON]
                      and atom.GetTotalNumHs(includeNeighbors=True) == 1):
                    sulfur_geometry = "thiol"
            if sulfur_geometry:
                if sulfur_geometry == "thiol":
                    hydrogen_bond = next(bond for bond in bonds
                                         if bond.GetOtherAtom(atom).GetAtomicNum() == HYDROGEN)
                    extension(hydrogen_bond.GetOtherAtom(atom), "chalcogen_sigma")
                vectors = [
                    _unit(_sub(_position(conformer, bond.GetOtherAtomIdx(index)), center))
                    for bond in bonds
                ]
                normal = _unit(_cross(vectors[0], vectors[1]))
                bisector = _unit(tuple(-sum(vector[i] for vector in vectors) for i in range(3)))
                axes.extend([
                    (f"{sulfur_geometry}_lone_pair_normal", index, center, normal),
                    (f"{sulfur_geometry}_lone_pair_normal", index, center, tuple(-x for x in normal)),
                    (f"{sulfur_geometry}_outward_bisector", index, center, bisector),
                ])
            elif (element == SULFUR and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 1
                  and len(bonds) == 1 and neighbor_elements == [CARBON]
                  and bonds[0].GetBondType() == Chem.BondType.SINGLE):
                logger.warning(
                    "not yet supported: thiol plane and S-H sigma extension at atom %s need explicit hydrogen coordinates",
                    index,
                )
        elif element == NITROGEN and atom.GetFormalCharge() == 0:
            triple = [b for b in bonds if b.GetBondType() == Chem.BondType.TRIPLE
                      and b.GetOtherAtom(atom).GetAtomicNum() == CARBON]
            if len(triple) == 1 and len(bonds) == 1:
                extension(triple[0].GetOtherAtom(atom), "nitrogen_triple_sigma")
            elif atom.IsInRing() and atom.GetTotalNumHs(includeNeighbors=True) == 0 and len(bonds) == 2 and (
                atom.GetIsAromatic() or any(b.GetBondType() == Chem.BondType.DOUBLE for b in bonds)
            ):
                vectors = [_unit(_sub(_position(conformer, b.GetOtherAtomIdx(index)), center)) for b in bonds]
                direction = _unit(tuple(-sum(v[i] for v in vectors) for i in range(3)))
                axes.append(("nitrogen_ring_lone_pair", index, center, direction))
            elif len(bonds) == 2:
                double_c = [b for b in bonds if b.GetBondType() == Chem.BondType.DOUBLE
                            and b.GetOtherAtom(atom).GetAtomicNum() == CARBON]
                single_h = [b for b in bonds if b.GetBondType() == Chem.BondType.SINGLE
                            and b.GetOtherAtom(atom).GetAtomicNum() == HYDROGEN]
                if len(double_c) == len(single_h) == 1:
                    vectors = [_unit(_sub(_position(conformer, b.GetOtherAtomIdx(index)), center)) for b in bonds]
                    direction = _unit(tuple(-sum(v[i] for v in vectors) for i in range(3)))
                    axes.append(("nitrogen_imine_lone_pair", index, center, direction))

        if element in HYDROGEN_DONOR_ELEMENTS:
            for bond in bonds:
                hydrogen = bond.GetOtherAtom(atom)
                if hydrogen.GetAtomicNum() == HYDROGEN and bond.GetBondType() == Chem.BondType.SINGLE:
                    extension(atom, "hydrogen_sigma", hydrogen.GetIdx())
            if atom.GetTotalNumHs() > 0:
                logger.warning(
                    "not yet supported: hydrogen axes at atom %s need explicit hydrogen coordinates",
                    index,
                )
        return axes


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
