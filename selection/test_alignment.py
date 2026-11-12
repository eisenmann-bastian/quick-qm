import csv
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from selection.alignment import (
    AlignmentCalculator, AlignmentCsv, AlignmentSdf, INPUT_COLUMNS, SDF_AXES_PROPERTY, SDF_SOURCE_PROPERTY,
)
from selection.matching import MoleculeSelector
from selection.output import CsvResultWriter


ROOT = Path(__file__).resolve().parent.parent


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.sdf = self.directory / "molecule.sdf"
        self.csv = self.directory / "matches.csv"

    def molecule(self, smiles, explicit_hydrogens=True):
        mol = Chem.MolFromSmiles(smiles)
        if explicit_hydrogens:
            mol = Chem.AddHs(mol)
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=17), 0)
        return mol

    def prepare(self, mol, *, name=None, smarts=None):
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(mol)
        writer.close()
        selected = MoleculeSelector().select(self.sdf, name=name, smarts=smarts)
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            CsvResultWriter().write(selected, destination, sdf_file=self.sdf)
        return AlignmentCsv().read(self.csv)

    def assert_unit(self, vector):
        self.assertAlmostEqual(dot(vector, vector), 1.0, places=7)

    def test_carbonyl_sigma_and_lone_pair_angles(self):
        for smiles in ("CC(=O)C", "CC(=S)C", "CC(=[Se])C"):
            with self.subTest(smiles=smiles):
                mol = self.molecule(smiles)
                rows = self.prepare(mol, smarts="[#8,#16,#34]")
                self.assertEqual([r["axis_type"] for r in rows],
                                 ["chalcogen_sigma", "chalcogen_lone_pair", "chalcogen_lone_pair"])
                axes = [json.loads(r["axis"]) for r in rows]
                for axis in axes:
                    self.assert_unit(axis)
                self.assertAlmostEqual(dot(axes[0], axes[1]), 0.5, places=7)
                self.assertAlmostEqual(dot(axes[0], axes[2]), 0.5, places=7)
                self.assertAlmostEqual(dot(axes[1], axes[2]), -0.5, places=7)
                self.assertEqual(rows[0]["center_atom_index"], 2)
                conformer = mol.GetConformer()
                original = conformer.GetAtomPosition(2)
                for actual, expected in zip(json.loads(rows[0]["center"]), original):
                    self.assertAlmostEqual(actual, expected, places=3)

    def test_carboxylate_geometry_is_resonance_invariant(self):
        mol = self.molecule("CC(=O)[O-]")
        first = self.prepare(mol, name="carboxylate")
        original = next(MoleculeSelector._read_sdf(self.sdf))
        conformer = original.GetConformer()
        carbon = conformer.GetAtomPosition(1)
        p1, p2 = conformer.GetAtomPosition(2), conformer.GetAtomPosition(3)
        v1 = (p1.x - carbon.x, p1.y - carbon.y, p1.z - carbon.z)
        v2 = (p2.x - carbon.x, p2.y - carbon.y, p2.z - carbon.z)
        plane_normal = (v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2],
                        v1[0]*v2[1]-v1[1]*v2[0])
        self.assertEqual(len(first), 6)
        for oxygen in (2, 3):
            rows = [r for r in first if r["center_atom_index"] == oxygen]
            self.assertEqual([r["axis_type"] for r in rows],
                             ["carboxylate_sigma", "carboxylate_lone_pair", "carboxylate_lone_pair"])
            axes = [json.loads(row["axis"]) for row in rows]
            for row, axis in zip(rows, axes):
                self.assert_unit(axis)
                self.assertAlmostEqual(dot(axis, plane_normal), 0, places=7)
                self.assertEqual(row["group_type"], "carboxylate")
                self.assertEqual(json.loads(row["group_atom_indices"]), [1, 2, 3])
                self.assertEqual(row["group_charge"], -1)
            self.assertAlmostEqual(dot(axes[0], axes[1]), 0.5, places=7)
            self.assertAlmostEqual(dot(axes[0], axes[2]), 0.5, places=7)
            self.assertAlmostEqual(dot(axes[1], axes[2]), -0.5, places=7)
        swapped = Chem.RWMol(mol)
        swapped.GetBondBetweenAtoms(1, 2).SetBondType(Chem.BondType.SINGLE)
        swapped.GetBondBetweenAtoms(1, 3).SetBondType(Chem.BondType.DOUBLE)
        swapped.GetAtomWithIdx(2).SetFormalCharge(-1)
        swapped.GetAtomWithIdx(3).SetFormalCharge(0)
        swapped = swapped.GetMol()
        Chem.SanitizeMol(swapped)
        second = self.prepare(swapped, name="carboxylate")
        def keyed(rows):
            return {(r["center_atom_index"], r["axis_type"], r["axis"])
                    for r in rows}
        self.assertEqual(keyed(first), keyed(second))

    def test_carboxylate_axes_apply_to_legacy_selectors(self):
        mol = self.molecule("CC(=O)[O-]")
        all_rows = self.prepare(mol, name="carboxylate")
        carbonyl = self.prepare(mol, name="carbonyl")
        charge = self.prepare(mol, name="negative charge")
        self.assertEqual(len(carbonyl), 3)
        self.assertEqual(len(charge), 3)
        for subset in (carbonyl, charge):
            self.assertEqual({row["axis"] for row in subset},
                             {row["axis"] for row in all_rows
                              if row["center_atom_index"] == subset[0]["center_atom_index"]})
            self.assertTrue(all(row["group_charge"] == -1 for row in subset))

    def test_formate_carboxylate_without_explicit_hydrogen(self):
        mol = self.molecule("O=C[O-]", explicit_hydrogens=False)
        rows = self.prepare(mol, name="carboxylate")
        self.assertEqual(len(rows), 6)
        self.assertEqual({row["center_atom_index"] for row in rows}, {0, 2})

    def test_carboxylate_sdf_metadata_and_original_bond_orders(self):
        self.prepare(self.molecule("CC(=O)[O-]"), name="carboxylate")
        processor = AlignmentCsv()
        rows = processor.read(self.csv, include_unmatched=True)
        annotated = AlignmentSdf.build(processor.records, rows)[0]
        original = next(MoleculeSelector._read_sdf(self.sdf))
        self.assertEqual(Chem.MolToSmiles(annotated, isomericSmiles=True).count("[O-]"), 1)
        for bond in original.GetBonds():
            other = annotated.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            self.assertEqual(bond.GetBondType(), other.GetBondType())
        axes = json.loads(annotated.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 6)
        self.assertEqual(annotated.GetNumAtoms(), original.GetNumAtoms() + 6)
        for entry in axes:
            source = entry["sources"][0]
            self.assertEqual(source["group_type"], "carboxylate")
            self.assertEqual(source["group_atom_indices"], [1, 2, 3])
            self.assertEqual(source["group_charge"], -1)
        output = self.directory / "carboxylate.sdf"
        with output.open("w", encoding="utf-8", newline="") as destination:
            AlignmentSdf.write([annotated], destination)
        restored = next(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(json.loads(restored.GetProp(SDF_AXES_PROPERTY))), 6)
        for index in range(original.GetNumAtoms()):
            self.assertEqual(restored.GetAtomWithIdx(index).GetFormalCharge(),
                             original.GetAtomWithIdx(index).GetFormalCharge())
        for bond in original.GetBonds():
            other = restored.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            self.assertEqual(bond.GetBondType(), other.GetBondType())

    def test_reference_batch_carboxylate_pipeline(self):
        selected = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "selector.py"),
             "--csv", str(ROOT / "selection" / "reference_molecules.csv")],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(selected.returncode, 0, selected.stderr)
        self.csv.write_text(selected.stdout, encoding="utf-8")
        output = self.directory / "reference.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        records = list(MoleculeSelector._read_sdf(output))
        acetate = next(m for m in records if m.GetProp("_Name") == "175")
        axes = json.loads(acetate.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 6)
        self.assertEqual({source["center_atom_index"] for axis in axes for source in axis["sources"]}, {0, 1})
        self.assertTrue(all(any(s["selector"] == "carboxylate" for s in axis["sources"]) for axis in axes))
        self.assertTrue(all(s["group_charge"] == -1 for axis in axes for s in axis["sources"]))

    def test_combined_cli_path_only_batch_all_formats(self):
        sdf = ROOT / "selection" / "reference_molecules" / "175.sdf"
        self.csv.write_text(f"sdf\n{sdf}\n", encoding="utf-8")
        command = [sys.executable, str(ROOT / "selection" / "selector.py"), "--csv", str(self.csv)]
        selection = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(selection.returncode, 0, selection.stderr)
        rows = list(csv.DictReader(io.StringIO(selection.stdout)))
        self.assertIn("carboxylate", {row["selector"] for row in rows})
        carboxylate = next(row for row in rows if row["selector"] == "carboxylate")
        self.assertEqual(set(json.loads(carboxylate["site_atom_indices"])), {0, 1})
        self.assertFalse(os.path.isabs(carboxylate["sdf_file"]))
        alignment = subprocess.run(command + ["--format", "alignment"],
                                   cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(alignment.returncode, 0, alignment.stderr)
        axes = list(csv.DictReader(io.StringIO(alignment.stdout)))
        carboxylate_axes = [row for row in axes if row["selector"] == "carboxylate"]
        self.assertEqual(len(carboxylate_axes), 6)
        self.assertEqual({row["center_atom_index"] for row in carboxylate_axes}, {"0", "1"})
        self.assertTrue(all(row["group_charge"] == "-1" for row in carboxylate_axes))
        output = self.directory / "combined.sdf"
        sdf_result = subprocess.run(command + ["--format", "sdf", "-o", str(output)],
                                    cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(sdf_result.returncode, 0, sdf_result.stderr)
        annotated = list(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(annotated), 1)
        metadata = json.loads(annotated[0].GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(metadata), 6)
        self.assertTrue(all(any(s["selector"] == "carboxylate" for s in axis["sources"])
                            for axis in metadata))
        again = subprocess.run(command + ["--format", "sdf", "-o", str(output)],
                               cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertEqual(len(list(MoleculeSelector._read_sdf(output))), 1)

    def test_combined_cli_path_only_multirecord_and_unmatched(self):
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(self.molecule("CCl"))
        writer.write(self.molecule("CC"))
        writer.close()
        self.csv.write_text("sdf\nmolecule.sdf\n", encoding="utf-8")
        output = self.directory / "multiple.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "selector.py"), "--csv", str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        records = list(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(records), 2)
        self.assertEqual(len(json.loads(records[0].GetProp(SDF_AXES_PROPERTY))), 1)
        self.assertEqual(json.loads(records[1].GetProp(SDF_AXES_PROPERTY)), [])

    def test_combined_cli_rejects_bad_batch_without_output(self):
        self.prepare(self.molecule("CS"), name="thiol")
        self.csv.write_text("sdf\nmolecule.sdf\nmissing.sdf\n", encoding="utf-8")
        for output_format in ("selection", "alignment", "sdf"):
            with self.subTest(format=output_format):
                output = self.directory / f"bad-{output_format}"
                result = subprocess.run(
                    [sys.executable, str(ROOT / "selection" / "selector.py"), "--csv", str(self.csv),
                     "--format", output_format, "-o", str(output)],
                    cwd=ROOT, capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("CSV line 3", result.stderr)
                self.assertFalse(output.exists())

    def test_combined_cli_complete_reference_batch(self):
        output = self.directory / "reference-combined.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "selector.py"), "--csv",
             str(ROOT / "selection" / "reference_molecules.csv"),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        molecules = list(MoleculeSelector._read_sdf(output))
        self.assertEqual([mol.GetProp("_Name") for mol in molecules], ["241", "175", "878", "1068", "6582"])
        axes = {mol.GetProp("_Name"): json.loads(mol.GetProp(SDF_AXES_PROPERTY)) for mol in molecules}
        self.assertEqual(len(axes["241"]), 2)
        self.assertEqual(len(axes["175"]), 6)
        self.assertEqual(len(axes["878"]), 6)
        self.assertEqual(len(axes["1068"]), 5)
        thiol = next(mol for mol in molecules if mol.GetProp("_Name") == "878")
        original = next(MoleculeSelector._read_sdf(ROOT / "selection" / "reference_molecules" / "878.sdf"))
        sulfur = original.GetConformer().GetAtomPosition(0)
        hydrogen = original.GetConformer().GetAtomPosition(5)
        sh = (sulfur.x-hydrogen.x, sulfur.y-hydrogen.y, sulfur.z-hydrogen.z)
        length = math.sqrt(dot(sh, sh))
        target = tuple(x / length for x in sh)
        entry = next(item for item in axes["878"]
                     if item["sources"][0]["axis_type"] == "chalcogen_sigma"
                     and dot(item["axis"], target) > 0.999999)
        self.assertEqual(entry["sources"][0]["center_atom_index"], 0)
        marker = thiol.GetConformer().GetAtomPosition(entry["pseudo_atom_index"])
        for value, center, axis in zip(marker, sulfur, target):
            self.assertAlmostEqual(value, center + axis, delta=0.00006)

    def test_halogen_and_hydrogen_sigma_extend_bonds(self):
        mol = self.molecule("CCl")
        for smarts, axis_type, center in (("[Cl]", "halogen_sigma", 1),
                                          ("[H]", "hydrogen_sigma", 2)):
            with self.subTest(smarts=smarts):
                rows = self.prepare(mol, smarts=smarts)
                row = next(r for r in rows if r["center_atom_index"] == center)
                self.assertEqual(row["axis_type"], axis_type)
                axis = json.loads(row["axis"])
                self.assert_unit(axis)
                conf = mol.GetConformer()
                a, b = conf.GetAtomPosition(center), conf.GetAtomPosition(0)
                displacement = [a.x - b.x, a.y - b.y, a.z - b.z]
                self.assertAlmostEqual(dot(axis, displacement) / math.sqrt(dot(displacement, displacement)),
                                       1.0, places=6)

    def test_heteroatom_donor_uses_explicit_hydrogen_center(self):
        mol = self.molecule("CO")
        rows = self.prepare(mol, name="hydrogen bond donor")
        self.assertEqual([row["axis_type"] for row in rows], ["chalcogen_sigma", "hydrogen_sigma"])
        hydrogen = rows[1]["center_atom_index"]
        self.assertEqual(mol.GetAtomWithIdx(hydrogen).GetSymbol(), "H")

    def test_nitrile_ring_and_imine_lone_pairs(self):
        for smiles, smarts, kind in (
            ("CC#N", "[N]", "nitrogen_triple_sigma"),
            ("c1ccncc1", "[n]", "nitrogen_ring_lone_pair"),
            ("CC=N", "[N]", "nitrogen_imine_lone_pair"),
        ):
            with self.subTest(smiles=smiles):
                mol = self.molecule(smiles)
                rows = self.prepare(mol, smarts=smarts)
                mol = next(MoleculeSelector._read_sdf(self.sdf))
                row = next(r for r in rows if r["axis_type"] == kind)
                axis = json.loads(row["axis"])
                self.assert_unit(axis)
                n = mol.GetAtomWithIdx(row["center_atom_index"])
                conf = mol.GetConformer()
                origin = conf.GetAtomPosition(n.GetIdx())
                neighbors = list(n.GetNeighbors())
                vectors = []
                for neighbor in neighbors:
                    p = conf.GetAtomPosition(neighbor.GetIdx())
                    d = (p.x-origin.x, p.y-origin.y, p.z-origin.z)
                    norm = math.sqrt(dot(d, d))
                    vectors.append(tuple(x / norm for x in d))
                if len(vectors) == 1:
                    self.assertAlmostEqual(dot(axis, vectors[0]), -1, places=6)
                else:
                    self.assertAlmostEqual(dot(axis, vectors[0]), dot(axis, vectors[1]), places=6)
                    self.assertLess(dot(axis, vectors[0]), 0)

    def test_aromatic_ring_has_opposite_normals_at_centroid(self):
        mol = self.molecule("c1ccccc1")
        rows = self.prepare(mol, name="aromatic ring")
        self.assertEqual(len(rows), 2)
        axes = [json.loads(row["axis"]) for row in rows]
        self.assertEqual([row["axis_type"] for row in rows], ["aromatic_normal"] * 2)
        self.assertIsNone(rows[0]["center_atom_index"])
        self.assertEqual(rows[0]["center"], rows[1]["center"])
        self.assertAlmostEqual(dot(axes[0], axes[1]), -1.0, places=7)
        self.assert_unit(axes[0])
        center = json.loads(rows[0]["center"])
        conf = mol.GetConformer()
        for index in range(6):
            point = conf.GetAtomPosition(index)
            distance = [point.x-center[0], point.y-center[1], point.z-center[2]]
            self.assertAlmostEqual(dot(axes[0], distance), 0, places=2)

    def test_fused_rings_each_have_two_normals(self):
        rows = self.prepare(self.molecule("c1ccc2ccccc2c1"), name="aromatic ring")
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["match_index"] for row in rows}, {"0", "1"})
        for axis in (json.loads(row["axis"]) for row in rows):
            self.assert_unit(axis)

    def test_unsupported_and_no_matches(self):
        with self.assertLogs("molecule_aligner.selection.alignment", level="WARNING") as logs:
            rows = self.prepare(self.molecule("C[N+](C)(C)C"), name="positive charge")
        self.assertEqual(rows, [])
        self.assertIn("not yet supported", "\n".join(logs.output))
        self.assertEqual(self.prepare(self.molecule("CC"), name="carbonyl"), [])

    def test_hydrogen_without_coordinates_is_not_invented(self):
        mol = self.molecule("CO", explicit_hydrogens=False)
        with self.assertLogs("molecule_aligner.selection.alignment", level="WARNING") as logs:
            rows = self.prepare(mol, name="hydrogen bond donor")
        self.assertEqual([row["axis_type"] for row in rows], ["chalcogen_sigma"])
        self.assertIn("explicit hydrogen coordinates", "\n".join(logs.output))

    def test_singly_bonded_chalcogen_axes_extend_each_heavy_bond(self):
        examples = (
            ("CSC", "thioether", 2),
            ("COC", "ether", 2),
            ("C[Se]C", "chalcogen bond donor", 2),
            ("C[Te]C", "chalcogen bond donor", 2),
            ("CO", "alcohol", 1),
            ("COSC", None, 2),
        )
        for smiles, name, count in examples:
            with self.subTest(smiles=smiles):
                molecule = self.molecule(smiles)
                rows = self.prepare(molecule, name=name, smarts="[S]" if name is None else None)
                molecule = next(MoleculeSelector._read_sdf(self.sdf))
                rows = [row for row in rows if row["axis_type"] == "chalcogen_sigma"]
                self.assertEqual(len(rows), count)
                index = rows[0]["center_atom_index"]
                atom = molecule.GetAtomWithIdx(index)
                conformer = molecule.GetConformer()
                center = conformer.GetAtomPosition(index)
                expected = []
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if bond.GetBondType() != Chem.BondType.SINGLE or neighbor.GetAtomicNum() <= 1:
                        continue
                    p = conformer.GetAtomPosition(neighbor.GetIdx())
                    vector = (center.x - p.x, center.y - p.y, center.z - p.z)
                    length = math.sqrt(dot(vector, vector))
                    expected.append(tuple(component / length for component in vector))
                self.assertEqual(len(expected), count)
                for row, vector in zip(rows, expected):
                    axis = json.loads(row["axis"])
                    self.assert_unit(axis)
                    for actual, target in zip(axis, vector):
                        self.assertAlmostEqual(actual, target, places=7)
                    self.assertEqual(json.loads(row["center"]), list(center))

    def test_thioether_lone_pair_directions(self):
        rows = self.prepare(self.molecule("CSC"), name="thioether")
        self.assertEqual([row["axis_type"] for row in rows], [
            "chalcogen_sigma", "chalcogen_sigma",
            "thioether_lone_pair_normal", "thioether_lone_pair_normal",
            "thioether_outward_bisector",
        ])
        mol = next(MoleculeSelector._read_sdf(self.sdf))
        conformer = mol.GetConformer()
        center = conformer.GetAtomPosition(1)
        neighbors = []
        for index in (0, 2):
            position = conformer.GetAtomPosition(index)
            vector = (position.x - center.x, position.y - center.y, position.z - center.z)
            length = math.sqrt(dot(vector, vector))
            neighbors.append(tuple(value / length for value in vector))
        axes = [json.loads(row["axis"]) for row in rows]
        for row, axis in zip(rows, axes):
            self.assertEqual(row["center_atom_index"], 1)
            self.assertEqual(json.loads(row["center"]), list(center))
            self.assert_unit(axis)
        self.assertAlmostEqual(dot(axes[2], axes[3]), -1, places=7)
        for normal in axes[2:4]:
            for vector in neighbors:
                self.assertAlmostEqual(dot(normal, vector), 0, places=7)
        self.assertAlmostEqual(dot(axes[4], axes[2]), 0, places=7)
        self.assertAlmostEqual(dot(axes[4], neighbors[0]), dot(axes[4], neighbors[1]), places=7)
        self.assertLess(dot(axes[4], neighbors[0]), 0)

    def test_thiol_lone_pair_and_hydrogen_directions(self):
        rows = self.prepare(self.molecule("CS"), name="thiol")
        self.assertEqual([row["axis_type"] for row in rows], [
            "chalcogen_sigma", "chalcogen_sigma", "thiol_lone_pair_normal", "thiol_lone_pair_normal",
            "thiol_outward_bisector", "hydrogen_sigma",
        ])
        mol = next(MoleculeSelector._read_sdf(self.sdf))
        sulfur = mol.GetAtomWithIdx(1)
        neighbors = {atom.GetSymbol(): atom.GetIdx() for atom in sulfur.GetNeighbors()}
        conformer = mol.GetConformer()
        center = conformer.GetAtomPosition(1)
        vectors = []
        for symbol in ("C", "H"):
            p = conformer.GetAtomPosition(neighbors[symbol])
            vector = (p.x-center.x, p.y-center.y, p.z-center.z)
            magnitude = math.sqrt(dot(vector, vector))
            vectors.append(tuple(value / magnitude for value in vector))
        axes = [json.loads(row["axis"]) for row in rows]
        for row, axis in zip(rows[:5], axes[:5]):
            self.assertEqual(row["center_atom_index"], 1)
            self.assertEqual(json.loads(row["center"]), list(center))
            self.assert_unit(axis)
        self.assertAlmostEqual(dot(axes[0], vectors[0]), -1, places=7)
        self.assertAlmostEqual(dot(axes[1], vectors[1]), -1, places=7)
        self.assertAlmostEqual(dot(axes[2], axes[3]), -1, places=7)
        for normal in axes[2:4]:
            for vector in vectors:
                self.assertAlmostEqual(dot(normal, vector), 0, places=7)
        self.assertAlmostEqual(dot(axes[4], axes[2]), 0, places=7)
        self.assertAlmostEqual(dot(axes[4], vectors[0]), dot(axes[4], vectors[1]), places=7)
        self.assertLess(dot(axes[4], vectors[0]), 0)
        self.assertEqual(rows[5]["center_atom_index"], neighbors["H"])
        self.assertAlmostEqual(dot(axes[5], vectors[1]), 1, places=7)

    def test_thiol_implicit_hydrogen_warns_and_preserves_known_axis(self):
        mol = self.molecule("CS", explicit_hydrogens=False)
        with self.assertLogs("molecule_aligner.selection.alignment", level="WARNING") as logs:
            rows = self.prepare(mol, name="thiol")
        self.assertEqual([row["axis_type"] for row in rows], ["chalcogen_sigma"])
        self.assertIn("thiol plane", "\n".join(logs.output))
        self.assertIn("S-H sigma extension", "\n".join(logs.output))
        self.assertIn("explicit hydrogen coordinates", "\n".join(logs.output))

    def test_thiol_degenerate_plane_is_an_error(self):
        mol = self.molecule("CS")
        hydrogen = next(n.GetIdx() for n in mol.GetAtomWithIdx(1).GetNeighbors() if n.GetSymbol() == "H")
        conformer = mol.GetConformer()
        conformer.SetAtomPosition(0, (-1.0, 0.0, 0.0))
        conformer.SetAtomPosition(1, (0.0, 0.0, 0.0))
        conformer.SetAtomPosition(hydrogen, (1.0, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "degenerate"):
            self.prepare(mol, name="thiol")

    def test_thiol_sdf_markers_deduplicate_across_selectors(self):
        self.prepare(self.molecule("CS"), name="thiol")
        with self.csv.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=INPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(rows + [dict(rows[0], selector="hydrogen bond donor")])
        output = self.directory / "thiol.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        mol = next(MoleculeSelector._read_sdf(output))
        axes = json.loads(mol.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 6)
        for entry in axes:
            self.assertEqual([s["selector"] for s in entry["sources"]], ["thiol", "hydrogen bond donor"])
            position = mol.GetConformer().GetAtomPosition(entry["pseudo_atom_index"])
            self.assertEqual(mol.GetAtomWithIdx(entry["pseudo_atom_index"]).GetAtomicNum(), 0)
            for value, center, axis in zip(position, entry["center"], entry["axis"]):
                self.assertAlmostEqual(value, center + axis, delta=0.00006)
            distance = math.sqrt(sum((value - center) ** 2 for value, center in zip(position, entry["center"])))
            self.assertAlmostEqual(distance, 1.0, delta=0.0001)

    def test_thioether_sdf_has_five_markers_and_deduplicates_selectors(self):
        self.prepare(self.molecule("CSC"), name="thioether")
        with self.csv.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=INPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(rows + [dict(rows[0], selector="chalcogen bond donor")])
        output = self.directory / "thioether.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        mol = next(MoleculeSelector._read_sdf(output))
        axes = json.loads(mol.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 5)
        for entry in axes:
            self.assertEqual([s["selector"] for s in entry["sources"]],
                             ["thioether", "chalcogen bond donor"])
            position = mol.GetConformer().GetAtomPosition(entry["pseudo_atom_index"])
            self.assertEqual(mol.GetAtomWithIdx(entry["pseudo_atom_index"]).GetAtomicNum(), 0)
            for value, center, axis in zip(position, entry["center"], entry["axis"]):
                self.assertAlmostEqual(value, center + axis, delta=0.00006)
            self.assertAlmostEqual(
                math.sqrt(sum((value - center) ** 2 for value, center in zip(position, entry["center"]))),
                1.0, delta=0.0001,
            )

    def test_thioether_extra_axes_only_apply_to_carbon_sulfur_carbon(self):
        for smiles in ("CS", "COSC", "CS(=O)C", "COC"):
            with self.subTest(smiles=smiles):
                rows = self.prepare(self.molecule(smiles), smarts="[S,O]")
                self.assertFalse(any(row["axis_type"].startswith("thioether_") for row in rows))

    def test_degenerate_thioether_plane_is_an_error(self):
        mol = self.molecule("CSC")
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, (-1.0, 0.0, 0.0))
        conf.SetAtomPosition(1, (0.0, 0.0, 0.0))
        conf.SetAtomPosition(2, (1.0, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "degenerate"):
            self.prepare(mol, name="thioether")

    def test_invalid_csv_indices_and_2d_geometry_fail(self):
        mol = Chem.MolFromSmiles("CCl")
        AllChem.Compute2DCoords(mol)
        with self.assertRaisesRegex(ValueError, "3D conformer"):
            self.prepare(mol, smarts="[Cl]")
        mol = self.molecule("CCl")
        self.prepare(mol, smarts="[Cl]")
        with self.csv.open(newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source)
            rows = list(reader)
        rows[0]["site_atom_indices"] = "[99]"
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=INPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(ValueError, "CSV line 2: site_atom_indices"):
            AlignmentCsv().read(self.csv)

    def test_degenerate_geometry_is_an_error(self):
        mol = self.molecule("CCl")
        conf = mol.GetConformer()
        conf.SetAtomPosition(1, conf.GetAtomPosition(0))
        with self.assertRaisesRegex(ValueError, "degenerate"):
            self.prepare(mol, smarts="[Cl]")

    def test_multiple_sdf_records_select_correct_conformer(self):
        first = self.molecule("CC")
        second = self.molecule("CCl")
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(first)
        writer.write(second)
        writer.close()
        selected = MoleculeSelector().select(self.sdf, smarts="[Cl]")
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            CsvResultWriter().write(selected, destination, sdf_file=self.sdf)
        rows = AlignmentCsv().read(self.csv)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["record_index"], "1")
        self.assertEqual(rows[0]["axis_type"], "halogen_sigma")

    def test_cli_fails_without_partial_output(self):
        mol = Chem.MolFromSmiles("CCl")
        AllChem.Compute2DCoords(mol)
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(mol)
        writer.close()
        selected = MoleculeSelector().select(self.sdf, smarts="[Cl]")
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            CsvResultWriter().write(selected, destination, sdf_file=self.sdf)
        output = self.directory / "axes.csv"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("3D conformer", result.stderr)
        self.assertFalse(output.exists())

    def test_cli_writes_relative_paths_and_refuses_overwrite(self):
        self.prepare(self.molecule("CCl"), smarts="[Cl]")
        output = self.directory / "axes.csv"
        command = [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
                   "--output", str(output)]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open(newline="", encoding="utf-8") as source:
            row = next(csv.DictReader(source))
        self.assertFalse(os.path.isabs(row["sdf_file"]))
        self.assertEqual((ROOT / row["sdf_file"]).resolve(), self.sdf.resolve())
        before = output.read_text()
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output.read_text(), before)

    def test_sdf_output_deduplicates_axes_and_preserves_provenance(self):
        mol = self.molecule("CC(=O)C")
        mol.SetProp("_Name", "acetone")
        mol.SetProp("custom", "preserved")
        self.prepare(mol, name="carbonyl")
        with self.csv.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        duplicate = dict(rows[0], selector="hydrogen bond acceptor")
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=INPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(rows + [duplicate])

        output = self.directory / "axes.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        records = list(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(records), 1)
        annotated = records[0]
        self.assertEqual(annotated.GetProp("_Name"), "acetone")
        self.assertEqual(annotated.GetProp("custom"), "preserved")
        self.assertEqual(annotated.GetNumAtoms(), mol.GetNumAtoms() + 3)
        self.assertEqual(annotated.GetNumBonds(), mol.GetNumBonds())
        provenance = json.loads(annotated.GetProp(SDF_SOURCE_PROPERTY))
        self.assertEqual((ROOT / provenance["sdf_file"]).resolve(), self.sdf.resolve())
        self.assertEqual(provenance["record_index"], 0)
        axes = json.loads(annotated.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 3)
        conf = annotated.GetConformer()
        for entry in axes:
            idx = entry["pseudo_atom_index"]
            self.assertEqual(annotated.GetAtomWithIdx(idx).GetAtomicNum(), 0)
            self.assertEqual(annotated.GetAtomWithIdx(idx).GetDegree(), 0)
            self.assertEqual([s["selector"] for s in entry["sources"]],
                             ["carbonyl", "hydrogen bond acceptor"])
            position = conf.GetAtomPosition(idx)
            for value, center, axis in zip(position, entry["center"], entry["axis"]):
                self.assertAlmostEqual(value, center + axis, delta=0.00006)
            self.assertAlmostEqual(
                math.sqrt(sum((value - center) ** 2 for value, center in zip(position, entry["center"]))),
                1.0, delta=0.0001,
            )
        original = next(MoleculeSelector._read_sdf(self.sdf))
        for index in range(original.GetNumAtoms()):
            self.assertEqual(tuple(original.GetConformer().GetAtomPosition(index)),
                             tuple(conf.GetAtomPosition(index)))

        again = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(again.returncode, 0)
        self.assertEqual(len(list(MoleculeSelector._read_sdf(output))), 1)

    def test_sdf_output_ring_centroid_and_unmatched_multirecord(self):
        ring, unmatched = self.molecule("c1ccccc1"), self.molecule("CC")
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(ring)
        writer.write(unmatched)
        writer.close()
        selected = MoleculeSelector().select(self.sdf, name="aromatic ring")
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            CsvResultWriter().write(selected, destination, sdf_file=self.sdf)
        output = self.directory / "rings.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("no supported alignment axes", result.stderr)
        records = list(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].GetNumAtoms(), ring.GetNumAtoms() + 2)
        self.assertEqual(records[1].GetNumAtoms(), unmatched.GetNumAtoms())
        self.assertEqual(json.loads(records[1].GetProp(SDF_AXES_PROPERTY)), [])
        axes = json.loads(records[0].GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 2)
        for axis in axes:
            self.assertIsNone(axis["sources"][0]["center_atom_index"])
            position = records[0].GetConformer().GetAtomPosition(axis["pseudo_atom_index"])
            for value, center, direction in zip(position, axis["center"], axis["axis"]):
                self.assertAlmostEqual(value, center + direction, delta=0.00006)

    def test_sdf_error_does_not_create_output(self):
        mol = Chem.MolFromSmiles("CCl")
        AllChem.Compute2DCoords(mol)
        writer = Chem.SDWriter(str(self.sdf))
        writer.write(mol)
        writer.close()
        selected = MoleculeSelector().select(self.sdf, smarts="[Cl]")
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            CsvResultWriter().write(selected, destination, sdf_file=self.sdf)
        output = self.directory / "axes.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("3D conformer", result.stderr)
        self.assertFalse(output.exists())

    def test_sdf_output_includes_each_input_file(self):
        first = self.molecule("CCl")
        second = self.molecule("CBr")
        second_path = self.directory / "second.sdf"
        for path, mol in ((self.sdf, first), (second_path, second)):
            writer = Chem.SDWriter(str(path))
            writer.write(mol)
            writer.close()
        all_rows = []
        for path, smarts in ((self.sdf, "[Cl]"), (second_path, "[Br]")):
            selected = MoleculeSelector().select(path, smarts=smarts)
            stream = io.StringIO()
            CsvResultWriter().write(selected, stream, sdf_file=path)
            all_rows.extend(csv.DictReader(io.StringIO(stream.getvalue())))
        with self.csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=INPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(all_rows)
        output = self.directory / "multiple.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        records = list(MoleculeSelector._read_sdf(output))
        self.assertEqual(len(records), 2)
        for record, path in zip(records, (self.sdf, second_path)):
            source = json.loads(record.GetProp(SDF_SOURCE_PROPERTY))
            self.assertEqual((ROOT / source["sdf_file"]).resolve(), path.resolve())
            self.assertEqual(len(json.loads(record.GetProp(SDF_AXES_PROPERTY))), 1)

    def test_sdf_rejects_reserved_property_without_overwriting(self):
        mol = self.molecule("CCl")
        mol.SetProp(SDF_AXES_PROPERTY, "existing data")
        self.prepare(mol, smarts="[Cl]")
        output = self.directory / "collision.sdf"
        result = subprocess.run(
            [sys.executable, str(ROOT / "selection" / "alignment.py"), str(self.csv),
             "--format", "sdf", "-o", str(output)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("already has property", result.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
