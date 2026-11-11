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

from selection.alignment import AlignmentCalculator, AlignmentCsv, INPUT_COLUMNS
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
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["axis_type"], "hydrogen_sigma")
        hydrogen = rows[0]["center_atom_index"]
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
        with self.assertLogs("quick_qm.selection.alignment", level="WARNING") as logs:
            rows = self.prepare(self.molecule("CSC"), name="thioether")
        self.assertEqual(rows, [])
        self.assertIn("not yet supported", "\n".join(logs.output))
        self.assertEqual(self.prepare(self.molecule("CC"), name="carbonyl"), [])

    def test_hydrogen_without_coordinates_is_not_invented(self):
        mol = self.molecule("CO", explicit_hydrogens=False)
        with self.assertLogs("quick_qm.selection.alignment", level="WARNING") as logs:
            rows = self.prepare(mol, name="hydrogen bond donor")
        self.assertEqual(rows, [])
        self.assertIn("explicit hydrogen coordinates", "\n".join(logs.output))

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


if __name__ == "__main__":
    unittest.main()
