import csv
import json
import math
import unittest
from unittest.mock import patch

import numpy as np
import yaml
from rdkit import Chem

from aligner import test_aligner
from aligner.torsion_optimization import TorsionOptimizer
from molecular_selector.matching import MolecularSelector


class TorsionTests(unittest.TestCase):
    setUp = test_aligner.AlignerTests.setUp
    write_rows = test_aligner.AlignerTests.write_rows
    run_cli = test_aligner.AlignerTests.run_cli
    filters = test_aligner.AlignerTests.filters
    def run_optimization(self, name, config, *extra):
        path = self.directory / f"{name}.yaml"
        path.write_text(yaml.safe_dump({"torsion_optimization": config}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--distance", 3, "--x-distance", 1,
                              "--config", path, "--output", name, *extra)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / f"{name}.csv").open(newline="") as source:
            row = next(csv.DictReader(source))
        pair = next(MolecularSelector._read_sdf(self.directory / row["sdf_file"]))
        return row, pair

    def test_distance_scores_and_transform(self):
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        flexible = next(MolecularSelector._read_sdf(self.paths[1]))
        for score_type in ("average_distance", "minimum_distance"):
            for objective in ("minimize", "maximize"):
                with self.subTest(score=score_type, objective=objective):
                    name = f"{score_type}-{objective}"
                    row, pair = self.run_optimization(name, {
                        "enabled": True, "score": score_type, "objective": objective,
                        "steps_per_iteration": 12, "iterations": 1,
                        "exclude_rigid_site": True, "include_vdw": True,
                    }, "--flexible-orientation", "toward-center")
                    matrix = np.array(json.loads(row["rotation"]))
                    translation = np.array(json.loads(row["translation"]))
                    rcoords = np.array(pair.GetConformer().GetPositions()[:rigid.GetNumAtoms()])
                    fcoords = np.array(pair.GetConformer().GetPositions()[rigid.GetNumAtoms():])
                    for index, original in enumerate(flexible.GetConformer().GetPositions()):
                        np.testing.assert_allclose(fcoords[index], matrix @ original + translation, atol=0.00006)
                    center = np.array(json.loads(row["flexible_aligned_center"]))
                    np.testing.assert_allclose(fcoords[1], center, atol=0.00006)
                    direction = np.array(json.loads(row["rigid_center"])) - center
                    np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]),
                                               direction / np.linalg.norm(direction), atol=1e-12)
                    ri = [i for i in range(rigid.GetNumAtoms()) if i != 1]
                    fi = [i for i in range(flexible.GetNumAtoms()) if i != 1]
                    distances = np.linalg.norm(rcoords[ri, None, :] - fcoords[None, fi, :], axis=2)
                    pt = Chem.GetPeriodicTable()
                    distances -= (np.array([pt.GetRvdw(rigid.GetAtomWithIdx(i).GetAtomicNum()) for i in ri])[:, None]
                                  + np.array([pt.GetRvdw(flexible.GetAtomWithIdx(i).GetAtomicNum()) for i in fi]))
                    expected = np.mean(distances) if score_type == "average_distance" else np.min(distances)
                    self.assertAlmostEqual(float(row["torsion_score"]), expected, delta=0.0002)
                    self.assertEqual(pair.GetProp("quick_qm_torsion_score_type"), score_type)
                    self.assertEqual(row["torsion_objective"], objective)

    def test_mmff_interaction_matches_independent_rdkit_calculation(self):
        row, pair = self.run_optimization("energy", {
            "enabled": True, "score": "mmff94_interaction", "objective": "minimize",
            "steps_per_iteration": 6, "iterations": 1,
        })
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        flexible = next(MolecularSelector._read_sdf(self.paths[1]))
        from rdkit.Chem import AllChem

        def energy(mol, inter):
            props = AllChem.MMFFGetMoleculeProperties(mol)
            self.assertIsNotNone(props)
            return AllChem.MMFFGetMoleculeForceField(
                mol, props, ignoreInterfragInteractions=not inter, nonBondedThresh=1e6).CalcEnergy()

        actual = energy(pair, True) - energy(rigid, False) - energy(flexible, False)
        self.assertAlmostEqual(float(row["torsion_score"]), actual, delta=0.03)
        self.assertEqual(row["torsion_score_type"], "mmff94_interaction")

    def test_missing_mmff_parameters_are_an_error(self):
        molecule = next(MolecularSelector._read_sdf(self.paths[0]))
        with patch("aligner.torsion_optimization.AllChem.MMFFGetMoleculeProperties", return_value=None):
            with self.assertRaisesRegex(ValueError, "MMFF94 parameters unavailable"):
                TorsionOptimizer.mmff_energy(molecule, False)

    def test_score_optimization_improves_baseline_and_respects_overrides(self):
        for objective in ("minimize", "maximize"):
            row, _ = self.run_optimization(objective, {
                "enabled": True, "score": "average_distance", "objective": objective,
                "steps_per_iteration": 12, "iterations": 2,
            })
            baseline, baseline_pair = self.run_optimization(f"base-{objective}", {
                "enabled": True, "score": "average_distance", "objective": objective,
                "steps_per_iteration": 3, "iterations": 1,
            }, "--no-optimize-torsion")
            self.assertEqual(baseline["torsion_score"], "")
            coords = baseline_pair.GetConformer().GetPositions()
            rigid_count = next(MolecularSelector._read_sdf(self.paths[0])).GetNumAtoms()
            flexible_indices = [i for i in range(len(coords) - rigid_count) if i != 1]
            baseline_score = np.mean(np.linalg.norm(coords[:rigid_count, None, :] -
                                                   coords[None, rigid_count + np.array(flexible_indices), :], axis=2))
            score = float(row["torsion_score"])
            self.assertTrue(math.isfinite(score))
            self.assertTrue(score < baseline_score - 1e-4 if objective == "minimize"
                            else score > baseline_score + 1e-4)
        config = self.directory / "invalid.yaml"
        config.write_text(yaml.safe_dump({"torsion_optimization": {"enabled": True, "score": "invalid"}}),
                          encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("torsion_optimization.score", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs.csv").exists())

    def test_site_exclusions_and_missing_site_errors(self):
        row = type("Row", (), {"line": 2, "center_index": 1,
                               "fields": {"site_atom_indices": "[1,2]", "atom_index_base": "0"}})()
        mol = next(MolecularSelector._read_sdf(self.paths[0]))
        indices = TorsionOptimizer.excluded(row, mol, True, [0])
        self.assertEqual(indices, list(range(3, mol.GetNumAtoms())))
        row.fields["site_atom_indices"] = "[999]"
        with self.assertRaisesRegex(ValueError, "outside"):
            TorsionOptimizer.excluded(row, mol, True, [])
        row.fields["site_atom_indices"] = ""
        self.assertNotIn(1, TorsionOptimizer.excluded(row, mol, True, []))
        row.center_index = None
        with self.assertRaisesRegex(ValueError, "required"):
            TorsionOptimizer.excluded(row, mol, True, [])
        with self.assertRaisesRegex(ValueError, "no atoms"):
            TorsionOptimizer.excluded(row, mol, False, list(range(mol.GetNumAtoms())))

    def test_selected_score_is_best_of_first_pass(self):
        baseline, pair = self.run_optimization("unoptimized", {"enabled": False})
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        coordinates = pair.GetConformer().GetPositions()
        rigid_count = rigid.GetNumAtoms()
        fixed = coordinates[:rigid_count]
        moving = coordinates[rigid_count:]
        center = np.array(json.loads(baseline["flexible_aligned_center"]))
        axis = np.array(json.loads(baseline["flexible_aligned_axis"]))
        mask = [i for i in range(len(moving)) if i != 1]
        scores = []
        for i in range(12):
            matrix = TorsionOptimizer.rotation(axis, 2 * math.pi * i / 12)
            rotated = (moving - center) @ matrix.T + center
            scores.append(np.min(np.linalg.norm(fixed[:, None, :] - rotated[None, mask, :], axis=2)))
        for objective, best in (("minimize", min), ("maximize", max)):
            row, _ = self.run_optimization(f"opt-{objective}", {
                "enabled": True, "objective": objective, "iterations": 1, "steps_per_iteration": 12,
            })
            self.assertAlmostEqual(float(row["torsion_score"]), best(scores), delta=0.0002)

    def test_polar_poses_and_failure_do_not_replace_existing_output(self):
        path = self.directory / "polar.yaml"
        path.write_text(yaml.safe_dump({
            "polar_coordinate_based_scan": {
                "radius": {"distance": 3},
                "phi": {"scan_start": 0, "scan_end": 180, "scan_step": 90},
                "theta": {"distance": 90},
            },
            "torsion_optimization": {"enabled": True, "iterations": 1, "steps_per_iteration": 3},
        }), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", path, "--output", "polar")
        self.assertEqual(result.returncode, 0, result.stderr)
        metadata_path = self.directory / "polar.csv"
        with metadata_path.open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["coordinate_system"] == "polar" and row["torsion_score"] for row in rows))
        before = metadata_path.read_bytes()
        output = self.directory / rows[0]["sdf_file"]
        sdf_before = output.read_bytes()
        result = self.run_cli(*self.filters(), "--config", path, "--output", "polar",
                              "--overwrite-output", "--torsion-exclude-flexible-atoms", "999")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("excluded atom index", result.stderr)
        self.assertEqual(metadata_path.read_bytes(), before)
        self.assertEqual(output.read_bytes(), sdf_before)


if __name__ == "__main__":
    unittest.main()
