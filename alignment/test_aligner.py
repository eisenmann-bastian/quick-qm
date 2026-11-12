import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from alignment import aligner
from alignment.aligner import PairAligner
from selection.matching import MoleculeSelector


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "alignment" / "aligner.py"


class AlignerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.input = self.directory / "axes.csv"
        self.paths = []
        for name, smiles in (("rigid", "CO"), ("flexible", "CCl")):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=19), 0)
            path = self.directory / f"{name}.sdf"
            writer = Chem.SDWriter(str(path))
            mol.SetProp("_Name", name)
            writer.write(mol)
            writer.close()
            self.paths.append(path)
        self.write_rows()

    def write_rows(self):
        fields = ["sdf_file", "record_index", "molecule_name", "selector", "axis_type",
                  "center_atom_index", "center", "axis"]
        rows = []
        for path, name, kind, axis in ((self.paths[0], "rigid", "oxygen", [1, 0, 0]),
                                        (self.paths[1], "flexible", "halogen", [0, 1, 0])):
            mol = next(MoleculeSelector._read_sdf(path))
            p = mol.GetConformer().GetAtomPosition(1)
            rows.append(dict(sdf_file=str(path), record_index=0, molecule_name=name, selector=kind,
                             axis_type="sigma", center_atom_index=1,
                             center=json.dumps([p.x, p.y, p.z]), axis=json.dumps(axis)))
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        return rows

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), str(self.input), *map(str, args)],
                              cwd=self.directory, capture_output=True, text=True)

    def filters(self):
        return ["--rigid", "^rigid$", "oxygen", "sigma",
                "--flexible", "^flexible$", "halogen", "sigma"]

    def test_pair_file_aligns_rotation_and_vdw_distance(self):
        result = self.run_cli(*self.filters())
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["sdf_record_index"], "0")
        expected = Chem.GetPeriodicTable().GetRvdw(8) + Chem.GetPeriodicTable().GetRvdw(17)
        self.assertAlmostEqual(float(row["distance_angstrom"]), expected)
        rigid_axis = np.array(json.loads(row["rigid_axis"]))
        flex_axis = np.array(json.loads(row["flexible_aligned_axis"]))
        np.testing.assert_allclose(flex_axis, -rigid_axis, atol=1e-12)
        rcenter = np.array(json.loads(row["rigid_center"]))
        fcenter = np.array(json.loads(row["flexible_aligned_center"]))
        np.testing.assert_allclose(fcenter, rcenter + expected * rigid_axis, atol=1e-12)
        pair = next(MoleculeSelector._read_sdf(self.directory / row["sdf_file"]))
        rigid = next(MoleculeSelector._read_sdf(self.paths[0]))
        flexible = next(MoleculeSelector._read_sdf(self.paths[1]))
        self.assertEqual(pair.GetNumAtoms(), rigid.GetNumAtoms() + flexible.GetNumAtoms())
        self.assertEqual(pair.GetNumBonds(), rigid.GetNumBonds() + flexible.GetNumBonds())
        self.assertEqual(len(Chem.GetMolFrags(pair)), 2)
        for index in range(rigid.GetNumAtoms()):
            np.testing.assert_allclose(pair.GetConformer().GetAtomPosition(index),
                                       rigid.GetConformer().GetAtomPosition(index), atol=1e-4)
        p = pair.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + 1)
        np.testing.assert_allclose(p, fcenter, atol=0.00006)
        self.assertEqual(pair.GetProp("quick_qm_pair_id"), row["pair_id"])

    def test_n_times_m_and_yaml_single_output(self):
        rows = self.write_rows()
        with self.input.open("a", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writerow(dict(rows[0], axis=json.dumps([-1, 0, 0])))
            writer.writerow(dict(rows[1], axis=json.dumps([0, -1, 0])))
        config = self.directory / "aligner.yaml"
        config.write_text(yaml.safe_dump({
            "rigid": [dict(molecule_name="rigid", selector="oxy", axis_type="sigma")],
            "flexible": [dict(molecule_name="flex", selector="hal", axis_type="sigma")],
            "output_format": "single", "output": "all.sdf", "distance": 4.5,
        }), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        molecules = list(MoleculeSelector._read_sdf(self.directory / "all.sdf"))
        with (self.directory / "all.csv").open(newline="", encoding="utf-8") as source:
            metadata = list(csv.DictReader(source))
        self.assertEqual(len(molecules), 4)
        self.assertEqual(len(metadata), 4)
        self.assertEqual([int(row["sdf_record_index"]) for row in metadata], [0, 1, 2, 3])
        self.assertTrue(all(float(row["distance_angstrom"]) == 4.5 for row in metadata))
        self.assertEqual({(row["rigid_csv_line"], row["flexible_csv_line"]) for row in metadata},
                         {("2", "3"), ("2", "5"), ("4", "3"), ("4", "5")})

    def test_rotation_parallel_and_antiparallel(self):
        for source, target in (([1, 0, 0], [1, 0, 0]), ([1, 0, 0], [-1, 0, 0]),
                               ([1, 0, 0], [0, 1, 0])):
            matrix = PairAligner.rotation(np.array(source, dtype=float), np.array(target, dtype=float))
            np.testing.assert_allclose(matrix @ source, target, atol=1e-12)
            np.testing.assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(matrix), 1)

    def test_aromatic_radius_and_fixed_override(self):
        rows = self.write_rows()
        for row in rows:
            row["center_atom_index"] = ""
            row["axis_type"] = "aromatic_normal"
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        args = ["--rigid", "rigid", ".*", "aromatic_normal",
                "--flexible", "flexible", ".*", "aromatic_normal"]
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 0)
        rows[1]["axis_type"] = "sigma"
        rows[1]["center_atom_index"] = "1"
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        result = self.run_cli("--rigid", "rigid", ".*", "aromatic_normal",
                              "--flexible", "flexible", ".*", "sigma", "-o", "mixed")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "mixed.csv").open(newline="") as source:
            self.assertAlmostEqual(float(next(csv.DictReader(source))["distance_angstrom"]),
                                   Chem.GetPeriodicTable().GetRvdw(17))
        result = self.run_cli("--rigid", "rigid", ".*", "aromatic_normal",
                              "--flexible", "flexible", ".*", "sigma", "-o", "fixed", "--distance", 7)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "fixed.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 7)

    def test_bad_input_and_existing_output_fail_without_overwrite(self):
        result = self.run_cli("--rigid", "[", ".*", ".*", "--flexible", ".*", ".*", ".*")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("regex", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs").exists())
        result = self.run_cli(*self.filters(), "--distance", "-1")
        self.assertNotEqual(result.returncode, 0)
        result = self.run_cli(*self.filters())
        self.assertEqual(result.returncode, 0, result.stderr)
        before = (self.directory / "aligned_pairs.csv").read_bytes()
        self.assertNotEqual(self.run_cli(*self.filters()).returncode, 0)
        self.assertEqual((self.directory / "aligned_pairs.csv").read_bytes(), before)

    def test_invalid_selected_axis_is_reported(self):
        rows = self.write_rows()
        rows[1]["axis"] = "[0,0,0]"
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        result = self.run_cli(*self.filters())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CSV line 3", result.stderr)
        self.assertIn("nonzero", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs").exists())

    def test_reference_csv_and_yaml_work_from_source_working_directory(self):
        source = ROOT / "alignment"
        output = self.directory / "reference"
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "test_out_aligned.csv", "--config",
             str(source / "aligner.config"), "--output", str(output)],
            cwd=source, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "reference.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertAlmostEqual(float(row["distance_angstrom"]), Chem.GetPeriodicTable().GetRvdw(1))
            np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]),
                                       -np.array(json.loads(row["rigid_axis"])), atol=1e-12)
        self.assertEqual(len(list(output.glob("*.sdf"))), 2)

    def test_config_list_filters_and_command_line_override(self):
        config = self.directory / "config.yaml"
        config.write_text(yaml.safe_dump({
            "rigid": {"molecule_name": "does-not-match", "selector": ".*", "axis_type": ".*"},
            "flexible": {"molecule_name": "flexible", "selector": ".*", "axis_type": ".*"},
            "distance": 1.0,
        }), encoding="utf-8")
        result = self.run_cli("--config", config, "--rigid", "rigid", ".*", ".*",
                              "--distance", "2.5", "-o", "override")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "override.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 2.5)

    def test_scan_step_endpoints_and_decimal_divisibility(self):
        for include_start in (True, False):
            for include_end in (True, False):
                name = f"step-{include_start}-{include_end}"
                result = self.run_cli(*self.filters(), "-o", name,
                                      "--alignment-mode", "distance-scan",
                                      "--scan-start", "0", "--scan-end", "0.3", "--scan-step", "0.1",
                                      "--include-start" if include_start else "--no-include-start",
                                      "--include-end" if include_end else "--no-include-end")
                self.assertEqual(result.returncode, 0, result.stderr)
                with (self.directory / f"{name}.csv").open(newline="") as source:
                    rows = list(csv.DictReader(source))
                expected = [0.0, 0.1, 0.2, 0.3][0 if include_start else 1:4 if include_end else 3]
                self.assertEqual([float(row["distance_angstrom"]) for row in rows], expected)
                self.assertEqual([int(row["pose_index"]) for row in rows], list(range(1, len(expected) + 1)))
                self.assertEqual({row["axis_pair_id"] for row in rows}, {"axis-pair-000001"})
                for row in rows:
                    pair = list(MoleculeSelector._read_sdf(self.directory / row["sdf_file"]))[int(row["sdf_record_index"])]
                    self.assertEqual(pair.GetProp("quick_qm_pose_index"), row["pose_index"])
                    self.assertEqual(pair.GetProp("quick_qm_distance_angstrom"), row["distance_angstrom"])
                    rigid = next(MoleculeSelector._read_sdf(self.paths[0]))
                    point = pair.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + 1)
                    np.testing.assert_allclose(point, json.loads(row["flexible_aligned_center"]), atol=0.00006)

    def test_scan_count_endpoint_grid_and_single_pose(self):
        for include_start, include_end, expected in (
            (True, True, [0, 0.6, 1.2, 1.8, 2.4, 3]),
            (True, False, [0, 0.6, 1.2, 1.8, 2.4]),
            (False, True, [0.6, 1.2, 1.8, 2.4, 3]),
            (False, False, [0.6, 1.2, 1.8, 2.4]),
        ):
            name = f"count-{include_start}-{include_end}"
            result = self.run_cli(*self.filters(), "-o", name, "--alignment-mode", "distance-scan",
                                  "--scan-start", 0, "--scan-end", 3, "--scan-count", 4,
                                  "--include-start" if include_start else "--no-include-start",
                                  "--include-end" if include_end else "--no-include-end")
            self.assertEqual(result.returncode, 0, result.stderr)
            with (self.directory / f"{name}.csv").open(newline="") as source:
                self.assertEqual([float(row["distance_angstrom"]) for row in csv.DictReader(source)], expected)
        result = self.run_cli(*self.filters(), "-o", "midpoint", "--alignment-mode", "distance-scan",
                              "--scan-start", 0, "--scan-end", 3, "--scan-count", 1,
                              "--no-include-start", "--no-include-end")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "midpoint.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 1.5)
        for start, end, expected in (
            (True, True, [0, 1.5, 3]),
            (True, False, [0, 1.5]),
            (False, True, [1.5, 3]),
        ):
            name = f"midpoint-{start}-{end}"
            result = self.run_cli(*self.filters(), "-o", name, "--alignment-mode", "distance-scan",
                                  "--scan-start", 0, "--scan-end", 3, "--scan-count", 1,
                                  "--include-start" if start else "--no-include-start",
                                  "--include-end" if end else "--no-include-end")
            self.assertEqual(result.returncode, 0, result.stderr)
            with (self.directory / f"{name}.csv").open(newline="") as source:
                self.assertEqual([float(row["distance_angstrom"]) for row in csv.DictReader(source)], expected)

    def test_invalid_scan_does_not_create_output(self):
        base = [*self.filters(), "--alignment-mode", "distance-scan"]
        cases = [
            (["--scan-start", "0", "--scan-end", "1", "--scan-step", "0.3"], "evenly divisible"),
            (["--scan-start", "0", "--scan-end", "1", "--scan-step", "0"], "positive"),
            (["--scan-start", "0", "--scan-end", "1", "--scan-step", "1",
              "--no-include-start", "--no-include-end"], "no poses"),
            (["--scan-start", "1", "--scan-end", "1", "--scan-count", "2"], "scan range"),
            (["--scan-start", "-1", "--scan-end", "1", "--scan-count", "2"], "scan range"),
            (["--scan-start", "nan", "--scan-end", "1", "--scan-count", "2"], "finite"),
            (["--scan-start", "0", "--scan-end", "1", "--scan-count", "0"], "positive integer"),
            (["--scan-start", "0", "--scan-end", "1"], "exactly one"),
            (["--scan-start", "0", "--scan-end", "1", "--scan-count", "2", "--distance", "2"], "cannot be used"),
        ]
        for options, message in cases:
            with self.subTest(options=options):
                result = self.run_cli(*base, *options)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse((self.directory / "aligned_pairs").exists())
                self.assertFalse((self.directory / "aligned_pairs.csv").exists())

    def test_scan_yaml_cli_override_and_cartesian_product(self):
        rows = self.write_rows()
        with self.input.open("a", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writerow(dict(rows[0], axis=json.dumps([-1, 0, 0])))
            writer.writerow(dict(rows[1], axis=json.dumps([0, -1, 0])))
        config = self.directory / "scan.yaml"
        config.write_text(yaml.safe_dump({
            "rigid": dict(molecule_name="rigid", selector=".*", axis_type=".*"),
            "flexible": dict(molecule_name="flexible", selector=".*", axis_type=".*"),
            "output_format": "single", "output": "scan.sdf", "alignment_mode": "distance-scan",
            "scan_start": 0, "scan_end": 3, "scan_step": 1,
            "include_start": False, "include_end": False,
        }), encoding="utf-8")
        result = self.run_cli("--config", config, "--scan-count", 2, "--include-start")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "scan.csv").open(newline="") as source:
            metadata = list(csv.DictReader(source))
        molecules = list(MoleculeSelector._read_sdf(self.directory / "scan.sdf"))
        self.assertEqual(len(metadata), 12)
        self.assertEqual(len(molecules), 12)
        self.assertEqual([int(row["sdf_record_index"]) for row in metadata], list(range(12)))
        for i in range(4):
            group = metadata[i * 3:(i + 1) * 3]
            self.assertEqual([float(row["distance_angstrom"]) for row in group], [0, 1, 2])
            self.assertEqual([int(row["pose_index"]) for row in group], [1, 2, 3])
            self.assertEqual({row["axis_pair_id"] for row in group}, {f"axis-pair-{i + 1:06d}"})
            self.assertEqual(molecules[i * 3].GetProp("quick_qm_axis_pair_id"), group[0]["axis_pair_id"])

    def test_scan_yaml_validation_and_step_override(self):
        config = self.directory / "scan.yaml"
        settings = {
            "rigid": dict(molecule_name="rigid", selector=".*", axis_type=".*"),
            "flexible": dict(molecule_name="flexible", selector=".*", axis_type=".*"),
            "alignment_mode": "distance-scan", "scan_start": 0, "scan_end": 1,
            "scan_count": 4, "include_start": "false",
        }
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("include_start must be a boolean", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs.csv").exists())

        settings["include_start"] = False
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config, "--scan-step", "0.5")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            self.assertEqual([float(row["distance_angstrom"]) for row in csv.DictReader(source)], [0.5, 1])

    def test_pair_files_group_scan_poses_per_axis_pair(self):
        rows = self.write_rows()
        with self.input.open("a", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writerow(dict(rows[0], axis=json.dumps([-1, 0, 0])))
            writer.writerow(dict(rows[1], axis=json.dumps([0, -1, 0])))
        result = self.run_cli(*self.filters(), "--alignment-mode", "distance-scan",
                              "--scan-start", 1, "--scan-end", 3, "--scan-count", 1)
        self.assertEqual(result.returncode, 0, result.stderr)
        directory = self.directory / "aligned_pairs"
        self.assertEqual(len(list(directory.glob("*.sdf"))), 4)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            metadata = list(csv.DictReader(source))
        self.assertEqual(len(metadata), 12)
        for i in range(4):
            group = metadata[i * 3:(i + 1) * 3]
            self.assertEqual({row["sdf_file"] for row in group}, {f"aligned_pairs/axis-pair-{i+1:06d}.sdf"})
            self.assertEqual([int(row["sdf_record_index"]) for row in group], [0, 1, 2])
            molecules = list(MoleculeSelector._read_sdf(self.directory / group[0]["sdf_file"]))
            self.assertEqual(len(molecules), 3)
            for record, row in zip(molecules, group):
                self.assertEqual(record.GetProp("quick_qm_pair_id"), row["pair_id"])
                self.assertEqual(record.GetProp("quick_qm_distance_angstrom"), row["distance_angstrom"])

    def test_overwrite_pair_files_preserves_unrelated_and_old_files(self):
        args = [*self.filters(), "--alignment-mode", "distance-scan",
                "--scan-start", 1, "--scan-end", 3, "--scan-count", 1]
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        directory = self.directory / "aligned_pairs"
        output = directory / "axis-pair-000001.sdf"
        original = output.read_bytes()
        extra = directory / "axis-pair-000002.sdf"
        extra.write_text("keep old file", encoding="utf-8")
        unrelated = directory / "notes.txt"
        unrelated.write_text("keep", encoding="utf-8")
        result = self.run_cli(*self.filters(), "--distance", 9)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output.read_bytes(), original)
        result = self.run_cli(*self.filters(), "--distance", 9, "--overwrite-output")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(list(MoleculeSelector._read_sdf(output))), 1)
        self.assertEqual(extra.read_text(encoding="utf-8"), "keep old file")
        self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep")
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 1)
        self.assertEqual(float(rows[0]["distance_angstrom"]), 9)

    def test_single_overwrite_yaml_and_cli_override(self):
        config = self.directory / "config.yaml"
        settings = {
            "rigid": dict(molecule_name="rigid", selector=".*", axis_type=".*"),
            "flexible": dict(molecule_name="flexible", selector=".*", axis_type=".*"),
            "output_format": "single", "output": "single.sdf", "overwrite_output": True,
        }
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        sdf = self.directory / "single.sdf"
        before = sdf.read_bytes()
        result = self.run_cli("--config", config, "--no-overwrite-output", "--distance", 9)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sdf.read_bytes(), before)
        result = self.run_cli("--config", config, "--distance", 9)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "single.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 9)
        settings["overwrite_output"] = "true"
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("overwrite_output must be a boolean", result.stderr)

    def test_overwrite_does_not_replace_input_files(self):
        before = self.input.read_bytes()
        result = self.run_cli(*self.filters(), "--output-format", "single",
                              "--output", self.input, "--overwrite-output")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.input.read_bytes(), before)
        before = self.paths[0].read_bytes()
        result = self.run_cli(*self.filters(), "--output-format", "single",
                              "--output", self.paths[0], "--overwrite-output")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("overwrite an input file", result.stderr)
        self.assertEqual(self.paths[0].read_bytes(), before)

    def test_calculation_failure_preserves_existing_output(self):
        result = self.run_cli(*self.filters())
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.directory / "aligned_pairs" / "axis-pair-000001.sdf"
        metadata = self.directory / "aligned_pairs.csv"
        before = output.read_bytes(), metadata.read_bytes()
        previous = Path.cwd()
        try:
            os.chdir(self.directory)
            with patch.object(PairAligner, "align", side_effect=ValueError("test calculation failure")):
                status = aligner.main([str(self.input), *self.filters(), "--overwrite-output"])
        finally:
            os.chdir(previous)
        self.assertEqual(status, 1)
        self.assertEqual((output.read_bytes(), metadata.read_bytes()), before)
        self.assertFalse(list(self.directory.glob(".aligner-*")))

    def test_replacement_failure_restores_old_output(self):
        result = self.run_cli(*self.filters())
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.directory / "aligned_pairs" / "axis-pair-000001.sdf"
        metadata = self.directory / "aligned_pairs.csv"
        before = output.read_bytes(), metadata.read_bytes()
        real_replace = os.replace
        replacements = 0

        def fail_third_replace(source, destination):
            nonlocal replacements
            replacements += 1
            if replacements == 3:
                raise OSError("test replacement failure")
            return real_replace(source, destination)

        previous = Path.cwd()
        try:
            os.chdir(self.directory)
            with patch.object(aligner.os, "replace", side_effect=fail_third_replace):
                status = aligner.main([str(self.input), *self.filters(), "--distance", "8", "--overwrite-output"])
        finally:
            os.chdir(previous)
        self.assertEqual(status, 1)
        self.assertEqual((output.read_bytes(), metadata.read_bytes()), before)
        self.assertFalse(list(self.directory.glob(".aligner-*")))


if __name__ == "__main__":
    unittest.main()
