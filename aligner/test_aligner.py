import csv
import json
import math
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

from aligner import aligner
from aligner.aligner import PairAligner
from aligner.reference_frame import ReferenceSelector
from molecular_selector.matching import MolecularSelector


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "aligner" / "aligner.py"


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
            mol = next(MolecularSelector._read_sdf(path))
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
        pair = next(MolecularSelector._read_sdf(self.directory / row["sdf_file"]))
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        flexible = next(MolecularSelector._read_sdf(self.paths[1]))
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
            "output_format": "single", "output": "all.sdf", "grid_based_scan": {"y": {"distance": 4.5}},
        }), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        molecules = list(MolecularSelector._read_sdf(self.directory / "all.sdf"))
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
        source = ROOT / "examples"
        output = self.directory / "reference"
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "test_out_aligned.csv", "--config",
             str(ROOT / "examples" / "configs" / "aligner.config"), "--output", str(output)],
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
            "grid_based_scan": {"y": {"distance": 1.0}},
        }), encoding="utf-8")
        result = self.run_cli("--config", config, "--rigid", "rigid", ".*", ".*",
                              "--distance", "2.5", "-o", "override")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "override.csv").open(newline="") as source:
            self.assertEqual(float(next(csv.DictReader(source))["distance_angstrom"]), 2.5)

    def test_include_exclude_filter_combinations(self):
        def filt(name, selector=".*", axis=".*"):
            return {"molecule_name": name, "selector": selector, "axis_type": axis}

        a = dict(molecule_name="alpha", selector="oxygen", axis_type="sigma")
        b = dict(molecule_name="beta", selector="oxygen", axis_type="sigma")
        c = dict(molecule_name="beta", selector="nitrogen", axis_type="sigma")
        cases = [
            ({}, [True, True, True]),
            ({"include": [], "exclude": []}, [True, True, True]),
            ({"include": [filt("^alpha$")]}, [True, False, False]),
            ({"exclude": [filt("^beta$", "nitrogen")]}, [True, True, False]),
            ({"include": [filt("alpha"), filt("beta")],
              "exclude": [filt("beta", "oxygen"), filt("alpha")]}, [False, False, True]),
            ({"include": [filt("alpha"), filt(".*")]}, [True, True, True]),
            ({"include": [{"molecule_name": "^beta$"}]}, [False, True, True]),
            ({"include": {"selector": "^oxygen$"}}, [True, True, False]),
            ({"exclude": [{"axis_type": "^sigma$"}]}, [False, False, False]),
            ({"include": [{"molecule_name": "^beta$", "selector": "^oxygen$"}]},
             [False, True, False]),
            ({"exclude": [{"molecule_name": "^beta$", "selector": "^oxygen$"}]},
             [True, False, True]),
            ({"include": [{"selector": "nitrogen"}, {"molecule_name": "^alpha$"}],
              "exclude": [{"molecule_name": "^beta$"}]}, [True, False, False]),
            ({"molecule_name": "^alpha$"}, [True, False, False]),
            ([{"selector": "^nitrogen$"}], [False, False, True]),
        ]
        for config, expected in cases:
            with self.subTest(config=config):
                matcher = aligner.AxisFilter(config)
                self.assertEqual([matcher.matches(row) for row in (a, b, c)], expected)
        for bad in ({"include": None}, {"exclude": ".*"}, {"include": [{}]},
                    {"include": [{"molecule_name": ".*", "unknown": ".*"}]},
                    {"exclude": [{"unknown": ".*"}]}, {"include": [filt("[")]},
                    {"include": [{"selector": None}]}, {"invalid": []}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                aligner.AxisFilter(bad)

    def test_filter_cli_yaml_overrides_and_default_cartesian_product(self):
        result = self.run_cli("--output", "all")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "all.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual({(r["rigid_csv_line"], r["flexible_csv_line"]) for r in rows},
                         {("2", "2"), ("2", "3"), ("3", "2"), ("3", "3")})
        config = self.directory / "filters.yaml"
        config.write_text(yaml.safe_dump({
            "rigid": {"include": [{"selector": "oxygen"}, {"molecule_name": "flexible"}],
                      "exclude": [{"molecule_name": "flexible"}]},
            "flexible": {"exclude": [{"molecule_name": "rigid", "axis_type": "sigma"}]},
        }), encoding="utf-8")
        result = self.run_cli("--config", config, "--output", "filtered")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "filtered.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual([(r["rigid_csv_line"], r["flexible_csv_line"]) for r in rows], [("2", "3")])
        result = self.run_cli("--config", config, "--output", "overridden",
                              "--rigid", ".*", ".*", ".*",
                              "--rigid-exclude", "rigid", ".*", ".*",
                              "--flexible-exclude", "flexible", ".*", ".*")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "overridden.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual([(r["rigid_csv_line"], r["flexible_csv_line"]) for r in rows], [("3", "2")])
        result = self.run_cli("--config", config, "--output", "empty",
                              "--rigid-exclude", ".*", ".*", ".*")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("0 rigid", result.stderr)
        self.assertFalse((self.directory / "empty.csv").exists())

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
                    pair = list(MolecularSelector._read_sdf(self.directory / row["sdf_file"]))[int(row["sdf_record_index"])]
                    self.assertEqual(pair.GetProp("quick_qm_pose_index"), row["pose_index"])
                    self.assertEqual(pair.GetProp("quick_qm_distance_angstrom"), row["distance_angstrom"])
                    rigid = next(MolecularSelector._read_sdf(self.paths[0]))
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
            "output_format": "single", "output": "scan.sdf",
            "grid_based_scan": {"y": {"scan_start": 0, "scan_end": 3, "scan_step": 1,
                                     "include_start": False, "include_end": False}},
        }), encoding="utf-8")
        result = self.run_cli("--config", config, "--scan-count", 2, "--include-start")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "scan.csv").open(newline="") as source:
            metadata = list(csv.DictReader(source))
        molecules = list(MolecularSelector._read_sdf(self.directory / "scan.sdf"))
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
            "grid_based_scan": {"y": {"scan_start": 0, "scan_end": 1,
                                     "scan_count": 4, "include_start": "false"}},
        }
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("include_start must be a boolean", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs.csv").exists())

        settings["grid_based_scan"]["y"]["include_start"] = False
        config.write_text(yaml.safe_dump(settings), encoding="utf-8")
        result = self.run_cli("--config", config, "--scan-step", "0.5")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            self.assertEqual([float(row["distance_angstrom"]) for row in csv.DictReader(source)], [0.5, 1])

    def test_grouped_scan_validation_and_fixed_override(self):
        config = self.directory / "grouped.yaml"
        invalid = (
            {"alignment_mode": "distance-scan"},
            {"x": {"distance": 1}},
            {"grid_based_scan": {"Y": {"distance": 1}}},
            {"grid_based_scan": []},
            {"grid_based_scan": {"x": {"alignment_mode": "distance-scan"}}},
            {"grid_based_scan": {"y": {"scan_start": -1, "scan_end": 1, "scan_count": 1}}},
            {"grid_based_scan": {"y": {"distance": 2, "scan_start": 1, "scan_end": 3, "scan_count": 1}}},
        )
        for value in invalid:
            with self.subTest(value=value):
                config.write_text(yaml.safe_dump(value), encoding="utf-8")
                result = self.run_cli(*self.filters(), "--config", config)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self.directory / "aligned_pairs.csv").exists())
        config.write_text(yaml.safe_dump({"grid_based_scan": {
            "y": {"scan_start": 1, "scan_end": 3, "scan_count": 1},
            "x": {"scan_start": -1, "scan_end": 1, "scan_step": 1},
        }}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--distance", 4, "--x-distance", -2)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 1)
        self.assertEqual(float(rows[0]["distance_angstrom"]), 4)
        self.assertEqual(float(rows[0]["x_offset_angstrom"]), -2)

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
            molecules = list(MolecularSelector._read_sdf(self.directory / group[0]["sdf_file"]))
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
        self.assertEqual(len(list(MolecularSelector._read_sdf(output))), 1)
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

    def test_reference_priority_and_degenerate_fallback(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conformer.SetAtomPosition(i, (float(i), float(i * i), 1.0))
        mol.AddConformer(conformer)
        row = aligner.AxisRow(7, {}, Path("test.sdf"), 0, 2,
                              np.array([2.0, 4.0, 1.0]), np.array([0.0, 0.0, 1.0]))
        candidates = sorted({other.GetIdx() for neighbor in mol.GetAtomWithIdx(2).GetNeighbors()
                             for other in neighbor.GetNeighbors() if other.GetIdx() != 2})
        a = ReferenceSelector().select(row, mol)
        b = ReferenceSelector({"priority": ["B", "A"]}).select(row, mol)
        self.assertEqual((a.method, a.atom_index), ("A", candidates[0]))
        self.assertEqual((b.method, b.atom_index), ("B", candidates[-1]))
        c = ReferenceSelector({"priority": ["A", "B", "C"],
                               "enabled": {"A": False, "B": False}}).select(row, mol)
        self.assertEqual((c.method, c.atom_index), ("C", None))
        np.testing.assert_allclose(np.dot(a.x, row.axis), 0, atol=1e-12)
        np.testing.assert_allclose(a.z, np.cross(a.x, row.axis), atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(a.x), 1)
        np.testing.assert_allclose(np.linalg.norm(a.z), 1)
        position = mol.GetConformer().GetAtomPosition(candidates[0])
        row_parallel = aligner.AxisRow(8, {}, Path("test.sdf"), 0, 2,
                                       np.array([position.x, position.y, position.z - 1]),
                                       np.array([0.0, 0.0, 1.0]))
        self.assertEqual(ReferenceSelector({"priority": ["A", "B"]}).select(row_parallel, mol).method, "B")
        center = aligner.AxisRow(9, {}, Path("test.sdf"), 0, None,
                                 np.zeros(3), np.array([0.0, 1.0, 0.0]))
        with self.assertRaisesRegex(ValueError, "CSV line 9.*no usable reference"):
            ReferenceSelector().select(center, mol)
        with self.assertRaisesRegex(ValueError, "no usable reference"):
            ReferenceSelector({"enabled": {"A": False, "B": False, "C": False}}).select(row, mol)

    def test_lateral_cartesian_grid_geometry_and_provenance(self):
        result = self.run_cli(*self.filters(), "--distance", 2,
                              "--x-alignment-mode", "distance-scan", "--x-scan-start", -1,
                              "--x-scan-end", 1, "--x-scan-count", 1,
                              "--z-alignment-mode", "distance-scan", "--z-scan-start", -2,
                              "--z-scan-end", 2, "--z-scan-count", 1)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 9)
        self.assertEqual([(float(r["x_offset_angstrom"]), float(r["z_offset_angstrom"])) for r in rows],
                         [(x, z) for x in (-1, 0, 1) for z in (-2, 0, 2)])
        molecules = list(MolecularSelector._read_sdf(self.directory / rows[0]["sdf_file"]))
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        for row, mol in zip(rows, molecules):
            xaxis, yaxis, zaxis = (np.array(json.loads(row[key])) for key in ("x_axis", "rigid_axis", "z_axis"))
            np.testing.assert_allclose(np.cross(xaxis, yaxis), zaxis, atol=1e-12)
            center = np.array(json.loads(row["rigid_center"]))
            target = center + 2 * yaxis + float(row["x_offset_angstrom"]) * xaxis + float(row["z_offset_angstrom"]) * zaxis
            np.testing.assert_allclose(json.loads(row["flexible_aligned_center"]), target, atol=1e-12)
            point = mol.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + 1)
            np.testing.assert_allclose(point, target, atol=0.00006)
            np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]), -yaxis, atol=1e-12)
            self.assertEqual(mol.GetProp("quick_qm_reference_method"), row["reference_method"])
            self.assertEqual(mol.GetProp("quick_qm_x_offset_angstrom"), row["x_offset_angstrom"])

    def test_toward_center_rotation_preserves_centers_and_internal_geometry(self):
        config = self.directory / "orientation.yaml"
        config.write_text(yaml.safe_dump({"flexible_orientation": "toward-center"}), encoding="utf-8")
        args = [*self.filters(), "--distance", 2, "--x-alignment-mode", "distance-scan",
                "--x-scan-start", -1, "--x-scan-end", 1, "--x-scan-count", 1,
                "--z-distance", 0.5]
        result = self.run_cli(*args, "--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        molecules = list(MolecularSelector._read_sdf(self.directory / rows[0]["sdf_file"]))
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        flexible = next(MolecularSelector._read_sdf(self.paths[1]))
        self.assertEqual(len(rows), 3)
        for row, molecule in zip(rows, molecules):
            rcenter = np.array(json.loads(row["rigid_center"]))
            fcenter = np.array(json.loads(row["flexible_aligned_center"]))
            direction = rcenter - fcenter
            direction /= np.linalg.norm(direction)
            matrix = np.array(json.loads(row["rotation"]))
            translation = np.array(json.loads(row["translation"]))
            np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]), direction, atol=1e-12)
            np.testing.assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(matrix), 1)
            self.assertEqual(row["flexible_orientation"], "toward-center")
            self.assertEqual(molecule.GetProp("quick_qm_flexible_orientation"), "toward-center")
            for i in range(flexible.GetNumAtoms()):
                original = np.array(flexible.GetConformer().GetAtomPosition(i))
                transformed = np.array(molecule.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + i))
                np.testing.assert_allclose(transformed, matrix @ original + translation, atol=0.00006)
            np.testing.assert_allclose(molecule.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + 1),
                                       fcenter, atol=0.00006)
        result = self.run_cli(*args, "--config", config, "--flexible-orientation", "stable", "-o", "stable")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "stable.csv").open(newline="") as source:
            stable_rows = list(csv.DictReader(source))
        for stable, changed in zip(stable_rows, rows):
            np.testing.assert_allclose(json.loads(stable["flexible_aligned_center"]),
                                       json.loads(changed["flexible_aligned_center"]), atol=1e-12)
            np.testing.assert_allclose(json.loads(stable["flexible_aligned_axis"]),
                                       -np.array(json.loads(stable["rigid_axis"])), atol=1e-12)
            self.assertEqual(stable["flexible_orientation"], "stable")

    def test_toward_center_zero_and_on_axis_cases_and_invalid_mode(self):
        for distance in (0, 3):
            name = f"zero-{distance}"
            result = self.run_cli(*self.filters(), "--distance", distance, "-o", name,
                                  "--flexible-orientation", "toward-center")
            self.assertEqual(result.returncode, 0, result.stderr)
            with (self.directory / f"{name}.csv").open(newline="") as source:
                row = next(csv.DictReader(source))
            np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]),
                                       -np.array(json.loads(row["rigid_axis"])), atol=1e-12)
        config = self.directory / "bad-orientation.yaml"
        config.write_text(yaml.safe_dump({"flexible_orientation": "invalid"}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("flexible_orientation", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs.csv").exists())

    def test_polar_cardinal_positions_and_orientation(self):
        config = self.directory / "polar.yaml"
        config.write_text(yaml.safe_dump({"polar_coordinate_based_scan": {
            "radius": {"distance": 2},
            "phi": {"scan_start": 0, "scan_end": 360, "scan_step": 90},
            "theta": {"distance": 90},
        }, "flexible_orientation": "toward-center"}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("merged 1 polar positions", result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 4)
        mols = list(MolecularSelector._read_sdf(self.directory / rows[0]["sdf_file"]))
        rigid = next(MolecularSelector._read_sdf(self.paths[0]))
        for row, molecule, expected in zip(rows, mols, [(2, 0), (0, 2), (-2, 0), (0, -2)]):
            self.assertEqual(row["coordinate_system"], "polar")
            self.assertEqual(row["polar_phi_unit"], "degrees")
            self.assertAlmostEqual(float(row["polar_radius_angstrom"]), 2)
            self.assertAlmostEqual(float(row["x_offset_angstrom"]), expected[0], places=10)
            self.assertAlmostEqual(float(row["z_offset_angstrom"]), expected[1], places=10)
            self.assertAlmostEqual(float(row["distance_angstrom"]), 0, places=10)
            center = np.array(json.loads(row["rigid_center"]))
            target = np.array(json.loads(row["flexible_aligned_center"]))
            direction = center - target
            np.testing.assert_allclose(json.loads(row["flexible_aligned_axis"]), direction / np.linalg.norm(direction), atol=1e-12)
            np.testing.assert_allclose(molecule.GetConformer().GetAtomPosition(rigid.GetNumAtoms() + 1),
                                       target, atol=0.00006)
            self.assertEqual(molecule.GetProp("quick_qm_polar_phi"), row["polar_phi"])

    def test_polar_radians_poles_zero_radius_and_single_output(self):
        config = self.directory / "polar.yaml"
        config.write_text(yaml.safe_dump({"polar_coordinate_based_scan": {
            "radius": {"scan_start": 0, "scan_end": 1, "scan_step": 1},
            "phi": {"scan_start": 0, "scan_end": 180, "scan_step": 90},
            "phi_unit": "degrees",
            "theta": {"scan_start": 0, "scan_end": math.pi, "scan_count": 1},
            "theta_unit": "radians",
        }}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--output-format", "single",
                              "--output", "polar.sdf", "--polar-merge-tolerance", 0.01)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("merged 12 polar positions", result.stderr)
        with (self.directory / "polar.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 6)
        self.assertEqual([int(row["sdf_record_index"]) for row in rows], list(range(6)))
        self.assertEqual(len(list(MolecularSelector._read_sdf(self.directory / "polar.sdf"))), 6)
        self.assertEqual({row["polar_theta_unit"] for row in rows}, {"radians"})
        self.assertEqual(sum(float(row["polar_radius_angstrom"]) == 0 for row in rows), 1)

    def test_polar_cli_scan_and_merge_tolerance(self):
        args = [*self.filters(), "--coordinate-system", "polar", "--polar-radius", 2,
                "--polar-phi-scan-start", 0, "--polar-phi-scan-end", 2,
                "--polar-phi-scan-step", 1, "--polar-theta", 90]
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual([float(r["polar_phi"]) for r in rows], [0, 2])
        self.assertIn("merged 1 polar positions", result.stderr)
        result = self.run_cli(*args, "--polar-merge-tolerance", 0.01, "-o", "fine")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "fine.csv").open(newline="") as source:
            self.assertEqual(len(list(csv.DictReader(source))), 3)

    def test_polar_invalid_options_fail_without_output(self):
        cases = [
            {"polar_coordinate_based_scan": {"phi_unit": "grads"}},
            {"polar_coordinate_based_scan": {"phi": {"distance": 361}}},
            {"polar_coordinate_based_scan": {"theta": {"scan_start": 0, "scan_end": 181, "scan_count": 1, "include_end": False}}},
            {"polar_coordinate_based_scan": {"radius": {"distance": -1}}},
            {"polar_coordinate_based_scan": {"merge_tolerance": 0}},
            {"distance_scan": {"y": {"distance": 2}}},
        ]
        config = self.directory / "bad-polar.yaml"
        for value in cases:
            with self.subTest(value=value):
                config.write_text(yaml.safe_dump(value), encoding="utf-8")
                result = self.run_cli(*self.filters(), "--config", config)
                self.assertNotEqual(result.returncode, 0, result.stderr)
                self.assertFalse((self.directory / "aligned_pairs.csv").exists())
        result = self.run_cli(*self.filters(), "--coordinate-system", "polar", "--distance", 3)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("conflicts", result.stderr)

    def test_multiple_mixed_scans_merge_and_record_provenance(self):
        config = self.directory / "multiple.yaml"
        config.write_text(yaml.safe_dump({"scans": [
            {"grid_based_scan": {"y": {"scan_start": 2, "scan_end": 2.04, "scan_step": 0.02}}},
            {"polar_coordinate_based_scan": {"radius": {"distance": 2}, "theta": {"distance": 0}}},
            {"grid_based_scan": {"y": {"distance": 2.045}}},
            {"grid_based_scan": {"y": {"distance": 3}}},
        ]}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--output", "multiple")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "multiple.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual([row["scan_index"] for row in rows], ["1", "1", "1", "4"])
        self.assertEqual([row["pose_index"] for row in rows], ["1", "2", "3", "4"])
        self.assertEqual([float(row["distance_angstrom"]) for row in rows], [2, 2.02, 2.04, 3])
        self.assertIn("merged 1 positions from earlier scans", result.stderr)
        molecules = list(MolecularSelector._read_sdf(self.directory / rows[0]["sdf_file"]))
        self.assertEqual([m.GetProp("quick_qm_scan_index") for m in molecules], ["1", "1", "1", "4"])

    def test_multiple_scan_tolerance_and_cli_overrides(self):
        config = self.directory / "multiple.yaml"
        config.write_text(yaml.safe_dump({
            "scan_merge_tolerance": 0.1,
            "scans": [
                {"grid_based_scan": {"y": {"distance": 2}, "x": {"distance": 0}}},
                {"grid_based_scan": {"y": {"distance": 3}, "x": {"distance": 0.06}}},
                {"polar_coordinate_based_scan": {"radius": {"distance": 3},
                                                   "phi": {"distance": 0}, "theta": {"distance": 90}}},
                {"polar_coordinate_based_scan": {"radius": {"distance": 3},
                                                   "phi": {"distance": 90}, "theta": {"distance": 90}}},
            ],
        }), encoding="utf-8")
        args = [*self.filters(), "--config", config, "--distance", 4, "--polar-radius", 5]
        result = self.run_cli(*args, "--output", "coarse")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "coarse.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual([r["scan_index"] for r in rows], ["1", "3", "4"])
        self.assertEqual(float(rows[0]["distance_angstrom"]), 4)
        self.assertEqual([float(r["polar_radius_angstrom"]) for r in rows[1:]], [5, 5])
        result = self.run_cli(*args, "--scan-merge-tolerance", 0.01, "--output", "fine")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "fine.csv").open(newline="") as source:
            self.assertEqual([r["scan_index"] for r in csv.DictReader(source)], ["1", "2", "3", "4"])

    def test_multiple_scans_top_level_lists_and_validation(self):
        config = self.directory / "multiple.yaml"
        config.write_text(yaml.safe_dump({
            "grid_based_scan": [{"y": {"distance": 1}}, {"y": {"distance": 2}}],
            "polar_coordinate_based_scan": [{"radius": {"distance": 3}, "theta": {"distance": 0}}],
        }, sort_keys=False), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            self.assertEqual([r["scan_index"] for r in csv.DictReader(source)], ["1", "2", "3"])
        for invalid in (
            {"scans": []}, {"scans": [{}]}, {"scans": [{"grid_based_scan": {}, "polar_coordinate_based_scan": {}}]},
            {"scans": [{"grid_based_scan": []}]}, {"grid_based_scan": []},
            {"scans": [{"grid_based_scan": {}}], "grid_based_scan": {}},
            {"scan_merge_tolerance": 0}, {"scan_merge_tolerance": "nan"},
        ):
            with self.subTest(invalid=invalid):
                config.write_text(yaml.safe_dump(invalid), encoding="utf-8")
                result = self.run_cli(*self.filters(), "--config", config, "--output", "invalid")
                self.assertNotEqual(result.returncode, 0, result.stderr)
                self.assertFalse((self.directory / "invalid.csv").exists())

    def test_polar_cli_radians_and_sampling_override_per_pair(self):
        rows = self.write_rows()
        with self.input.open("a", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writerow(dict(rows[0], axis=json.dumps([-1, 0, 0])))
            writer.writerow(dict(rows[1], axis=json.dumps([0, -1, 0])))
        config = self.directory / "polar.yaml"
        config.write_text(yaml.safe_dump({"polar_coordinate_based_scan": {
            "radius": {"distance": 2},
            "phi": {"scan_start": 0, "scan_end": math.pi, "scan_count": 1},
            "theta": {"distance": 90},
            "phi_unit": "radians",
        }}), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--polar-phi-scan-step",
                              str(math.pi), "--no-polar-phi-include-end")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 4)
        self.assertEqual(len({row["axis_pair_id"] for row in rows}), 4)
        self.assertTrue(all(row["polar_phi_unit"] == "radians" and float(row["polar_phi"]) == 0 for row in rows))

    def test_lateral_yaml_priority_cli_overrides_and_3d_scan(self):
        config = self.directory / "lateral.yaml"
        config.write_text(yaml.safe_dump({
            "grid_based_scan": {"x": {"scan_start": -1, "scan_end": 1, "scan_step": 1},
                              "z": {"distance": -0.5}},
            "reference_axis": {"priority": ["B", "A", "C"], "enabled": {"B": True}},
        }), encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--alignment-mode", "distance-scan",
                              "--scan-start", 1, "--scan-end", 2, "--scan-step", 1,
                              "--x-scan-count", 1, "--no-x-include-start", "--z-distance", 0.5)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "aligned_pairs.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["reference_method"] for row in rows}, {"B"})
        self.assertEqual([(float(r["distance_angstrom"]), float(r["x_offset_angstrom"]),
                           float(r["z_offset_angstrom"])) for r in rows],
                         [(1, 0, 0.5), (1, 1, 0.5), (2, 0, 0.5), (2, 1, 0.5)])

    def test_invalid_reference_and_lateral_scan_do_not_write_output(self):
        config = self.directory / "invalid.yaml"
        for value in ({"priority": ["A", "A"]}, {"priority": ["X"]}, {"enabled": {"A": "yes"}}):
            config.write_text(yaml.safe_dump({"reference_axis": value}), encoding="utf-8")
            result = self.run_cli(*self.filters(), "--config", config, "--x-distance", 1)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("reference_axis", result.stderr)
            self.assertFalse((self.directory / "aligned_pairs.csv").exists())
        config.write_text(yaml.safe_dump({"reference_axis": {"enabled": {"A": False, "B": False, "C": False}}}),
                          encoding="utf-8")
        result = self.run_cli(*self.filters(), "--config", config, "--x-distance", 1)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CSV line 2: no usable reference", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs").exists())
        result = self.run_cli(*self.filters(), "--x-alignment-mode", "distance-scan",
                              "--x-scan-start", -1, "--x-scan-end", 1, "--x-scan-step", 0.3)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("evenly divisible", result.stderr)
        self.assertFalse((self.directory / "aligned_pairs.csv").exists())

    def test_ring_origin_fallback_and_missing_reference_preserves_results(self):
        rows = self.write_rows()
        rows[0].update(center_atom_index="", axis_type="aromatic_normal",
                       center=json.dumps([1, 0, 0]), axis=json.dumps([0, 1, 0]))
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        filters = ["--rigid", "rigid", ".*", "aromatic_normal",
                   "--flexible", "flexible", ".*", "sigma"]
        result = self.run_cli(*filters, "--x-distance", 1, "--output-format", "single", "-o", "ring.sdf")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / "ring.csv").open(newline="") as source:
            row = next(csv.DictReader(source))
        self.assertEqual(row["reference_method"], "C")
        self.assertEqual(row["reference_atom_index"], "")
        np.testing.assert_allclose(json.loads(row["x_axis"]), [-1, 0, 0])
        np.testing.assert_allclose(json.loads(row["z_axis"]), [0, 0, -1])
        before = (self.directory / "ring.sdf").read_bytes(), (self.directory / "ring.csv").read_bytes()
        rows[0]["center"] = json.dumps([0, 0, 0])
        with self.input.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        result = self.run_cli(*filters, "--x-distance", 1, "--output-format", "single",
                              "-o", "ring.sdf", "--overwrite-output")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no usable reference", result.stderr)
        self.assertEqual(((self.directory / "ring.sdf").read_bytes(),
                          (self.directory / "ring.csv").read_bytes()), before)
        result = self.run_cli(*filters, "--output-format", "single", "-o", "y-only.sdf")
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
