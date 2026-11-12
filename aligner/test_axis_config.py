import csv
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

from aligner.axes import SDF_AXES_PROPERTY
from aligner.axis_config import AXIS_NAMES, DEFAULT_CONFIG, AxisConfig
from molecular_selector.matching import MolecularSelector


ROOT = Path(__file__).resolve().parent.parent


class AxisConfigurationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.config = self.path / "axes.config"
        with DEFAULT_CONFIG.open(encoding="utf-8") as source:
            self.settings = yaml.safe_load(source)

    def save(self):
        self.config.write_text(yaml.safe_dump(self.settings), encoding="utf-8")

    def command(self, *args):
        return subprocess.run(
            [sys.executable, str(ROOT / "aligner" / "selector.py"), *map(str, args)],
            cwd=ROOT, capture_output=True, text=True,
        )

    def test_default_config_is_complete_and_enabled(self):
        config = AxisConfig.load()
        self.assertEqual(set(config.definitions), set(AXIS_NAMES.values()))
        for name in AXIS_NAMES.values():
            self.assertTrue(config.enabled(name, 1))
            self.assertTrue(config.enabled(name, 2))
            self.assertTrue(config.enabled(name, 99))

    def test_aromatic_direction_and_sdf_marker_are_filtered(self):
        self.settings["axes"]["aromatic_ring_normal"]["directions"][2] = False
        self.save()
        ring = ROOT / "examples" / "reference_molecules" / "241.sdf"
        common = [ring, "--name", "aromatic ring", "--axes-config", self.config]
        result = self.command(*common, "--format", "alignment")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["axis_type"], "aromatic_normal")
        self.assertEqual(rows[0]["axis_name"], "aromatic_ring_normal_1")
        output = self.path / "ring.sdf"
        result = self.command(*common, "--format", "sdf", "-o", output)
        self.assertEqual(result.returncode, 0, result.stderr)
        mol = next(MolecularSelector._read_sdf(output))
        axes = json.loads(mol.GetProp(SDF_AXES_PROPERTY))
        self.assertEqual(len(axes), 1)
        self.assertEqual(axes[0]["sources"][0]["axis_name"], "aromatic_ring_normal_1")

    def test_thiol_direction_and_other_axes_remain(self):
        self.settings["axes"]["chalcogen_bond_extension"]["directions"][2] = False
        self.settings["axes"]["thiol_plane_normal"]["directions"][1] = False
        self.save()
        thiol = ROOT / "examples" / "reference_molecules" / "878.sdf"
        result = self.command(thiol, "--name", "thiol", "--format", "alignment",
                              "--axes-config", self.config)
        self.assertEqual(result.returncode, 0, result.stderr)
        names = [row["axis_name"] for row in csv.DictReader(io.StringIO(result.stdout))]
        self.assertEqual(names, ["chalcogen_bond_extension_1", "thiol_plane_normal_2",
                                 "thiol_outward_bisector_1", "hydrogen_donor_extension_1"])

    def test_group_switch_disables_both_normals(self):
        self.settings["axes"]["thiol_plane_normal"]["enabled"] = False
        self.save()
        thiol = ROOT / "examples" / "reference_molecules" / "878.sdf"
        result = self.command(thiol, "--name", "thiol", "--format", "alignment",
                              "--axes-config", self.config)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        self.assertEqual(len(rows), 4)
        self.assertNotIn("thiol_plane_normal", [r["axis_name"].rsplit("_", 1)[0] for r in rows])

    def test_selection_output_does_not_apply_axis_filter(self):
        thiol = ROOT / "examples" / "reference_molecules" / "878.sdf"
        result = self.command(thiol, "--name", "thiol")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(list(csv.DictReader(io.StringIO(result.stdout)))), 1)
        result = self.command(thiol, "--name", "thiol", "--axes-config", self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires --format", result.stderr)

    def test_disabled_all_preserves_unmatched_record(self):
        for entry in self.settings["axes"].values():
            entry["enabled"] = False
        self.save()
        ring = ROOT / "examples" / "reference_molecules" / "241.sdf"
        output = self.path / "none.sdf"
        result = self.command(ring, "--name", "aromatic ring", "--format", "sdf",
                              "--axes-config", self.config, "-o", output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("no supported alignment axes", result.stderr)
        result_mol = next(MolecularSelector._read_sdf(output))
        original = next(MolecularSelector._read_sdf(ring))
        self.assertEqual(result_mol.GetNumAtoms(), original.GetNumAtoms())
        self.assertEqual(json.loads(result_mol.GetProp(SDF_AXES_PROPERTY)), [])

    def test_existing_alignment_cli_uses_same_config(self):
        self.settings["axes"]["aromatic_ring_normal"]["directions"][1] = False
        self.save()
        ring = ROOT / "examples" / "reference_molecules" / "241.sdf"
        selected = self.command(ring, "--name", "aromatic ring")
        self.assertEqual(selected.returncode, 0, selected.stderr)
        input_path = self.path / "selection.csv"
        input_path.write_text(selected.stdout, encoding="utf-8")
        result = subprocess.run(
            [sys.executable, str(ROOT / "aligner" / "alignment.py"), str(input_path),
             "--axes-config", str(self.config)],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        self.assertEqual([row["axis_name"] for row in rows], ["aromatic_ring_normal_2"])

    def test_bad_config_fails_before_output(self):
        ring = ROOT / "examples" / "reference_molecules" / "241.sdf"
        output = self.path / "failure.sdf"
        for content in (
            "axes: [",
            "axes: {}\n",
            "axes:\n  42: {enabled: true}\n",
            "axes:\n  x: {enabled: true}\n  x: {enabled: false}\n",
        ):
            with self.subTest(content=content):
                self.config.write_text(content, encoding="utf-8")
                result = self.command(ring, "--name", "aromatic ring", "--format", "sdf",
                                      "--axes-config", self.config, "-o", output)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())
        self.settings["axes"]["aromatic_ring_normal"]["directions"][2] = "false"
        self.save()
        with self.assertRaisesRegex(ValueError, "directions"):
            AxisConfig.load(self.config)
        self.settings["axes"]["aromatic_ring_normal"]["directions"] = {0: False}
        self.save()
        with self.assertRaisesRegex(ValueError, "directions"):
            AxisConfig.load(self.config)
        self.settings["axes"]["aromatic_ring_normal"]["directions"] = {1: False}
        self.settings["axes"]["not_an_axis"] = {"enabled": True}
        self.save()
        with self.assertRaisesRegex(ValueError, "unknown"):
            AxisConfig.load(self.config)


if __name__ == "__main__":
    unittest.main()
