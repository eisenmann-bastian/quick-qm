import csv
import logging
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from molecular_selector import MolecularSelector

ROOT = Path(__file__).resolve().parent.parent


class CrossRepositoryTests(unittest.TestCase):
    def test_cli_logging_does_not_configure_selector_namespace(self):
        from alignment import aligner

        selector_logger = logging.getLogger("molecular_selector")
        previous = (selector_logger.level, selector_logger.handlers[:], selector_logger.propagate)
        with tempfile.TemporaryDirectory() as directory:
            status = aligner.main([str(Path(directory) / "missing.csv")])
        self.assertEqual(status, 1)
        self.assertEqual(
            (selector_logger.level, selector_logger.handlers, selector_logger.propagate),
            previous,
        )

    def test_all_cli_errors_use_shared_project_namespace(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = str(Path(directory) / "missing.csv")
            commands = [
                [sys.executable, str(ROOT / "alignment" / "axes.py"), missing],
                [sys.executable, str(ROOT / "alignment" / "aligner.py"), missing],
                [sys.executable, str(ROOT / "selection" / "selector.py"), "--csv", missing],
            ]
            for command in commands:
                with self.subTest(command=command):
                    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 1, result.stderr)
                    self.assertIn("ERROR molecular_aligner.", result.stderr)
                    self.assertEqual(result.stdout, "")

    def test_selection_axes_and_alignment(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            selections = output / "selections.csv"
            axes = output / "axes.csv"
            structures = output / "pairs.sdf"
            commands = [
                [sys.executable, "-m", "molecular_selector", "--csv",
                 str(ROOT / "selection" / "reference_molecules.csv"), "-o", str(selections)],
                [sys.executable, str(ROOT / "alignment" / "axes.py"), str(selections), "-o", str(axes)],
                [sys.executable, str(ROOT / "alignment" / "aligner.py"), str(axes),
                 "--rigid", "^241$", "^aromatic ring$", "^aromatic_normal$",
                 "--flexible", "^878$", "^thiol$", "^hydrogen_sigma$",
                 "--output-format", "single", "--output", str(structures)],
            ]
            for command in commands:
                result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            with structures.with_suffix(".csv").open(newline="", encoding="utf-8") as source:
                rows = list(csv.DictReader(source))
            molecules = list(MolecularSelector.read_sdf(structures))
            self.assertGreater(len(rows), 0)
            self.assertEqual(len(molecules), len(rows))
            self.assertTrue(all(row["rigid_source"] and row["flexible_source"] for row in rows))
