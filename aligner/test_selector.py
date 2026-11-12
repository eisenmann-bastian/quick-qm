import contextlib
import csv
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem

from aligner import selector
from molecular_selector.batch import BatchSelector
from molecular_selector.definitions import NamedSelector, SelectorRegistry
from molecular_selector.matching import MolecularSelector
from molecular_selector.output import CsvResultWriter


class SelectorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "molecules.sdf"

    def write_molecules(self, smiles):
        writer = Chem.SDWriter(str(self.path))
        for index, value in enumerate(smiles):
            molecule = Chem.MolFromSmiles(value)
            molecule.SetProp("_Name", f"molecule-{index}")
            writer.write(molecule)
        writer.close()

    def test_components_share_custom_registry(self):
        self.write_molecules(["CSC"])
        registry = SelectorRegistry({"sulfur": NamedSelector("[S:1]", ("sulphur",))})
        matcher = MolecularSelector(registry)
        input_file = self.path.parent / "custom.csv"
        input_file.write_text("file,sulphur\nmolecules.sdf,x\n", encoding="utf-8")
        result = BatchSelector(matcher).select(input_file)
        output = io.StringIO()
        CsvResultWriter().write(result, output)
        rows = list(csv.DictReader(io.StringIO(output.getvalue())))
        self.assertEqual(rows[0]["selector"], "sulfur")
        self.assertEqual(json.loads(rows[0]["site_atom_indices"]), [1])
        with self.assertRaisesRegex(ValueError, "unknown selector"):
            matcher.select(self.path, name="thioether")

    def test_thioether_highlights_only_sulfur(self):
        self.write_molecules(["CSC", "CSCCSC", "CC"])
        result = selector.select_sdf(self.path, name="ThioEther")
        self.assertEqual(result["atom_index_base"], 0)
        self.assertEqual(
            result["molecules"][0]["matches"],
            [{"atom_indices": [0, 1, 2], "site_atom_indices": [1]}],
        )
        self.assertEqual(
            [match["site_atom_indices"] for match in result["molecules"][1]["matches"]],
            [[1], [4]],
        )
        self.assertEqual(result["molecules"][2]["matches"], [])
        self.assertEqual(selector.select_sdf(self.path, name="sulfide")["selector"]["name"], "thioether")

    def test_custom_smarts_highlights_full_match(self):
        self.write_molecules(["CC(=O)N"])
        result = selector.select_sdf(self.path, smarts="[C](=[O])[N]")
        match = result["molecules"][0]["matches"][0]
        self.assertEqual(match["atom_indices"], [1, 2, 3])
        self.assertEqual(match["site_atom_indices"], [1, 2, 3])

    def test_hydrogen_bond_sites_exclude_amide_nitrogen_as_acceptor(self):
        self.write_molecules(["CC(=O)N", "C[NH3+]", "C=O"])
        donor = selector.select_sdf(self.path, name="hydrogen-bond donor")
        acceptor = selector.select_sdf(self.path, name="hydrogen bond acceptor")
        self.assertEqual(donor["molecules"][0]["matches"][0]["site_atom_indices"], [3])
        self.assertEqual(acceptor["molecules"][0]["matches"][0]["site_atom_indices"], [2])
        self.assertEqual(acceptor["molecules"][1]["matches"], [])
        self.assertEqual(acceptor["molecules"][2]["matches"][0]["site_atom_indices"], [1])

    def test_all_named_smarts_compile_and_are_exposed(self):
        for name, definition in selector.SELECTORS.items():
            with self.subTest(name=name):
                if definition.smarts is not None:
                    self.assertIsNotNone(Chem.MolFromSmarts(definition.smarts))
                self.assertEqual(selector._named_selector(name)[0], name)

    def test_interaction_site_selectors(self):
        examples = {
            "halogen bond donor": ("CCl", "[Cl-]", 1),
            "halogen bond acceptor": ("C=O", "C[NH3+]", 1),
            "positive charge": ("C[NH3+]", "CN", 1),
            "negative charge": ("CC(=O)[O-]", "CC(=O)O", 3),
            "chalcogen bond donor": ("CSC", "CS", 1),
            "chalcogen bond acceptor": ("C=O", "C[NH3+]", 1),
            "metal center": ("[Zn+2]", "C[NH3+]", 0),
            "metal coordination ligand": ("C=O", "C[NH3+]", 1),
        }
        for name, (positive, negative, site) in examples.items():
            with self.subTest(name=name):
                self.write_molecules([positive, negative])
                records = selector.select_sdf(self.path, name=name)["molecules"]
                self.assertEqual(records[0]["matches"][0]["site_atom_indices"], [site])
                self.assertEqual(records[1]["matches"], [])

    def test_halogen_and_chalcogen_donor_scope(self):
        self.write_molecules(["CCl", "CBr", "CI", "CF", "[Cl-]", "CSC", "C[Se]C", "CS", "CS(=O)C"])
        halogens = selector.select_sdf(self.path, name="halogen-bond donor")["molecules"]
        self.assertEqual([bool(record["matches"]) for record in halogens[:5]], [True, True, True, False, False])
        donors = selector.select_sdf(self.path, name="chalcogen donor")["molecules"]
        self.assertEqual([bool(record["matches"]) for record in donors[5:]], [True, True, False, False])

    def test_aromatic_ring_sites_include_heterocyclic_and_fused_rings(self):
        self.write_molecules(["c1ccccc1", "c1ccncc1", "c1ccc2ccccc2c1", "C1CCCCC1", "c1ccccc1-c2ccccc2"])
        result = selector.select_sdf(self.path, name="pi-pi")
        self.assertEqual(result["selector"]["name"], "aromatic ring")
        self.assertIsNone(result["selector"]["smarts"])
        rings = [record["matches"] for record in result["molecules"]]
        self.assertEqual([len(x) for x in rings], [1, 1, 2, 0, 2])
        for molecule in (rings[0], rings[1], rings[2], rings[4]):
            for match in molecule:
                self.assertEqual(len(match["site_atom_indices"]), 6)
                self.assertEqual(match["site_atom_indices"], match["atom_indices"])
        self.assertEqual(len(set(rings[2][0]["site_atom_indices"]) & set(rings[2][1]["site_atom_indices"])), 2)

    def test_functional_group_positive_and_negative_examples(self):
        examples = {
            "thioether": ("CSC", "CS"),
            "thiol": ("CS", "CSC"),
            "carbonyl": ("CC=O", "CCO"),
            "alcohol": ("CCO", "Oc1ccccc1"),
            "phenol": ("Oc1ccccc1", "CCO"),
            "ether": ("COC", "CO"),
            "amine": ("CN", "CNC=O"),
            "amide": ("CNC=O", "CN"),
            "carboxylic acid": ("CC(=O)O", "CC(=O)OC"),
            "ester": ("CC(=O)OC", "CC(=O)O"),
            "aldehyde": ("CC=O", "CC(=O)C"),
            "ketone": ("CC(=O)C", "CC=O"),
        }
        for name, (positive, negative) in examples.items():
            with self.subTest(name=name):
                self.write_molecules([positive, negative])
                records = selector.select_sdf(self.path, name=name)["molecules"]
                self.assertTrue(records[0]["matches"])
                self.assertFalse(records[1]["matches"])

    def test_carboxylate_selects_both_oxygens_in_both_resonance_forms(self):
        self.write_molecules([
            "CC(=O)[O-]", "CC([O-])=O", "O=C[O-]",
            "CC(=O)O", "CC(=O)OC", "C(=O)([O-])OC", "[O-]C(=O)[O-]",
            "O=C([O-])C(=O)[O-]",
        ])
        records = selector.select_sdf(self.path, name="carboxylate")["molecules"]
        for record in records[:3]:
            self.assertEqual(len(record["matches"]), 1)
            self.assertEqual(len(record["matches"][0]["site_atom_indices"]), 2)
        self.assertEqual(records[0]["matches"][0]["site_atom_indices"], [2, 3])
        self.assertEqual(set(records[1]["matches"][0]["site_atom_indices"]), {2, 3})
        for record in records[3:7]:
            self.assertEqual(record["matches"], [])
        self.assertEqual(len(records[7]["matches"]), 2)

    def test_custom_smarts_returns_more_than_default_match_limit(self):
        self.write_molecules(["C." * 1004 + "C"])
        result = selector.select_sdf(self.path, smarts="[#6]")
        self.assertEqual(len(result["molecules"][0]["matches"]), 1005)

    def test_explicit_hydrogens_keep_original_atom_indices(self):
        writer = Chem.SDWriter(str(self.path))
        writer.write(Chem.AddHs(Chem.MolFromSmiles("CSC")))
        writer.close()
        result = selector.select_sdf(self.path, name="thioether")
        self.assertEqual(result["molecules"][0]["matches"][0]["site_atom_indices"], [1])
        hydrogens = selector.select_sdf(self.path, smarts="[H]")
        self.assertEqual(len(hydrogens["molecules"][0]["matches"]), 6)

    def test_invalid_queries_and_files_are_errors(self):
        self.write_molecules(["C"])
        for kwargs in ({}, {"name": "unknown"}, {"smarts": "["}, {"smarts": ""},
                       {"name": "thiol", "smarts": "C"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                selector.select_sdf(self.path, **kwargs)
        with self.assertRaises(ValueError):
            selector.select_sdf(self.path.parent / "missing.sdf", name="thiol")
        self.path.write_text("", encoding="utf-8")
        with self.assertRaises((ValueError, OSError)):
            selector.select_sdf(self.path, name="thiol")

    def test_invalid_record_after_valid_record_is_not_silently_skipped(self):
        self.write_molecules(["CSC"])
        with self.path.open("a", encoding="utf-8") as file:
            file.write("bad\n$$$$\n")
        with self.assertRaisesRegex(ValueError, "record 2"):
            selector.select_sdf(self.path, name="thioether")

    def test_incomplete_record_is_rejected(self):
        self.write_molecules(["CSC"])
        with self.path.open("a", encoding="utf-8") as file:
            file.write("unfinished\n")
        with self.assertRaisesRegex(ValueError, "incomplete SDF record 2"):
            selector.select_sdf(self.path, name="thioether")

    def test_cli_csv_and_list(self):
        self.write_molecules(["CSC"])
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main([str(self.path), "--name", "thioether"]), 0)
        rows = list(csv.DictReader(io.StringIO(output.getvalue())))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["site_atom_indices"], "[1]")
        self.assertEqual(rows[0]["atom_indices"], "[0,1,2]")
        self.assertEqual(rows[0]["atom_index_base"], "0")
        self.assertEqual(rows[0]["sdf_file"], os.path.relpath(self.path))
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main(["--list-selectors"]), 0)
        self.assertIn("thioether:", output.getvalue())

    def test_csv_multiple_selectors_and_relative_paths(self):
        self.write_molecules(["CSC", "c1ccccc1"])
        csv_file = self.path.parent / "input.csv"
        csv_file.write_text(
            "\ufeffpath,ThioEther,pi-pi,hydrogen bond donor\n"
            "molecules.sdf,x,Y,0\n"
            f"{self.path},,x,\n",
            encoding="utf-8",
        )
        result = selector.select_csv(csv_file)
        self.assertEqual([row["csv_line"] for row in result["rows"]], [2, 3])
        first = result["rows"][0]
        self.assertEqual(Path(first["sdf_file"]).resolve(), self.path.resolve())
        self.assertEqual([r["selector"]["name"] for r in first["results"]], ["thioether", "aromatic ring"])
        self.assertEqual(first["results"][0]["molecules"][0]["matches"][0]["site_atom_indices"], [1])
        self.assertEqual(len(first["results"][1]["molecules"][1]["matches"]), 1)
        self.assertEqual(len(result["rows"][1]["results"]), 1)

    def test_single_column_batch_uses_all_selectors_and_reads_sdf_once(self):
        self.write_molecules(["CSC", "CC(=O)[O-]"])
        csv_file = self.path.parent / "input.csv"
        csv_file.write_text("sdf\nmolecules.sdf\n", encoding="utf-8")
        matcher = MolecularSelector()
        with patch.object(matcher, "_read_sdf", wraps=matcher._read_sdf) as read_sdf:
            result = BatchSelector(matcher).select(csv_file)
        read_sdf.assert_called_once()
        results = result["rows"][0]["results"]
        self.assertEqual([item["selector"]["name"] for item in results], list(matcher.registry.definitions))
        thioether = next(item for item in results if item["selector"]["name"] == "thioether")
        self.assertEqual(thioether["molecules"][0]["matches"][0]["site_atom_indices"], [1])
        carboxylate = next(item for item in results if item["selector"]["name"] == "carboxylate")
        self.assertEqual(set(carboxylate["molecules"][1]["matches"][0]["site_atom_indices"]), {2, 3})

    def test_csv_marks_letters_not_nonletters(self):
        self.write_molecules(["CSC"])
        csv_file = self.path.parent / "input.csv"
        csv_file.write_text("path;thiol;thioether;amine\nmolecules.sdf;1;é;+\n", encoding="utf-8")
        result = selector.select_csv(csv_file, delimiter=";")
        self.assertEqual([r["selector"]["name"] for r in result["rows"][0]["results"]], ["thioether"])

    def test_csv_reports_invalid_headers_rows_and_inputs(self):
        self.write_molecules(["CSC"])
        csv_file = self.path.parent / "input.csv"
        cases = [
            ("", "header"),
            ("path\nmolecules.sdf,extra\n", "expected 1 columns"),
            ("path\nmissing.sdf\n", "CSV line 2"),
            ("path\n", "no data"),
            ("path,unknown\nmolecules.sdf,x\n", "unknown selector"),
            ("path,thioether,sulfide\nmolecules.sdf,x,\n", "duplicate CSV"),
            ("path,thiol\nmolecules.sdf\n", "expected 2 columns"),
            ("path,thiol\nmolecules.sdf,x,extra\n", "expected 2 columns"),
            ("path,thiol\n,x\n", "path is empty"),
            ("path,thiol\nmolecules.sdf,1\n", "no selector"),
            ("path,thiol\nmissing.sdf,x\n", "CSV line 2"),
            ("path,thiol\n", "no data"),
        ]
        for content, message in cases:
            with self.subTest(content=content):
                csv_file.write_text(content, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    selector.select_csv(csv_file)
        with self.assertRaisesRegex(ValueError, "delimiter"):
            selector.select_csv(csv_file, delimiter="::")

    def test_csv_cli_and_single_input_are_mutually_exclusive(self):
        self.write_molecules(["CSC"])
        csv_file = self.path.parent / "input.csv"
        csv_file.write_text("file,thioether\nmolecules.sdf,x\n", encoding="utf-8")
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main(["--csv", str(csv_file)]), 0)
        rows = list(csv.DictReader(io.StringIO(output.getvalue())))
        self.assertEqual(rows[0]["selector"], "thioether")
        self.assertEqual(rows[0]["csv_line"], "2")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            selector.main([str(self.path), "--csv", str(csv_file)])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            selector.main(["--csv", str(csv_file), "--name", "thiol"])

    def test_csv_failure_does_not_print_partial_output(self):
        self.write_molecules(["CSC"])
        csv_file = self.path.parent / "input.csv"
        csv_file.write_text("file,thioether\nmolecules.sdf,x\nmissing.sdf,x\n", encoding="utf-8")
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main(["--csv", str(csv_file)]), 1)
        self.assertEqual(output.getvalue(), "")

    def test_csv_output_multiple_matches_and_missing_matches(self):
        self.write_molecules(["CSCCSC", "CC"])
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main([str(self.path), "--name", "thioether"]), 0)
        rows = list(csv.DictReader(io.StringIO(output.getvalue())))
        self.assertEqual(len(rows), 3)
        self.assertEqual([row["record_index"] for row in rows], ["0", "0", "1"])
        self.assertEqual([row["match_index"] for row in rows], ["0", "1", ""])
        self.assertEqual([row["site_atom_indices"] for row in rows], ["[1]", "[4]", ""])
        self.assertEqual([json.loads(row["atom_indices"]) for row in rows[:2]], [[0, 1, 2], [3, 4, 5]])
        self.assertEqual(rows[2]["atom_indices"], "")
        self.assertTrue(all(row["selector"] == "thioether" for row in rows))

    def test_csv_output_custom_smarts_and_aromatic_ring(self):
        self.write_molecules(["c1ccccc1"])
        for args, expected_smarts, expected_selector in (
            (["--smarts", "[c]"], "[c]", ""),
            (["--name", "pi-pi"], "", "aromatic ring"),
        ):
            output = io.StringIO()
            with self.subTest(args=args), contextlib.redirect_stdout(output):
                self.assertEqual(selector.main([str(self.path), *args]), 0)
            rows = list(csv.DictReader(io.StringIO(output.getvalue())))
            self.assertEqual(rows[0]["selector"], expected_selector)
            self.assertEqual(rows[0]["smarts"], expected_smarts)

    def test_csv_output_file_does_not_overwrite_existing_files(self):
        self.write_molecules(["CSC"])
        destination = self.path.parent / "result.csv"
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main([str(self.path), "--name", "thioether", "-o", str(destination)]), 0)
        self.assertEqual(output.getvalue(), "")
        with destination.open(newline="", encoding="utf-8") as file:
            self.assertEqual(list(csv.DictReader(file))[0]["site_atom_indices"], "[1]")
        before = destination.read_bytes()
        self.assertEqual(selector.main([str(self.path), "--name", "thioether", "-o", str(destination)]), 1)
        self.assertEqual(destination.read_bytes(), before)

    def test_batch_csv_output_file_preserves_selectors_and_input(self):
        self.write_molecules(["CSC"])
        input_file = self.path.parent / "input.csv"
        output_file = self.path.parent / "output.csv"
        content = "file,thioether,pi-pi\nmolecules.sdf,x,x\n"
        input_file.write_text(content, encoding="utf-8")
        self.assertEqual(selector.main(["--csv", str(input_file), "-o", str(output_file)]), 0)
        with output_file.open(newline="", encoding="utf-8") as file:
            rows = list(csv.DictReader(file))
        self.assertEqual([row["selector"] for row in rows], ["thioether", "aromatic ring"])
        self.assertEqual([row["csv_line"] for row in rows], ["2", "2"])
        self.assertEqual(rows[1]["match_index"], "")
        self.assertEqual(selector.main(["--csv", str(input_file), "-o", str(input_file)]), 1)
        self.assertEqual(input_file.read_text(encoding="utf-8"), content)

    def test_batch_output_paths_are_relative_to_working_directory(self):
        self.write_molecules(["CSC"])
        input_file = self.path.parent / "input.csv"
        output_file = self.path.parent / "output.csv"
        input_file.write_text(
            f"file,thioether\nmolecules.sdf,x\n{self.path},x\n",
            encoding="utf-8",
        )
        self.assertEqual(selector.main(["--csv", str(input_file), "-o", str(output_file)]), 0)
        with output_file.open(newline="", encoding="utf-8") as file:
            rows = list(csv.DictReader(file))
        self.assertTrue(all(not os.path.isabs(row["sdf_file"]) for row in rows))
        self.assertTrue(all(Path(row["sdf_file"]).resolve() == self.path.resolve() for row in rows))

    def test_single_relative_sdf_output_path(self):
        self.write_molecules(["CSC"])
        relative = os.path.relpath(self.path)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(selector.main([relative, "--name", "thioether"]), 0)
        row = next(csv.DictReader(io.StringIO(output.getvalue())))
        self.assertEqual(row["sdf_file"], relative)
        self.assertFalse(os.path.isabs(row["sdf_file"]))

    def test_failed_batch_does_not_create_output(self):
        self.write_molecules(["CSC"])
        input_file = self.path.parent / "input.csv"
        output_file = self.path.parent / "output.csv"
        input_file.write_text("file,thiol\nmolecules.sdf,x\nmissing.sdf,x\n", encoding="utf-8")
        self.assertEqual(selector.main(["--csv", str(input_file), "-o", str(output_file)]), 1)
        self.assertFalse(output_file.exists())


if __name__ == "__main__":
    unittest.main()
