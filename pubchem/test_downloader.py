import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import downloader
from core.logging import configure_logging, get_logger


def sdf(cid):
    return f"{cid}\n  PubChem\n\nM  END\n> <PUBCHEM_COMPOUND_CID>\n{cid}\n\n$$$$\n"


class DownloaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.csv_file = self.root / "input.csv"
        self.output = self.root / "output"
        logger_patch = patch.object(downloader, "logger")
        self.logger = logger_patch.start()
        self.addCleanup(logger_patch.stop)
        logging_patch = patch.object(downloader, "configure_logging")
        logging_patch.start()
        self.addCleanup(logging_patch.stop)

    def errors(self):
        return "\n".join(call.args[0] % call.args[1:] for call in self.logger.error.call_args_list)

    def run_csv(self, content, **kwargs):
        self.csv_file.write_text(content, encoding="utf-8")
        return downloader.download_csv(self.csv_file, self.output, **kwargs)

    @patch.object(downloader.pubchempy, "get_cids")
    @patch.object(downloader.pubchempy, "get_sdf")
    def test_direct_identifier_downloads_and_priority(self, get_sdf, get_cids):
        get_sdf.side_effect = [sdf(241), sdf(175), sdf(878)]
        failures = self.run_csv(
            "CID,name,inchi,smiles,formula\n"
            "241,ignored,ignored,ignored,ignored\n"
            ",acetate,ignored,ignored,ignored\n"
            ",,InChI=example,ignored,ignored\n",
            record_type="2d",
        )
        self.assertEqual(failures, 0)
        self.assertEqual(
            [(call.args, call.kwargs) for call in get_sdf.call_args_list],
            [(("241",), {"namespace": "cid", "record_type": "2d"}),
             (("acetate",), {"namespace": "name", "record_type": "2d"}),
             (("InChI=example",), {"namespace": "inchi", "record_type": "2d"})],
        )
        get_cids.assert_not_called()
        self.assertEqual((self.output / "row-2-cid-241.sdf").read_text(), sdf(241))
        self.assertEqual(len(list(self.output.glob("*.sdf"))), 3)

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(702))
    @patch.object(downloader.pubchempy, "get_cids", return_value=[702])
    def test_smiles_identity_default(self, get_cids, get_sdf):
        self.assertEqual(self.run_csv("smiles\nOCC\n"), 0)
        get_cids.assert_called_once_with("OCC", namespace="smiles", searchtype="identity")
        get_sdf.assert_called_once_with("702", namespace="cid", record_type="3d")

    @patch.object(downloader.pubchempy, "get_sdf", side_effect=lambda cid, **kwargs: sdf(cid))
    @patch.object(downloader.pubchempy, "get_cids", return_value=[241, 702, 241])
    def test_result_selection_for_formula_and_smiles(self, get_cids, get_sdf):
        for header in ("formula", "smiles"):
            for selection, expected, failures in (
                ("exact", [], 1), ("first", ["241"], 0), ("all", ["241", "702"], 0),
            ):
                with self.subTest(header=header, selection=selection):
                    output = self.root / f"{header}-{selection}"
                    self.csv_file.write_text(f"{header}\nC\n", encoding="utf-8")
                    get_cids.reset_mock()
                    get_sdf.reset_mock()
                    self.assertEqual(downloader.download_csv(
                        self.csv_file, output, result_selection=selection,
                        filter_mode="similarity", record_type="2d",
                    ), failures)
                    kwargs = {"namespace": header}
                    if header == "smiles":
                        kwargs["searchtype"] = "similarity"
                    get_cids.assert_called_once_with("C", **kwargs)
                    self.assertEqual([call.args[0] for call in get_sdf.call_args_list], expected)
                    self.assertEqual(len(list(output.glob("*.sdf"))), len(expected))

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(241))
    @patch.object(downloader.pubchempy, "get_cids", return_value=[241])
    def test_each_filter_mode_is_sent_to_api(self, get_cids, get_sdf):
        for mode in downloader.FILTER_MODES:
            with self.subTest(mode=mode):
                identifier = downloader.Identifier("smiles", "smiles")
                self.assertEqual(downloader._resolve_cids(identifier, "C", mode, "exact"), ["241"])
                self.assertEqual(get_cids.call_args.kwargs["searchtype"], mode)

    @patch.object(downloader.pubchempy, "get_sdf", side_effect=[None, sdf(702)])
    @patch.object(downloader.pubchempy, "get_cids", return_value=[241, 702])
    def test_all_continues_after_failed_structure(self, get_cids, get_sdf):
        self.assertEqual(self.run_csv("smiles\nC\n", result_selection="all"), 1)
        self.assertIn("CID 241", self.errors())
        self.assertTrue((self.output / "row-2-cid-702.sdf").exists())

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(702))
    @patch.object(downloader.pubchempy, "get_cids", return_value=[])
    def test_empty_search_does_not_stop_later_rows(self, get_cids, get_sdf):
        self.assertEqual(self.run_csv("smiles,cid\nC,\n,702\n"), 1)
        self.assertIn("returned no CIDs", self.errors())
        self.assertTrue((self.output / "row-3-cid-702.sdf").exists())

    @patch.object(downloader.pubchempy, "get_sdf")
    def test_direct_multiple_and_invalid_records_are_rejected(self, get_sdf):
        for content in (None, "", "error", sdf(241) + sdf(702), "a\n$$$$\n",
                        sdf(241).replace("241\n\n$$$$", "702\n\n$$$$")):
            with self.subTest(content=content):
                get_sdf.return_value = content
                self.assertEqual(self.run_csv("cid\n241\n"), 1)
                self.assertEqual(list(self.output.glob("*.sdf")), [])

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(241))
    def test_existing_file_is_not_overwritten(self, get_sdf):
        self.output.mkdir()
        existing = self.output / "row-2-cid-241.sdf"
        existing.write_text("existing", encoding="utf-8")
        self.assertEqual(self.run_csv("cid\n241\n"), 1)
        self.assertEqual(existing.read_text(), "existing")

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(241))
    def test_header_delimiter_options_and_row_errors(self, get_sdf):
        self.assertEqual(self.run_csv(" CID ;name\n241;ignored\n", delimiter=";"), 0)
        for content in ("other\nx\n", "cid,CID\n1,2\n", ""):
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.run_csv(content)
        for kwargs in ({"delimiter": ";;"}, {"record_type": "invalid"},
                       {"result_selection": "invalid"}, {"filter_mode": "invalid"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.run_csv("cid\n1\n", **kwargs)
        self.assertEqual(self.run_csv("cid\ninvalid\n\n,extra\n"), 2)

    @patch.object(downloader.pubchempy, "get_sdf", side_effect=OSError("network error"))
    def test_network_failure_reported(self, get_sdf):
        self.assertEqual(self.run_csv("name\nbenzene\n"), 1)
        self.assertIn("network error", self.errors())

    @patch.object(downloader.pubchempy, "get_cids", side_effect=downloader.pubchempy.PubChemPyError("search failed"))
    def test_api_error_reported(self, get_cids):
        self.assertEqual(self.run_csv("formula\nCH4\n"), 1)
        self.assertIn("search failed", self.errors())

    @patch.object(downloader.pubchempy, "get_cids", return_value=[241, "bad"])
    @patch.object(downloader.pubchempy, "get_sdf")
    def test_invalid_search_result_does_not_download(self, get_sdf, get_cids):
        self.assertEqual(self.run_csv("formula\nC\n", result_selection="all"), 1)
        self.assertIn("invalid PubChem CID", self.errors())
        get_sdf.assert_not_called()

    @patch.object(downloader.pubchempy, "get_sdf", return_value=sdf(241))
    @patch.object(downloader.pubchempy, "get_cids", return_value=[241, 702])
    def test_cli_options(self, get_cids, get_sdf):
        self.csv_file.write_text("smiles\nC\n", encoding="utf-8")
        status = downloader.main([
            str(self.csv_file), "-o", str(self.output), "--filter-mode", "superstructure",
            "--result-selection", "first",
        ])
        self.assertEqual(status, 0)
        get_cids.assert_called_once_with("C", namespace="smiles", searchtype="superstructure")
        get_sdf.assert_called_once_with("241", namespace="cid", record_type="3d")

    @patch.object(downloader.pubchempy, "get_cids", return_value=[241, 702])
    def test_cli_defaults_reject_multiple_matches(self, get_cids):
        self.csv_file.write_text("smiles\nC\n", encoding="utf-8")
        status = downloader.main([str(self.csv_file), "-o", str(self.output)])
        self.assertEqual(status, 1)
        self.assertIn("matched 2 CIDs", self.errors())
        get_cids.assert_called_once_with("C", namespace="smiles", searchtype="identity")


class SharedLoggingTests(unittest.TestCase):
    def test_configuration_is_shared_and_idempotent(self):
        import io

        parent = logging.getLogger("quick_qm")
        old_handlers, old_level, old_propagate = parent.handlers[:], parent.level, parent.propagate
        self.addCleanup(self.restore_logger, parent, old_handlers, old_level, old_propagate)
        parent.handlers = []
        stream = io.StringIO()
        configure_logging(stream=stream)
        configure_logging(stream=stream)
        get_logger("another_tool").info("hello")
        self.assertEqual(stream.getvalue().count("hello"), 1)
        self.assertIn("INFO quick_qm.another_tool:", stream.getvalue())

    @staticmethod
    def restore_logger(logger, handlers, level, propagate):
        for handler in logger.handlers:
            if handler not in handlers:
                handler.close()
        logger.handlers = handlers
        logger.setLevel(level)
        logger.propagate = propagate


if __name__ == "__main__":
    unittest.main()
