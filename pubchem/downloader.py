#!/usr/bin/env python3
"""Download PubChem SDFs from CSV identifiers."""

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import pubchempy

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.logging import configure_logging, get_logger

logger = get_logger("pubchem.downloader")
FILTER_MODES = ("identity", "substructure", "superstructure", "similarity")
RESULT_SELECTIONS = ("exact", "first", "all")


@dataclass(frozen=True)
class Identifier:
    header: str
    namespace: str


class CsvHeaders:
    """Identifier headers and their PubChemPy namespaces, in priority order."""

    identifiers = (
        Identifier("cid", "cid"),
        Identifier("name", "name"),
        Identifier("inchi", "inchi"),
        Identifier("smiles", "smiles"),
        Identifier("formula", "formula"),
    )

    def __init__(self, fieldnames):
        if not fieldnames:
            raise ValueError("CSV must have a header")
        self.fieldnames = [header.strip().lower() for header in fieldnames]
        if len(set(self.fieldnames)) != len(self.fieldnames):
            raise ValueError("CSV contains duplicate headers")
        supported = [identifier.header for identifier in self.identifiers]
        if not set(self.fieldnames).intersection(supported):
            raise ValueError(f"CSV must include an identifier column: {', '.join(supported)}")

    def select(self, row):
        for identifier in self.identifiers:
            value = row.get(identifier.header)
            if value and value.strip():
                return identifier, value.strip()
        raise ValueError("row has no populated identifier")


def _normalize_cid(result):
    cid = str(result).strip()
    if not cid.isascii() or not cid.isdigit() or int(cid) < 1:
        raise ValueError(f"invalid PubChem CID: {result!r}")
    return str(int(cid))


def _resolve_cids(identifier, value, filter_mode, result_selection):
    kwargs = {"searchtype": filter_mode} if identifier.namespace == "smiles" else {}
    results = pubchempy.get_cids(value, namespace=identifier.namespace, **kwargs)
    if not isinstance(results, list):
        raise ValueError(f"invalid {identifier.header} search result")
    if not results:
        raise ValueError(f"{identifier.header} returned no CIDs")
    cids = list(dict.fromkeys(_normalize_cid(cid) for cid in results))
    if result_selection == "exact" and len(cids) != 1:
        raise ValueError(f"{identifier.header} matched {len(cids)} CIDs; expected exactly one")
    return cids[:1] if result_selection == "first" else cids


def _read_sdf(identifier, namespace, record_type):
    sdf = pubchempy.get_sdf(identifier, namespace=namespace, record_type=record_type)
    if not isinstance(sdf, str) or not sdf.strip():
        raise ValueError(f"no SDF returned for {namespace} {identifier}")
    records = re.split(r"(?m)^\$\$\$\$\s*$", sdf)
    if records[-1].strip() or len(records) != 2 or not records[0].strip():
        raise ValueError(f"expected exactly one complete SDF record for {namespace} {identifier}")
    matches = re.findall(r"(?m)^>\s*<PUBCHEM_COMPOUND_CID>[^\r\n]*\r?\n(\d+)\s*$", records[0])
    if len(matches) != 1:
        raise ValueError(f"SDF for {namespace} {identifier} is missing a unique PubChem CID")
    cid = _normalize_cid(matches[0])
    if namespace == "cid" and cid != _normalize_cid(identifier):
        raise ValueError(f"SDF CID {cid} does not match requested CID {identifier}")
    return cid, sdf


def _save_sdf(output_dir, line, cid, sdf):
    destination = output_dir / f"row-{line}-cid-{cid}.sdf"
    with destination.open("x", encoding="utf-8", newline="") as output:
        output.write(sdf)
    logger.info("Row %s: saved %s", line, destination)


def download_csv(
    csv_file, output_dir, *, delimiter=",", record_type=None,
    result_selection="exact", filter_mode="identity",
):
    """Download rows and return the number of failed rows."""
    if len(delimiter) != 1:
        raise ValueError("delimiter must be a single character")
    if record_type not in (None, "2d", "3d"):
        raise ValueError("record_type must be '2d', '3d', or None")
    if result_selection not in RESULT_SELECTIONS:
        raise ValueError(f"result_selection must be one of {', '.join(RESULT_SELECTIONS)}")
    if filter_mode not in FILTER_MODES:
        raise ValueError(f"filter_mode must be one of {', '.join(FILTER_MODES)}")

    output_dir = Path(output_dir)
    failures = 0
    with Path(csv_file).open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter=delimiter, strict=True)
        headers = CsvHeaders(reader.fieldnames)
        reader.fieldnames = headers.fieldnames
        output_dir.mkdir(parents=True, exist_ok=True)

        for row in reader:
            line = reader.line_num
            try:
                if None in row:
                    raise ValueError("row has more values than the header")
                identifier, value = headers.select(row)
                if identifier.namespace in ("cid", "name", "inchi"):
                    if identifier.namespace == "cid":
                        value = _normalize_cid(value)
                    cid, sdf = _read_sdf(value, identifier.namespace, record_type or "3d")
                    _save_sdf(output_dir, line, cid, sdf)
                    continue
                cids = _resolve_cids(identifier, value, filter_mode, result_selection)
            except (ValueError, OSError, pubchempy.PubChemPyError) as exc:
                failures += 1
                logger.error("Row %s: %s", line, exc)
                continue
            row_failed = False
            for cid in cids:
                try:
                    _, sdf = _read_sdf(cid, "cid", record_type or "3d")
                    _save_sdf(output_dir, line, cid, sdf)
                except (ValueError, OSError, pubchempy.PubChemPyError) as exc:
                    row_failed = True
                    logger.error("Row %s, CID %s: %s", line, cid, exc)
            if row_failed:
                failures += 1
    return failures


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", type=Path, help="CSV file with identifier headers")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("sdf"))
    parser.add_argument("--delimiter", default=",", help="CSV delimiter (default: comma)")
    parser.add_argument("--record-type", choices=("2d", "3d"), default="3d")
    parser.add_argument(
        "--result-selection", choices=RESULT_SELECTIONS, default="exact",
        help="for formula and SMILES: require one CID (default), or select first or all",
    )
    parser.add_argument(
        "--filter-mode", choices=FILTER_MODES, default="identity",
        help="PubChem SMILES search type (default: identity)",
    )
    args = parser.parse_args(argv)
    configure_logging()
    try:
        failures = download_csv(
            args.csv_file, args.output_dir, delimiter=args.delimiter, record_type=args.record_type,
            result_selection=args.result_selection, filter_mode=args.filter_mode,
        )
    except (ValueError, OSError, csv.Error, pubchempy.PubChemPyError) as exc:
        logger.error("%s", exc)
        return 1
    if failures:
        logger.error("%s row(s) failed.", failures)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
