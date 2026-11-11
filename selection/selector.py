#!/usr/bin/env python3
"""Find functional groups or SMARTS matches in an SDF file."""

import argparse
import csv
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.logging import configure_logging, get_logger
from selection.batch import BatchSelector
from selection.definitions import SELECTORS, NamedSelector, SelectorRegistry
from selection.matching import MoleculeSelector
from selection.output import CsvResultWriter

logger = get_logger("selection.selector")

_registry = SelectorRegistry()
_molecule_selector = MoleculeSelector(_registry)
_batch_selector = BatchSelector(_molecule_selector)
_csv_writer = CsvResultWriter()

_normalize_name = _registry.normalize
_named_selector = _registry.resolve
_read_sdf = _molecule_selector._read_sdf
select_sdf = _molecule_selector.select
select_csv = _batch_selector.select
write_results_csv = _csv_writer.write


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf_file", type=Path, nargs="?", help="input SDF file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--name", help="named functional group or interaction site")
    group.add_argument("--smarts", help="SMARTS query")
    group.add_argument("--list-selectors", action="store_true", help="list named selectors")
    group.add_argument("--csv", type=Path, help="batch CSV with SDF paths and named selector columns")
    parser.add_argument("--delimiter", default=",", help="batch CSV delimiter (default: comma)")
    parser.add_argument("-o", "--output", type=Path, help="write CSV to a new file instead of stdout")
    args = parser.parse_args(argv)
    configure_logging()
    if args.list_selectors:
        if args.sdf_file is not None:
            parser.error("--list-selectors does not accept an SDF file")
        if args.output is not None:
            parser.error("--output is not supported with --list-selectors")
        for name, definition in SELECTORS.items():
            print(f"{name}: {definition.smarts if definition.smarts is not None else 'RDKit aromatic rings'}")
        return 0
    if args.csv is not None:
        if args.sdf_file is not None:
            parser.error("--csv does not accept an SDF file")
    else:
        if args.sdf_file is None:
            parser.error("an SDF file is required")
        if args.name is None and args.smarts is None:
            parser.error("specify --name or --smarts")
    try:
        if args.csv is not None:
            result = select_csv(args.csv, delimiter=args.delimiter)
        else:
            result = select_sdf(args.sdf_file, name=args.name, smarts=args.smarts)
        if args.output is None:
            write_results_csv(result, sys.stdout, sdf_file=args.sdf_file)
        else:
            with args.output.open("x", encoding="utf-8", newline="") as destination:
                write_results_csv(result, destination, sdf_file=args.sdf_file)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
