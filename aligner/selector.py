#!/usr/bin/env python3
"""Find molecular sites and optionally derive alignment axes or an annotated SDF."""

import argparse
import csv
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from molecular_utilities.logging import configure_logging, get_logger
from aligner import LOGGER_NAMESPACE
from aligner.axis_config import DEFAULT_CONFIG, AxisConfig
from molecular_selector.batch import BatchSelector
from molecular_selector.definitions import SELECTORS, NamedSelector, SelectorRegistry
from molecular_selector.matching import MolecularSelector
from molecular_selector.output import CsvResultWriter
from aligner.selection_pipeline import SelectionPipeline

logger = get_logger("selection.selector", namespace=LOGGER_NAMESPACE)

_registry = SelectorRegistry()
_molecule_selector = MolecularSelector(_registry)
_batch_selector = BatchSelector(_molecule_selector)
_csv_writer = CsvResultWriter()
_pipeline = SelectionPipeline(_molecule_selector)

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
    group.add_argument("--csv", type=Path, help="batch CSV with one SDF path column")
    parser.add_argument("--delimiter", default=",", help="batch CSV delimiter (default: comma)")
    parser.add_argument("--format", choices=("selection", "alignment", "sdf"), default="selection",
                        help="output selection CSV (default), alignment CSV, or annotated SDF")
    parser.add_argument("--axes-config", type=Path, default=None,
                        help="YAML axis configuration for alignment or SDF output")
    parser.add_argument("-o", "--output", type=Path, help="write to a new file instead of stdout")
    args = parser.parse_args(argv)
    configure_logging(namespace=LOGGER_NAMESPACE)
    if args.list_selectors:
        if args.sdf_file is not None:
            parser.error("--list-selectors does not accept an SDF file")
        if args.output is not None:
            parser.error("--output is not supported with --list-selectors")
        if args.format != "selection":
            parser.error("--format is not supported with --list-selectors")
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
    if args.format == "sdf" and args.output is None:
        parser.error("--format sdf requires --output")
    if args.axes_config is not None and args.format == "selection":
        parser.error("--axes-config requires --format alignment or sdf")
    try:
        axis_config = AxisConfig.load(args.axes_config or DEFAULT_CONFIG) if args.format != "selection" else None
        if args.csv is not None:
            result = select_csv(args.csv, delimiter=args.delimiter)
        else:
            result = select_sdf(args.sdf_file, name=args.name, smarts=args.smarts)
        result = _pipeline.prepare(result, args.format, sdf_file=args.sdf_file, axis_config=axis_config)
        if args.output is None:
            _pipeline.write(result, sys.stdout, args.format, sdf_file=args.sdf_file)
        else:
            with args.output.open("x", encoding="utf-8", newline="") as destination:
                _pipeline.write(result, destination, args.format, sdf_file=args.sdf_file)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
