#!/usr/bin/env python3
"""Find functional groups or SMARTS matches in an SDF file."""

import argparse
import csv
import io
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Lipinski

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.logging import configure_logging, get_logger

logger = get_logger("selection.selector")


@dataclass(frozen=True)
class NamedSelector:
    smarts: str | None
    aliases: tuple[str, ...] = ()


SELECTORS = {
    "thioether": NamedSelector("[#6]-[S;X2;H0:1]-[#6]", ("sulfide",)),
    "thiol": NamedSelector("[S;X2;H1;+0:1]-[#6]"),
    "carbonyl": NamedSelector("[C]=[O:1]"),
    "alcohol": NamedSelector("[O;X2;H1;+0:1]-[C;X4]"),
    "phenol": NamedSelector("[O;X2;H1;+0:1]-[c]"),
    "ether": NamedSelector("[#6]-[O;X2;H0:1]-[#6]"),
    "amine": NamedSelector("[N;X3;+0;!$(N-[C,S,P]=[O,S,N]):1]"),
    "amide": NamedSelector("[N:1]-[C](=[O:2])"),
    "carboxylic acid": NamedSelector("[C](=[O:1])-[O;H1;+0:2]"),
    "ester": NamedSelector("[C](=[O:1])-[O;H0:2]-[#6]"),
    "aldehyde": NamedSelector("[C;H1](=[O:1])-[#6]"),
    "ketone": NamedSelector("[#6]-[C](=[O:1])-[#6]"),
    "hydrogen bond donor": NamedSelector(
        Chem.MolToSmarts(Lipinski.HDonorSmarts),
        ("h bond donor", "hbond donor"),
    ),
    "hydrogen bond acceptor": NamedSelector(
        Chem.MolToSmarts(Lipinski.HAcceptorSmarts),
        ("h bond acceptor", "hbond acceptor"),
    ),
    "halogen bond donor": NamedSelector(
        "[#6]-[Cl,Br,I;+0:1]", ("halogen donor",)
    ),
    "halogen bond acceptor": NamedSelector(
        Chem.MolToSmarts(Lipinski.HAcceptorSmarts), ("halogen acceptor",)
    ),
    "positive charge": NamedSelector(
        "[+1,+2,+3,+4:1]", ("cation", "ionic positive", "cation pi cation", "ionic interaction cation")
    ),
    "negative charge": NamedSelector(
        "[-1,-2,-3,-4:1]", ("anion", "ionic negative", "ionic interaction anion")
    ),
    "aromatic ring": NamedSelector(
        None, ("aromatic group", "pi pi", "pi pi ring", "pi pi interaction", "cation pi ring")
    ),
    "chalcogen bond donor": NamedSelector(
        "[#6]-[S,Se,Te;X2;H0;+0:1]-[#6]", ("chalcogen donor",)
    ),
    "chalcogen bond acceptor": NamedSelector(
        Chem.MolToSmarts(Lipinski.HAcceptorSmarts), ("chalcogen acceptor",)
    ),
    "metal center": NamedSelector(
        "[Li,Na,K,Mg,Ca,Mn,Fe,Co,Ni,Cu,Zn,Cd,Hg,Pt,Pd;+1,+2,+3,+4:1]",
        ("metal coordination center",),
    ),
    "metal coordination ligand": NamedSelector(
        Chem.MolToSmarts(Lipinski.HAcceptorSmarts), ("metal ligand",)
    ),
}


def _normalize_name(name):
    return " ".join(re.sub(r"[-_]+", " ", name).casefold().split())


def _named_selector(name):
    normalized = _normalize_name(name)
    for canonical, definition in SELECTORS.items():
        if normalized in (_normalize_name(canonical), *map(_normalize_name, definition.aliases)):
            return canonical, definition
    raise ValueError(f"unknown selector {name!r}; use --list-selectors to see available names")


def _read_sdf(path):
    with path.open("rb") as source:
        lines = []
        record_index = 0
        for line in source:
            lines.append(line)
            if line.rstrip(b"\r\n") != b"$$$$":
                continue
            supplier = Chem.ForwardSDMolSupplier(
                io.BytesIO(b"".join(lines)), removeHs=False, strictParsing=True
            )
            mol = next(supplier, None)
            if mol is None:
                raise ValueError(f"could not parse SDF record {record_index + 1}")
            yield mol
            record_index += 1
            lines = []
        if b"".join(lines).strip():
            raise ValueError(f"incomplete SDF record {record_index + 1}: missing $$$$ delimiter")


def select_sdf(sdf_file, *, name=None, smarts=None):
    """Return one entry per SDF record with zero-based atom-index matches."""
    if (name is None) == (smarts is None):
        raise ValueError("specify exactly one of name or smarts")
    if name is not None:
        selector_name, definition = _named_selector(name)
        query_smarts = definition.smarts
    else:
        selector_name = None
        query_smarts = smarts
    aromatic_ring = selector_name == "aromatic ring"
    if not aromatic_ring:
        if not isinstance(query_smarts, str) or not query_smarts.strip():
            raise ValueError("SMARTS must not be empty")
        query = Chem.MolFromSmarts(query_smarts)
        if query is None or query.GetNumAtoms() == 0:
            raise ValueError(f"invalid SMARTS: {query_smarts!r}")
        site_positions = [
            atom.GetIdx() for atom in query.GetAtoms() if atom.GetAtomMapNum() > 0
        ] if name is not None else list(range(query.GetNumAtoms()))
        if not site_positions:
            site_positions = list(range(query.GetNumAtoms()))

    path = Path(sdf_file)
    if not path.is_file():
        raise ValueError(f"SDF file does not exist: {path}")
    molecules = []
    for index, mol in enumerate(_read_sdf(path)):
        matches = []
        seen_sites = set()
        if aromatic_ring:
            atom_matches = (
                ring for ring in mol.GetRingInfo().AtomRings()
                if all(
                    mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)]).GetIsAromatic()
                    for i in range(len(ring))
                )
            )
        else:
            atom_matches = mol.GetSubstructMatches(query, uniquify=True, maxMatches=0)
        for atoms in atom_matches:
            sites = tuple(atoms) if aromatic_ring else tuple(atoms[position] for position in site_positions)
            key = tuple(sorted(sites))
            if name is not None and key in seen_sites:
                continue
            seen_sites.add(key)
            matches.append({"atom_indices": list(atoms), "site_atom_indices": list(sites)})
        molecules.append({
            "record_index": index,
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
            "matches": matches,
        })
    if not molecules:
        raise ValueError("SDF file contains no molecule records")
    return {
        "selector": {"name": selector_name, "smarts": query_smarts},
        "atom_index_base": 0,
        "molecules": molecules,
    }


def select_csv(csv_file, *, delimiter=","):
    """Run named selectors marked with letters in a CSV and return row results."""
    if len(delimiter) != 1:
        raise ValueError("CSV delimiter must be a single character")
    path = Path(csv_file)
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter=delimiter, strict=True)
        header = next(reader, None)
        if not header or len(header) < 2 or not header[0].strip():
            raise ValueError("CSV must have a path header followed by selector headers")
        selectors = []
        for column in header[1:]:
            canonical, _ = _named_selector(column)
            if canonical in selectors:
                raise ValueError(f"duplicate CSV selector column: {column!r}")
            selectors.append(canonical)
        for row in reader:
            line = reader.line_num
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"CSV line {line}: expected {len(header)} columns, got {len(row)}")
            value = row[0].strip()
            if not value:
                raise ValueError(f"CSV line {line}: SDF path is empty")
            selected = [
                name for name, cell in zip(selectors, row[1:])
                if any(character.isalpha() for character in cell)
            ]
            if not selected:
                raise ValueError(f"CSV line {line}: no selector is marked with a letter")
            sdf_path = Path(value)
            if not sdf_path.is_absolute():
                sdf_path = path.resolve().parent / sdf_path
            results = []
            for name in selected:
                try:
                    results.append(select_sdf(sdf_path, name=name))
                except (OSError, ValueError, RuntimeError) as exc:
                    raise ValueError(f"CSV line {line}, selector {name!r}: {exc}") from exc
            rows.append({"csv_line": line, "sdf_file": str(sdf_path), "results": results})
    if not rows:
        raise ValueError("CSV contains no data rows")
    return {"rows": rows}


def write_results_csv(result, destination, *, sdf_file=None):
    """Write one CSV row per match, preserving records without matches."""
    writer = csv.writer(destination, lineterminator="\n")
    writer.writerow([
        "csv_line", "sdf_file", "record_index", "molecule_name", "selector",
        "smarts", "atom_index_base", "match_index", "atom_indices", "site_atom_indices",
    ])
    if "rows" in result:
        groups = (
            (row["csv_line"], row["sdf_file"], item)
            for row in result["rows"] for item in row["results"]
        )
    else:
        groups = [(None, str(sdf_file), result)]
    for line, path, item in groups:
        for molecule in item["molecules"]:
            matches = molecule["matches"]
            common = [
                line, path, molecule["record_index"], molecule["name"],
                item["selector"]["name"], item["selector"]["smarts"], item["atom_index_base"],
            ]
            if not matches:
                writer.writerow(common + ["", "", ""])
            for index, match in enumerate(matches):
                writer.writerow(common + [
                    index,
                    json.dumps(match["atom_indices"], separators=(",", ":")),
                    json.dumps(match["site_atom_indices"], separators=(",", ":")),
                ])


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
