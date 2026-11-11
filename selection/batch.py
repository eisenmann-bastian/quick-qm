"""Batch selection from CSV input."""

import csv
from pathlib import Path

from .matching import MoleculeSelector


class BatchSelector:
    """Apply selectors to SDF paths, supporting path-only and legacy batch CSVs."""

    def __init__(self, molecule_selector=None):
        self.molecule_selector = molecule_selector if molecule_selector is not None else MoleculeSelector()

    def select(self, csv_file, *, delimiter=","):
        """Apply all selectors for one column, or marked selectors for legacy input."""
        if len(delimiter) != 1:
            raise ValueError("CSV delimiter must be a single character")
        path = Path(csv_file)
        rows = []
        with path.open(encoding="utf-8-sig", newline="") as source:
            reader = csv.reader(source, delimiter=delimiter, strict=True)
            header = next(reader, None)
            if not header or not header[0].strip():
                raise ValueError("CSV must have a path header")
            select_all = len(header) == 1
            selectors = []
            for column in header[1:]:
                canonical, _ = self.molecule_selector.registry.resolve(column)
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
                if not select_all and not selected:
                    raise ValueError(f"CSV line {line}: no selector is marked with a letter")
                sdf_path = Path(value)
                if not sdf_path.is_absolute():
                    sdf_path = path.resolve().parent / sdf_path
                if select_all:
                    try:
                        results = self.molecule_selector.select_all(sdf_path)
                    except (OSError, ValueError, RuntimeError) as exc:
                        raise ValueError(f"CSV line {line}, SDF {sdf_path}: {exc}") from exc
                else:
                    results = []
                    for name in selected:
                        try:
                            results.append(self.molecule_selector.select(sdf_path, name=name))
                        except (OSError, ValueError, RuntimeError) as exc:
                            raise ValueError(f"CSV line {line}, selector {name!r}: {exc}") from exc
                rows.append({"csv_line": line, "sdf_file": str(sdf_path), "results": results})
        if not rows:
            raise ValueError("CSV contains no data rows")
        return {"rows": rows}
