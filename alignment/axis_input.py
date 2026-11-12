"""Selector axis CSV input and filtering."""

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from molecular_selector.matching import MolecularSelector

FILTER_FIELDS = ("molecule_name", "selector", "axis_type")
REQUIRED_FIELDS = (*FILTER_FIELDS, "sdf_file", "record_index", "center_atom_index", "center", "axis")


@dataclass(frozen=True)
class AxisRow:
    line: int
    fields: dict
    path: Path
    record_index: int
    center_index: int | None
    center: np.ndarray
    axis: np.ndarray


class AxisInput:
    """Load selected axes and their referenced source molecules."""

    def __init__(self):
        self.molecules = {}

    @staticmethod
    def vector(value, field):
        try:
            data = json.loads(value)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{field} must be a JSON array of three numbers") from exc
        if (not isinstance(data, list) or len(data) != 3
                or any(type(x) not in (int, float) or not math.isfinite(x) for x in data)):
            raise ValueError(f"{field} must be a JSON array of three finite numbers")
        return np.asarray(data, dtype=float)

    def read(self, path, rigid, flexible):
        left, right = [], []
        with Path(path).open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, strict=True)
            if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise ValueError("CSV must have unique headers")
            missing = set(REQUIRED_FIELDS) - set(reader.fieldnames)
            if missing:
                raise ValueError(f"CSV missing required columns: {', '.join(sorted(missing))}")
            for row in reader:
                line = reader.line_num
                if None in row or any(row[key] is None for key in reader.fieldnames):
                    raise ValueError(f"CSV line {line}: incorrect number of columns")
                is_left = rigid.matches(row)
                is_right = flexible.matches(row)
                if not (is_left or is_right):
                    continue
                try:
                    sdf = Path(row["sdf_file"])
                    if not row["sdf_file"] or not sdf.is_file():
                        raise ValueError(f"SDF file does not exist: {sdf}")
                    record_index = int(row["record_index"])
                    if record_index < 0:
                        raise ValueError("record_index must be non-negative")
                    center_index = int(row["center_atom_index"]) if row["center_atom_index"] else None
                    center = self.vector(row["center"], "center")
                    axis = self.vector(row["axis"], "axis")
                    length = np.linalg.norm(axis)
                    if not math.isfinite(length) or length < 1e-12:
                        raise ValueError("axis must have nonzero finite length")
                    sdf = sdf.resolve()
                    if sdf not in self.molecules:
                        self.molecules[sdf] = list(MolecularSelector.read_sdf(sdf))
                    molecules = self.molecules[sdf]
                    if record_index >= len(molecules):
                        raise ValueError("record_index is outside the SDF")
                    mol = molecules[record_index]
                    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                        raise ValueError("SDF record must have 3D coordinates")
                    if center_index is not None:
                        if center_index < 0 or center_index >= mol.GetNumAtoms():
                            raise ValueError("center_atom_index is outside the SDF")
                        point = mol.GetConformer().GetAtomPosition(center_index)
                        if not np.allclose(center, (point.x, point.y, point.z), rtol=0, atol=1e-3):
                            raise ValueError("center does not match the SDF atom coordinates")
                    elif row["axis_type"] != "aromatic_normal":
                        raise ValueError("center_atom_index is required except for aromatic ring axes")
                    entry = AxisRow(line, dict(row), sdf, record_index, center_index, center, axis / length)
                except (ValueError, OSError, RuntimeError) as exc:
                    raise ValueError(f"CSV line {line}: {exc}") from exc
                if is_left:
                    left.append(entry)
                if is_right:
                    right.append(entry)
        if not left or not right:
            raise ValueError(f"filters matched {len(left)} rigid and {len(right)} flexible axes")
        return left, right

    def molecule(self, row):
        return self.molecules[row.path][row.record_index]


class AxisFilter:
    def __init__(self, config=None, include_override=None, exclude_override=None):
        if config is None:
            config = {}
        if (isinstance(config, list)
                or (isinstance(config, dict) and config and set(config) <= set(FILTER_FIELDS))):
            config = {"include": config}
        if not isinstance(config, dict) or set(config) - {"include", "exclude"}:
            raise ValueError("axis filters must contain only include and exclude")
        include = include_override if include_override is not None else config.get("include", [])
        exclude = exclude_override if exclude_override is not None else config.get("exclude", [])
        self.include = self.compile_filters(include)
        self.exclude = self.compile_filters(exclude)

    @staticmethod
    def compile_filters(expressions):
        if isinstance(expressions, dict):
            expressions = [expressions]
        if not isinstance(expressions, list):
            raise ValueError("include and exclude must be lists of filters")
        filters = []
        for item in expressions:
            if not isinstance(item, dict) or not item or set(item) - set(FILTER_FIELDS):
                raise ValueError("each filter must contain at least one of molecule_name, selector "
                                 "and axis_type, and no other fields")
            try:
                filters.append({key: re.compile(value) for key, value in item.items()})
            except (re.error, TypeError) as exc:
                raise ValueError(f"invalid filter regex: {exc}") from exc
        return filters

    def matches(self, row):
        def matches_any(filters):
            return any(all(pattern.search(row[key]) is not None for key, pattern in group.items())
                       for group in filters)

        return (not self.include or matches_any(self.include)) and not matches_any(self.exclude)

