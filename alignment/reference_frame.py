"""Build a local Cartesian frame around a rigid molecular axis."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ReferenceFrame:
    x: np.ndarray
    z: np.ndarray
    method: str
    atom_index: int | None
    point: np.ndarray


class ReferenceSelector:
    def __init__(self, config=None):
        if config is None:
            config = {}
        if not isinstance(config, dict) or set(config) - {"priority", "enabled"}:
            raise ValueError("reference_axis must contain only priority and enabled")
        priority = config.get("priority", ["A", "B", "C"])
        if (not isinstance(priority, list) or not priority
                or any(method not in ("A", "B", "C") for method in priority)
                or len(set(priority)) != len(priority)):
            raise ValueError("reference_axis priority must be a nonempty list of unique A, B, C")
        enabled = config.get("enabled", {})
        if (not isinstance(enabled, dict) or set(enabled) - {"A", "B", "C"}
                or any(type(value) is not bool for value in enabled.values())):
            raise ValueError("reference_axis enabled must map A, B, C to booleans")
        self.priority = [method for method in priority if enabled.get(method, True)]

    def select(self, row, molecule):
        candidates = set()
        if row.center_index is not None:
            atom = molecule.GetAtomWithIdx(row.center_index)
            for neighbor in atom.GetNeighbors():
                for other in neighbor.GetNeighbors():
                    if other.GetIdx() != row.center_index:
                        candidates.add(other.GetIdx())
        conformer = molecule.GetConformer()
        for method in self.priority:
            if method == "C":
                index = None
                point = np.zeros(3)
            elif candidates:
                index = min(candidates) if method == "A" else max(candidates)
                position = conformer.GetAtomPosition(index)
                point = np.array([position.x, position.y, position.z])
            else:
                continue
            vector = point - row.center
            projected = vector - np.dot(vector, row.axis) * row.axis
            length = np.linalg.norm(projected)
            if not np.isfinite(projected).all() or not np.isfinite(length) or length < 1e-10:
                continue
            x = projected / length
            z = np.cross(x, row.axis)
            return ReferenceFrame(x, z, method, index, point)
        raise ValueError(f"CSV line {row.line}: no usable reference axis from enabled methods")
