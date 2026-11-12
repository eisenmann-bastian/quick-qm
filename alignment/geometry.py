"""Rigid/flexible molecular alignment geometry."""

import math

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D


class PairAligner:
    @staticmethod
    def rotation(source, target):
        cross = np.cross(source, target)
        sine = np.linalg.norm(cross)
        cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
        if sine < 1e-12:
            if cosine > 0:
                return np.eye(3)
            basis = np.eye(3)[np.argmin(np.abs(source))]
            axis = np.cross(source, basis)
            axis /= np.linalg.norm(axis)
            return 2 * np.outer(axis, axis) - np.eye(3)
        axis = cross / sine
        skew = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + skew * sine + skew @ skew * (1 - cosine)

    @staticmethod
    def radius(row, mol):
        if row.center_index is None:
            return 0.0
        atom = mol.GetAtomWithIdx(row.center_index)
        if atom.GetAtomicNum() == 0:
            raise ValueError(f"CSV line {row.line}: a pseudo-atom has no van der Waals radius")
        radius = Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError(f"CSV line {row.line}: no valid van der Waals radius for {atom.GetSymbol()}")
        return float(radius)

    def align(self, rigid, flexible, rigid_mol, flexible_mol, distance, offset=None,
              orientation="stable"):
        if distance is None:
            distance = self.radius(rigid, rigid_mol) + self.radius(flexible, flexible_mol)
        matrix = self.rotation(flexible.axis, -rigid.axis)
        target = rigid.center + distance * rigid.axis
        if offset is not None:
            target = target + offset
        if orientation == "toward-center":
            direction = rigid.center - target
            length = np.linalg.norm(direction)
            if not math.isfinite(length):
                raise ValueError("non-finite center separation")
            if length > 1e-12:
                matrix = self.rotation(-rigid.axis, direction / length) @ matrix
        translation = target - matrix @ flexible.center
        moved = Chem.Mol(flexible_mol)
        conformer = moved.GetConformer()
        for index in range(moved.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            xyz = matrix @ np.array((point.x, point.y, point.z)) + translation
            if not np.isfinite(xyz).all():
                raise ValueError("non-finite coordinates in flexible molecule")
            conformer.SetAtomPosition(index, Point3D(*xyz))
        rigid_conf = rigid_mol.GetConformer()
        for index in range(rigid_mol.GetNumAtoms()):
            point = rigid_conf.GetAtomPosition(index)
            if not all(math.isfinite(value) for value in (point.x, point.y, point.z)):
                raise ValueError("non-finite coordinates in rigid molecule")
        pair = Chem.CombineMols(rigid_mol, moved)
        return pair, matrix, translation, target, distance


