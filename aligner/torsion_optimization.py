"""Optimize a placed flexible molecule around its alignment axis."""

import json
import math
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


@dataclass(frozen=True)
class OptimizationResult:
    angle: float
    score: float
    matrix: np.ndarray


class TorsionOptimizer:
    def __init__(self, config):
        if not isinstance(config, dict):
            raise ValueError("torsion_optimization must be a mapping")
        allowed = {"enabled", "score", "objective", "steps_per_iteration", "iterations",
                   "include_vdw", "exclude_flexible_site", "exclude_rigid_site",
                   "exclude_flexible_atoms", "exclude_rigid_atoms"}
        unknown = set(config) - allowed
        if unknown:
            raise ValueError(f"unknown torsion_optimization keys: {', '.join(sorted(map(str, unknown)))}")
        self.enabled = config.get("enabled", False)
        self.score = config.get("score", "minimum_distance")
        self.objective = config.get("objective", "maximize")
        self.steps = config.get("steps_per_iteration", 12)
        self.iterations = config.get("iterations", 3)
        self.include_vdw = config.get("include_vdw", False)
        self.exclude_flexible_site = config.get("exclude_flexible_site", True)
        self.exclude_rigid_site = config.get("exclude_rigid_site", False)
        self.exclude_flexible_atoms = config.get("exclude_flexible_atoms", [])
        self.exclude_rigid_atoms = config.get("exclude_rigid_atoms", [])
        for key in ("enabled", "include_vdw", "exclude_flexible_site", "exclude_rigid_site"):
            if type(getattr(self, key)) is not bool:
                raise ValueError(f"torsion_optimization.{key} must be a boolean")
        if self.score not in ("average_distance", "minimum_distance", "mmff94_interaction"):
            raise ValueError("torsion_optimization.score must be average_distance, minimum_distance or mmff94_interaction")
        if self.objective not in ("minimize", "maximize"):
            raise ValueError("torsion_optimization.objective must be minimize or maximize")
        if type(self.steps) is not int or self.steps < 3:
            raise ValueError("torsion_optimization.steps_per_iteration must be an integer >= 3")
        if type(self.iterations) is not int or self.iterations < 1:
            raise ValueError("torsion_optimization.iterations must be a positive integer")
        for key in ("exclude_flexible_atoms", "exclude_rigid_atoms"):
            value = getattr(self, key)
            if not isinstance(value, list) or any(type(i) is not int or i < 0 for i in value):
                raise ValueError(f"torsion_optimization.{key} must be a list of non-negative atom indices")
        if self.score == "mmff94_interaction" and self.include_vdw:
            raise ValueError("torsion_optimization.include_vdw applies only to distance scores")

    @staticmethod
    def excluded(row, molecule, include_site, explicit):
        indices = set(explicit)
        if include_site:
            raw = row.fields.get("site_atom_indices")
            if raw:
                try:
                    sites = json.loads(raw)
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"CSV line {row.line}: invalid site_atom_indices") from exc
                if (not isinstance(sites, list) or not sites
                        or any(type(i) is not int or i < 0 for i in sites)
                        or row.fields.get("atom_index_base", "0") != "0"):
                    raise ValueError(f"CSV line {row.line}: site_atom_indices must be a nonempty 0-based integer list")
                indices.update(sites)
            elif row.center_index is not None:
                indices.add(row.center_index)
            else:
                raise ValueError(f"CSV line {row.line}: site_atom_indices or a center atom is required for exclusion")
        if any(index >= molecule.GetNumAtoms() for index in indices):
            raise ValueError(f"CSV line {row.line}: excluded atom index outside the molecule")
        selected = [index for index in range(molecule.GetNumAtoms()) if index not in indices]
        if not selected:
            raise ValueError(f"CSV line {row.line}: atom exclusions leave no atoms to score")
        return selected

    @staticmethod
    def vdw(mol, indices):
        table = Chem.GetPeriodicTable()
        radii = []
        for index in indices:
            atom = mol.GetAtomWithIdx(index)
            if atom.GetAtomicNum() == 0:
                raise ValueError(f"pseudo-atom {index} has no van der Waals radius")
            radius = table.GetRvdw(atom.GetAtomicNum())
            if not math.isfinite(radius) or radius <= 0:
                raise ValueError(f"atom {index} has no valid van der Waals radius")
            radii.append(radius)
        return np.array(radii)

    @staticmethod
    def mmff_energy(molecule, interfragment):
        Chem.GetSymmSSSR(molecule)
        properties = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94")
        if properties is None:
            raise ValueError("MMFF94 parameters unavailable for a molecule")
        field = AllChem.MMFFGetMoleculeForceField(
            molecule, properties, ignoreInterfragInteractions=not interfragment, nonBondedThresh=1e6)
        if field is None:
            raise ValueError("could not initialize MMFF94 force field")
        energy = field.CalcEnergy()
        if not math.isfinite(energy):
            raise ValueError("MMFF94 returned a non-finite energy")
        return energy

    @staticmethod
    def rotation(axis, angle):
        x, y, z = axis
        skew = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
        return np.eye(3) + math.sin(angle) * skew + (1 - math.cos(angle)) * (skew @ skew)

    def optimize(self, rigid_row, flexible_row, rigid_mol, flexible_mol, pair, matrix, target):
        rigid_count = rigid_mol.GetNumAtoms()
        rigid_positions = np.asarray(pair.GetConformer().GetPositions()[:rigid_count])
        flexible_positions = np.asarray(pair.GetConformer().GetPositions()[rigid_count:])
        axis = matrix @ flexible_row.axis
        axis = axis / np.linalg.norm(axis)
        if self.score == "mmff94_interaction":
            isolated = self.mmff_energy(rigid_mol, False) + self.mmff_energy(flexible_mol, False)
        else:
            rigid_indices = self.excluded(rigid_row, rigid_mol, self.exclude_rigid_site, self.exclude_rigid_atoms)
            flexible_indices = self.excluded(flexible_row, flexible_mol, self.exclude_flexible_site,
                                             self.exclude_flexible_atoms)
            radii = (self.vdw(rigid_mol, rigid_indices)[:, None] +
                     self.vdw(flexible_mol, flexible_indices)[None, :]) if self.include_vdw else 0

        def score(rotation):
            positions = (flexible_positions - target) @ rotation.T + target
            if self.score == "mmff94_interaction":
                candidate = Chem.Mol(pair)
                conformer = candidate.GetConformer()
                for index, point in enumerate(positions, rigid_count):
                    conformer.SetAtomPosition(index, Point3D(*point))
                value = self.mmff_energy(candidate, True) - isolated
            else:
                distances = np.linalg.norm(rigid_positions[rigid_indices, None, :] -
                                           positions[None, flexible_indices, :], axis=2) - radii
                value = np.mean(distances) if self.score == "average_distance" else np.min(distances)
            if not math.isfinite(value):
                raise ValueError("torsion optimization produced a non-finite score")
            return float(value)

        best_angle = 0.0
        best_score = score(np.eye(3))
        for iteration in range(self.iterations):
            span = 2 * math.pi / (self.steps ** iteration)
            center = best_angle
            for step in range(self.steps):
                angle = center - span / 2 + span * step / self.steps
                value = score(self.rotation(axis, angle))
                if (value > best_score if self.objective == "maximize" else value < best_score):
                    best_angle, best_score = angle, value
        rotation = self.rotation(axis, best_angle)
        conf = pair.GetConformer()
        for index, point in enumerate((flexible_positions - target) @ rotation.T + target, rigid_count):
            conf.SetAtomPosition(index, Point3D(*point))
        return OptimizationResult(math.degrees(best_angle) % 360, best_score, rotation)
