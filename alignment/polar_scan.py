"""Spherical positions and distance-based duplicate suppression."""

import itertools
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class PolarPosition:
    radius: float
    phi: float
    theta: float
    x: float
    y: float
    z: float


class PositionMerger:
    def __init__(self, tolerance):
        self.tolerance = tolerance
        self.cells = {}

    def add(self, position):
        cell = tuple(math.floor(value / self.tolerance) for value in position)
        for delta in itertools.product((-1, 0, 1), repeat=3):
            neighbor = tuple(a + b for a, b in zip(cell, delta))
            for previous in self.cells.get(neighbor, ()):
                if math.dist(position, previous) <= self.tolerance:
                    return False
        self.cells.setdefault(cell, []).append(position)
        return True


def spherical_position(radius, phi, theta, phi_unit, theta_unit):
    p = math.radians(phi) if phi_unit == "degrees" else phi
    t = math.radians(theta) if theta_unit == "degrees" else theta
    return PolarPosition(radius, phi, theta,
                         radius * math.sin(t) * math.cos(p),
                         radius * math.cos(t),
                         radius * math.sin(t) * math.sin(p))
