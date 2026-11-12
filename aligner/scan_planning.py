"""Configuration parsing and distance scan planning."""

import argparse
import math
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path

import yaml

from aligner.axis_input import FILTER_FIELDS
from aligner.axis_config import UniqueKeyLoader


@dataclass(frozen=True)
class DistancePlan:
    mode: str
    start: Fraction | None = None
    step: Fraction | None = None
    first: int = 0
    count: int = 1
    distance: float | None = None

    @staticmethod
    def number(value, name):
        if type(value) not in (int, float, str):
            raise ValueError(f"{name} must be a finite number")
        try:
            number = Decimal(str(value))
        except InvalidOperation as exc:
            raise ValueError(f"{name} must be a finite number") from exc
        if not number.is_finite() or not math.isfinite(float(number)):
            raise ValueError(f"{name} must be a finite number")
        return Fraction(number)

    @classmethod
    def from_options(cls, settings, args, allow_negative=False, default_distance=None):
        def option(name, default=None):
            value = getattr(args, name)
            return value if value is not None else settings.get(name, default)

        mode = option("alignment_mode", "default")
        distance = option("distance")
        scan_names = ("scan_start", "scan_end", "include_start", "include_end", "scan_step", "scan_count")
        scan = {name: option(name) for name in scan_names}
        if mode not in ("default", "distance-scan"):
            raise ValueError("alignment_mode must be 'default' or 'distance-scan'")
        if mode == "default":
            if any(value is not None for value in scan.values()):
                raise ValueError("scan options require alignment_mode distance-scan")
            if distance is None:
                return cls(mode, distance=default_distance)
            value = cls.number(distance, "distance")
            if value < 0 and not allow_negative:
                raise ValueError("distance must be non-negative")
            return cls(mode, distance=float(value))
        if distance is not None:
            raise ValueError("distance cannot be used with distance-scan")
        if scan["scan_start"] is None or scan["scan_end"] is None:
            raise ValueError("distance-scan requires scan_start and scan_end")
        start = cls.number(scan["scan_start"], "scan_start")
        end = cls.number(scan["scan_end"], "scan_end")
        if (start < 0 and not allow_negative) or end <= start:
            raise ValueError("scan range must have scan_start < scan_end" +
                             ("" if allow_negative else " and 0 <= scan_start"))
        for name in ("include_start", "include_end"):
            if scan[name] is None:
                scan[name] = True
            elif type(scan[name]) is not bool:
                raise ValueError(f"{name} must be a boolean")
        # A CLI sampling choice overrides the alternative from the config.
        if args.scan_step is not None:
            scan["scan_count"] = None
        elif args.scan_count is not None:
            scan["scan_step"] = None
        if (scan["scan_step"] is None) == (scan["scan_count"] is None):
            raise ValueError("distance-scan requires exactly one of scan_step and scan_count")
        include_start, include_end = scan["include_start"], scan["include_end"]
        first = 0 if include_start else 1
        if scan["scan_step"] is not None:
            step = cls.number(scan["scan_step"], "scan_step")
            if step <= 0:
                raise ValueError("scan_step must be positive")
            intervals = (end - start) / step
            if intervals.denominator != 1:
                raise ValueError("scan range must be evenly divisible by scan_step")
            count = intervals.numerator + 1 - int(not include_start) - int(not include_end)
            if count <= 0:
                raise ValueError("scan produces no poses")
        else:
            interior_count = scan["scan_count"]
            if type(interior_count) is not int or interior_count < 1:
                raise ValueError("scan_count must be a positive integer")
            intervals = interior_count + 1
            count = interior_count + int(include_start) + int(include_end)
            step = (end - start) / intervals
        return cls(mode, start, step, first, count)

    def distances(self):
        if self.mode == "default":
            yield self.distance
            return
        for index in range(self.first, self.first + self.count):
            yield float(self.start + index * self.step)


def load_config(path):
    if path is None:
        return {}
    try:
        with Path(path).open(encoding="utf-8") as source:
            data = yaml.load(source, Loader=UniqueKeyLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"invalid configuration {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("configuration must be a YAML mapping")
    unknown = set(data) - {"rigid", "flexible", "output", "output_format", "overwrite_output",
                           "grid_based_scan", "polar_coordinate_based_scan", "scans", "scan_merge_tolerance",
                           "reference_axis", "flexible_orientation", "torsion_optimization"}
    if unknown:
        raise ValueError(f"unknown configuration keys: {', '.join(sorted(map(str, unknown)))}")
    return data


def triplets(values):
    return [dict(zip(FILTER_FIELDS, group)) for group in values]


def axis_plan(axis, scans, args, *, polar=False):
    scan_keys = ("scan_start", "scan_end", "include_start", "include_end", "scan_step", "scan_count")
    keys = ("distance", *scan_keys)
    values = scans.get(axis, {})
    section = "polar_coordinate_based_scan" if polar else "grid_based_scan"
    if not isinstance(values, dict) or set(values) - set(keys):
        raise ValueError(f"{section}.{axis} must be a mapping of distance and scan options")
    prefix = f"polar_{axis}_" if polar else ("" if axis == "y" else f"{axis}_")
    options = argparse.Namespace(**{key: getattr(args, prefix + key) for key in ("alignment_mode", *keys)})
    requested = bool(values) or any(getattr(options, key) is not None for key in ("alignment_mode", *keys))
    values = dict(values)
    cli_scan = any(getattr(options, key) is not None for key in scan_keys)
    if options.distance is not None and not cli_scan and options.alignment_mode != "distance-scan":
        values = {key: value for key, value in values.items() if key not in scan_keys}
    elif options.alignment_mode == "distance-scan" or cli_scan:
        values.pop("distance", None)
    if options.alignment_mode == "default":
        values = {key: value for key, value in values.items() if key not in scan_keys}
    values["alignment_mode"] = "distance-scan" if any(key in values for key in scan_keys) or cli_scan else "default"
    try:
        plan = DistancePlan.from_options(values, options, allow_negative=not polar and axis != "y",
                                         default_distance=None if axis in ("y", "radius") else 0.0)
    except ValueError as exc:
        raise ValueError(f"{section}.{axis}: {exc}") from exc
    return plan, requested


def polar_options(scans, args):
    if not isinstance(scans, dict) or set(scans) - {"radius", "phi", "theta", "phi_unit", "theta_unit", "merge_tolerance"}:
        raise ValueError("polar_coordinate_based_scan must contain only radius, phi, theta, angle units and merge_tolerance")
    units = []
    plans = []
    for angle in ("phi", "theta"):
        unit = getattr(args, f"polar_{angle}_unit") or scans.get(f"{angle}_unit", "degrees")
        if unit not in ("degrees", "radians"):
            raise ValueError(f"polar {angle}_unit must be degrees or radians")
        plan, _ = axis_plan(angle, scans, args, polar=True)
        maximum = (360 if angle == "phi" else 180) if unit == "degrees" else (2 * math.pi if angle == "phi" else math.pi)
        bounds = [plan.distance] if plan.mode == "default" else [float(plan.start), float(plan.start + plan.step * (plan.count - 1 + plan.first))]
        # Validate declared boundaries too, including excluded endpoints.
        if plan.mode == "distance-scan":
            raw_end = getattr(args, f"polar_{angle}_scan_end")
            bounds.append(float(DistancePlan.number(raw_end if raw_end is not None else scans[angle]["scan_end"], "scan_end")))
        if any(value is None or value < 0 or value > maximum + 1e-12 for value in bounds):
            raise ValueError(f"polar {angle} must be between 0 and {maximum} {unit}")
        units.append(unit)
        plans.append(plan)
    tolerance = args.polar_merge_tolerance if args.polar_merge_tolerance is not None else scans.get("merge_tolerance", 0.05)
    tolerance = float(DistancePlan.number(tolerance, "merge_tolerance"))
    if tolerance <= 0:
        raise ValueError("polar merge_tolerance must be positive")
    radius, _ = axis_plan("radius", scans, args, polar=True)
    return radius, plans[0], plans[1], units[0], units[1], tolerance


@dataclass(frozen=True)
class ScanPlan:
    system: str
    plans: tuple
    needs_frame: bool
    phi_unit: str = ""
    theta_unit: str = ""
    tolerance: float = 0.05


def scan_plans(settings, args):
    keys = {"grid_based_scan": "grid", "polar_coordinate_based_scan": "polar"}
    configured = []
    if "scans" in settings:
        if any(key in settings for key in keys):
            raise ValueError("scans cannot be combined with top-level scan sections")
        entries = settings["scans"]
        if not isinstance(entries, list) or not entries:
            raise ValueError("scans must be a nonempty list")
        for entry in entries:
            if not isinstance(entry, dict) or len(entry) != 1 or next(iter(entry)) not in keys:
                raise ValueError("each scans entry must contain exactly one grid_based_scan or polar_coordinate_based_scan")
            key, value = next(iter(entry.items()))
            configured.append((keys[key], value))
    else:
        for key, value in settings.items():
            if key in keys:
                values = value if isinstance(value, list) else [value]
                if not values:
                    raise ValueError(f"{key} must not be an empty list")
                configured.extend((keys[key], item) for item in values)

    grid_cli_names = ("alignment_mode", "distance", "scan_start", "scan_end", "include_start",
                      "include_end", "scan_step", "scan_count")
    grid_cli = any(getattr(args, prefix + name) is not None
                   for prefix in ("", "x_", "z_") for name in grid_cli_names)
    polar_cli = any(value is not None for name, value in vars(args).items() if name.startswith("polar_"))
    systems = {system for system, _ in configured}
    if not configured:
        if grid_cli and polar_cli:
            raise ValueError("grid and polar scan settings cannot be combined without configured scans")
        system = args.coordinate_system or ("polar" if polar_cli else "grid")
        configured = [(system, {})]
        systems = {system}
    if (grid_cli and "grid" not in systems) or (polar_cli and "polar" not in systems):
        raise ValueError("coordinate-system conflicts with configured scan settings")
    if args.coordinate_system and systems != {args.coordinate_system}:
        raise ValueError("coordinate-system conflicts with configured scan settings")

    result = []
    for system, config in configured:
        if system == "grid":
            if not isinstance(config, dict) or set(config) - {"x", "y", "z"}:
                raise ValueError("grid_based_scan must be a mapping containing only x, y and z")
            y, _ = axis_plan("y", config, args)
            x, x_requested = axis_plan("x", config, args)
            z, z_requested = axis_plan("z", config, args)
            result.append(ScanPlan(system, (y, x, z), x_requested or z_requested))
        else:
            radius, phi, theta, phi_unit, theta_unit, tolerance = polar_options(config, args)
            result.append(ScanPlan(system, (radius, phi, theta), True, phi_unit, theta_unit, tolerance))
    return result

