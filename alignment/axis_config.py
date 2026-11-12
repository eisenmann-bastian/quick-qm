"""Names and configuration for alignment directions."""

from pathlib import Path

import yaml


AXIS_NAMES = {
    "aromatic_normal": "aromatic_ring_normal",
    "carboxylate_sigma": "carboxylate_bond_extension",
    "carboxylate_lone_pair": "carboxylate_lone_pair",
    "chalcogen_sigma": "chalcogen_bond_extension",
    "chalcogen_lone_pair": "chalcogen_lone_pair",
    "halogen_sigma": "halogen_bond_extension",
    "hydrogen_sigma": "hydrogen_donor_extension",
    "nitrogen_triple_sigma": "nitrile_nitrogen_extension",
    "nitrogen_ring_lone_pair": "ring_nitrogen_lone_pair",
    "nitrogen_imine_lone_pair": "imine_nitrogen_lone_pair",
    "thioether_lone_pair_normal": "thioether_plane_normal",
    "thioether_outward_bisector": "thioether_outward_bisector",
    "thiol_lone_pair_normal": "thiol_plane_normal",
    "thiol_outward_bisector": "thiol_outward_bisector",
}
DEFAULT_CONFIG = Path(__file__).with_name("axes.config")


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject duplicate keys instead of silently overriding axis settings."""

    def construct_mapping(self, node, deep=False):
        mapping = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, (str, int)) or key in mapping:
                raise ValueError(f"invalid or duplicate YAML key: {key!r}")
            mapping[key] = self.construct_object(value_node, deep=deep)
        return mapping


class AxisConfig:
    """Validate axis settings and resolve individual direction switches."""

    def __init__(self, definitions):
        self.definitions = definitions

    @classmethod
    def load(cls, path=DEFAULT_CONFIG):
        try:
            with Path(path).open(encoding="utf-8") as source:
                document = yaml.load(source, Loader=UniqueKeyLoader)
        except yaml.YAMLError as exc:
            raise ValueError(f"invalid axis configuration {path}: {exc}") from exc
        if not isinstance(document, dict) or set(document) != {"axes"} or not isinstance(document["axes"], dict):
            raise ValueError("axis configuration must contain only an 'axes' mapping")
        definitions = document["axes"]
        if any(not isinstance(name, str) for name in definitions):
            raise ValueError("axis configuration names must be strings")
        expected = set(AXIS_NAMES.values())
        missing, extra = expected - definitions.keys(), definitions.keys() - expected
        if missing or extra:
            raise ValueError(f"axis configuration names mismatch: missing={sorted(missing)}, unknown={sorted(extra)}")
        for name, options in definitions.items():
            if not isinstance(options, dict) or set(options) - {"enabled", "directions"} or "enabled" not in options:
                raise ValueError(f"axis {name!r} must specify enabled and optional directions")
            if type(options["enabled"]) is not bool:
                raise ValueError(f"axis {name!r} enabled must be a boolean")
            directions = options.get("directions", {})
            if not isinstance(directions, dict) or any(
                type(index) is not int or index < 1 or type(value) is not bool
                for index, value in directions.items()
            ):
                raise ValueError(f"axis {name!r} directions must map positive integers to booleans")
        return cls(definitions)

    def enabled(self, name, direction):
        options = self.definitions[name]
        return options["enabled"] and options.get("directions", {}).get(direction, True)
