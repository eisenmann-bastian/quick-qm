"""Molecule parsing and functional-group matching."""

import io
from pathlib import Path

from rdkit import Chem

from .definitions import SelectorRegistry


class MoleculeSelector:
    """Apply a named selector or SMARTS to all records in an SDF."""

    def __init__(self, registry=None):
        self.registry = registry if registry is not None else SelectorRegistry()

    @staticmethod
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

    def select(self, sdf_file, *, name=None, smarts=None):
        """Return one entry per SDF record with zero-based atom-index matches."""
        if (name is None) == (smarts is None):
            raise ValueError("specify exactly one of name or smarts")
        if name is not None:
            selector_name, definition = self.registry.resolve(name)
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
        for index, mol in enumerate(self._read_sdf(path)):
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
