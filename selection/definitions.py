"""Named functional-group and interaction-site selector definitions."""

import re
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Lipinski


@dataclass(frozen=True)
class NamedSelector:
    smarts: str | None
    aliases: tuple[str, ...] = ()


CARBOXYLATE_SMARTS = "[C;+0;X3;$([C]-[#6]),$([C;H1])](=[O;+0;X1;H0:1])-[O-;X1;H0:2]"
CARBOXYLATE_QUERY = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)


def carboxylate_sites(mol):
    """Map each oxygen in a carboxylate group to its carbon and oxygen indices."""
    groups = {}
    for carbon, oxygen1, oxygen2 in mol.GetSubstructMatches(CARBOXYLATE_QUERY, maxMatches=0):
        group = (carbon, *sorted((oxygen1, oxygen2)))
        groups[oxygen1] = group
        groups[oxygen2] = group
    return groups


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
    "carboxylate": NamedSelector(CARBOXYLATE_SMARTS, ("carboxylate anion",)),
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


class SelectorRegistry:
    """Resolve canonical selector names and aliases."""

    def __init__(self, definitions=None):
        self.definitions = SELECTORS if definitions is None else definitions

    @staticmethod
    def normalize(name):
        return " ".join(re.sub(r"[-_]+", " ", name).casefold().split())

    def resolve(self, name):
        normalized = self.normalize(name)
        for canonical, definition in self.definitions.items():
            if normalized in (self.normalize(canonical), *map(self.normalize, definition.aliases)):
                return canonical, definition
        raise ValueError(f"unknown selector {name!r}; use --list-selectors to see available names")
