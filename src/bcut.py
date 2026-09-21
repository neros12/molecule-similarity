from math import dist
from typing import cast

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def _bcut_descriptors(smiles: str):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    descriptors = cast(list[float], rdMolDescriptors.BCUT2D(mol))

    return descriptors


def bcut_similarity(smiles1: str, smiles2: str):
    descriptors_a = _bcut_descriptors(smiles1)
    descriptors_b = _bcut_descriptors(smiles2)

    return 1.0 / (1.0 + dist(descriptors_a, descriptors_b))
