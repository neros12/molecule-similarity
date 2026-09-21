from typing import cast

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def atom_pair_similarity(smiles1: str, smiles2: str):
    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
    mol1.UpdatePropertyCache(strict=False)
    mol2.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol1)
    Chem.FastFindRings(mol2)

    generator = rdFingerprintGenerator.GetAtomPairGenerator(use2D=True)
    fp1 = cast(
        DataStructs.ULongSparseIntVect,
        generator.GetSparseCountFingerprint(mol1),
    )
    fp2 = cast(
        DataStructs.ULongSparseIntVect,
        generator.GetSparseCountFingerprint(mol2),
    )

    return DataStructs.DiceSimilarity(fp1, fp2)
