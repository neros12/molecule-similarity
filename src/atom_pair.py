from typing import cast

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def atom_pair_similarity(smiles1: str, smiles2: str):
    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)

    for molecule in (mol1, mol2):
        molecule.UpdatePropertyCache(strict=False)

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


if __name__ in "__main__":
    print(atom_pair_similarity("CCCC", "CCC"))
