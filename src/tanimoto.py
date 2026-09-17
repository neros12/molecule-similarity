# RDKit 스텁에 누락된 멤버 타입은 알려진 반환 타입으로 보완한다.
# pyright: reportUnknownMemberType=false

from typing import cast

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def tanimoto_similarity(
    smiles1: str,
    smiles2: str,
    *,
    radius: int = 2,
    fp_size: int = 2048,
    include_chirality: bool = True,
):
    if radius < 0:
        raise ValueError("radius must be nonnegative.")

    if fp_size <= 0:
        raise ValueError("fp_size must be positive.")

    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)

    for molecule in (mol1, mol2):
        molecule.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(molecule)
        if include_chirality:
            Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=include_chirality,
    )
    fingerprint_a = cast(DataStructs.ExplicitBitVect, generator.GetFingerprint(mol1))
    fingerprint_b = cast(DataStructs.ExplicitBitVect, generator.GetFingerprint(mol2))

    return DataStructs.TanimotoSimilarity(fingerprint_a, fingerprint_b)
