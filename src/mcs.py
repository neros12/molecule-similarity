from typing import cast

from rdkit import Chem
from rdkit.Chem import rdFMCS


def mcs_similarity(
    smiles1: str,
    smiles2: str,
    *,
    timeout: int = 10,
    include_chirality: bool = False,
):
    if timeout <= 0:
        raise ValueError("timeout must be positive.")

    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)

    for molecule in (mol1, mol2):
        molecule.UpdatePropertyCache(strict=False)
        if include_chirality:
            Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)

    result = rdFMCS.FindMCS(
        [mol1, mol2],
        maximizeBonds=False,
        timeout=timeout,
        atomCompare=cast(rdFMCS.AtomCompare, rdFMCS.AtomCompare.CompareElements),
        bondCompare=cast(rdFMCS.BondCompare, rdFMCS.BondCompare.CompareOrderExact),
        ringMatchesRingOnly=True,
        matchValences=True,
        matchChiralTag=include_chirality,
    )
    canceled = cast(bool, result.canceled)
    if canceled:
        raise TimeoutError(f"MCS search did not finish within {timeout} seconds.")

    common_atoms = cast(int, result.numAtoms)
    union_atoms = mol1.GetNumAtoms() + mol2.GetNumAtoms() - common_atoms

    return common_atoms / union_atoms
