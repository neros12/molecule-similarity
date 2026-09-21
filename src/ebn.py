from typing import cast

import numpy as np
from rdkit import Chem

FloatArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]
FloatArray2D = np.ndarray[tuple[int, int], np.dtype[np.float64]]


__all__ = ["ebn_similarity"]


_METAL_ATOMIC_NUMBERS = {
    3,
    4,
    11,
    12,
    13,
    *range(19, 32),
    *range(37, 51),
    *range(55, 85),
    *range(87, 117),
}


def _build_molecular_matrix(mol: Chem.Mol):
    elements = ["C", "N", "O", "Si", "P", "S", "F", "Cl", "Br", "I"]
    bond_types = [
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
        Chem.BondType.AROMATIC,
    ]

    # 열: 원소(0:10), 비고리 S/D/T(10:13), CP(13), aryl(14),
    # 고리 S/D/T(15:18), 기타 AROMATIC(18).
    element_columns = {element: index for index, element in enumerate(elements)}
    bond_columns = {bond_type: index for index, bond_type in enumerate(bond_types[:3])}
    cp_pattern = Chem.MolFromSmarts("C1[*]C[*]C1")
    aryl_pattern = Chem.MolFromSmarts("C1[*]C[*]C[*]1")
    special_bond_columns: dict[int, int] = {}
    for pattern, column in ((cp_pattern, 13), (aryl_pattern, 14)):
        for matched in mol.GetSubstructMatches(pattern, maxMatches=0):
            for pattern_bond in pattern.GetBonds():
                begin = matched[pattern_bond.GetBeginAtomIdx()]
                end = matched[pattern_bond.GetEndAtomIdx()]
                matched_bond = mol.GetBondBetweenAtoms(begin, end)
                if matched_bond is not None:
                    special_bond_columns.setdefault(matched_bond.GetIdx(), column)

    molecular_matrix: FloatArray2D = np.zeros((mol.GetNumAtoms(), 19), dtype=np.float64)
    for atom in mol.GetAtoms():
        atom_index = atom.GetIdx()
        element_column = element_columns.get(atom.GetSymbol())
        if element_column is not None:
            molecular_matrix[atom_index, element_column] = 1.0

        for bond in cast(tuple[Chem.Bond, ...], atom.GetBonds()):
            bond_column = special_bond_columns.get(bond.GetIdx())
            if bond_column is None:
                bond_type = bond.GetBondType()
                if bond_type == Chem.BondType.AROMATIC:
                    bond_column = 18
                else:
                    bond_offset = 15 if atom.IsInRing() else 10
                    bond_column = bond_offset + bond_columns.get(bond_type, 0)
            molecular_matrix[atom_index, bond_column] += 1.0

    return molecular_matrix


def _get_ebn_vector(mol: Chem.Mol) -> FloatArray1D:
    adjacency_matrix: FloatArray2D = Chem.GetAdjacencyMatrix(mol)

    molecular_matrix = _build_molecular_matrix(mol)

    v1 = molecular_matrix.sum(axis=0)
    neighbor_matrix = molecular_matrix.copy()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in _METAL_ATOMIC_NUMBERS:
            neighbor_matrix[atom.GetIdx(), 10:] = 0.0

    v2 = (adjacency_matrix @ neighbor_matrix).sum(axis=0)
    v2[10:] -= neighbor_matrix[:, 10:].sum(axis=0)

    return np.concat([v1, v2])


def ebn_similarity(smiles1: str, smiles2: str):
    """
    Element-Bond Neighborhood Similarity

    두 분자의 원소와 결합 이웃 정보를 기반으로 EBN 유사도를 계산하세요.

     EBN (Element-Bond Neighborhood)은 분자 전체의 원소 및 결합 구성과
     인접 원자의 특징을 함께 비교합니다. 이 구현의 계산 과정은 다음과
     같습니다.

    Args:
        smiles1: 비교할 첫 번째 분자.
        smiles2: 비교할 두 번째 분자.

    Returns:
        float: 두 EBN 벡터의 코사인 유사도를 세제곱한 값.
            어느 한 벡터의 노름이 0이면 0.0을 반환합니다.
    """
    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
    mol1.UpdatePropertyCache(strict=False)
    mol2.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol1)
    Chem.FastFindRings(mol2)

    vector1 = _get_ebn_vector(mol1)
    vector2 = _get_ebn_vector(mol2)

    norm_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    if norm_product == 0.0:
        return 0.0

    cosine_similarity = float(np.dot(vector1, vector2) / norm_product)

    return cosine_similarity**3
