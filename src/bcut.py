# RDKit 스텁에 누락된 멤버 타입은 알려진 반환 타입으로 보완한다.
# pyright: reportUnknownMemberType=false

from collections.abc import Sequence
from math import dist, isfinite
from typing import cast

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

BCUT_DESCRIPTOR_NAMES = (
    "mass_high",
    "mass_low",
    "charge_high",
    "charge_low",
    "logp_high",
    "logp_low",
    "mr_high",
    "mr_low",
)


def bcut_descriptors(smiles: str):
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    descriptors = cast(list[float], rdMolDescriptors.BCUT2D(molecule))

    return descriptors


def bcut_similarity(
    smiles1: str,
    smiles2: str,
    *,
    scales: Sequence[float] | None = None,
):
    if scales is None:
        scales = (1.0,) * 8

    if len(scales) != 8:
        raise ValueError("scales must contain eight values in BCUT descriptor order.")

    for scale in scales:
        if not isfinite(scale) or scale <= 0:
            raise ValueError("Every BCUT scale must be finite and positive.")

    descriptors_a = bcut_descriptors(smiles1)
    descriptors_b = bcut_descriptors(smiles2)
    scaled_a = [value / scale for value, scale in zip(descriptors_a, scales)]
    scaled_b = [value / scale for value, scale in zip(descriptors_b, scales)]

    return 1.0 / (1.0 + dist(scaled_a, scaled_b))
