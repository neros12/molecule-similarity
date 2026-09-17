# Molecular similarity

두 SMILES의 분자 유사도를 RDKit으로 계산하는 네 가지 모듈이다. 분자 유사도는
**분자 표현 방법과 비교 함수의 조합**으로 정의된다. 예를 들어 Morgan은 원자 주변
구조를 지문으로 표현하는 방법이고, Tanimoto는 두 지문을 비교하는 계수다.
BCUT은 계수가 아니라 연속형 서술자다.
[RDKit 개요](https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity)

| 모듈 | 분자 표현 | 비교 방법 | 적합한 비교 |
| --- | --- | --- | --- |
| `tanimoto.py` | Morgan bit fingerprint | Tanimoto | 원자 주변의 국소 구조 |
| `bcut.py` | BCUT2D 8개 서술자 | 스케일을 적용한 유클리드 거리 → 유사도 | 연결성과 원자 물성 |
| `atom_pair.py` | Atom Pair count fingerprint | Dice | 원자쌍의 유형·결합 거리·빈도 |
| `mcs.py` | 최대 공통 부분구조(MCS) | 공통 원자 수 기반 Jaccard | 공통 골격 |

## 조사 내용과 계산 정의

### Morgan–Tanimoto

반지름 2, 길이 2048 bit의 Morgan 지문을 기본으로 사용한다. 반지름은 원자 주변에서
탐색할 결합 거리를 뜻한다. 켜진 bit 수가 각각 `a`, `b`, 공통 bit 수가 `c`일 때
`T = c / (a + b - c)`다. 구조 탐색의 출발점으로 사용하기 쉽지만, 유한한 지문에는
해시 충돌이 생길 수 있다. 반지름과 bit 수가 달라지면 점수도 달라진다.
기본값 `include_chirality=True`로 입체화학 정보를 반영한다.
구현은 `GetMorganGenerator`, `GetFingerprint`, `TanimotoSimilarity`를 호출한다.
[지문 API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html),
[유사도 API](https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html)

### BCUT

RDKit의 `BCUT2D`는 원자 질량, Gasteiger 전하, Crippen logP 기여도, Crippen MR
기여도로 가중한 행렬의 최대·최소 고유값을 계산한다. 반환 순서는 다음과 같다.

```text
mass_high, mass_low, charge_high, charge_low,
logp_high, logp_low, mr_high, mr_low
```

[RDKit의 BCUT 구현](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Descriptors/BCUT.cpp)

이 프로젝트는 두 서술자 벡터 `x`, `y`를 아래처럼 비교한다.

```text
d = sqrt(sum(((x[i] - y[i]) / scales[i]) ** 2))
similarity = 1 / (1 + d)
```

이는 **프로젝트에서 선택한 유사도 정의**이며, 표준화된 BCUT 유사도 계수가 아니다.
`scales=None`이면 모든 축에 1을 사용하여 원래 단위로 비교한다. 따라서 변화 폭이
큰 축이 거리를 지배할 수 있다. 데이터셋 비교에는 기준 데이터에서 구한 축별
표준편차 등 고정된 양수 8개를 `scales`로 전달할 수 있다. 모든 비교에 같은 값을
사용해야 하며, 비교할 두 분자만으로 매번 스케일을 다시 구하지 않는다.
표준편차가 0인 축은 사용 전에 별도 정책으로 양수 스케일을 정해야 한다.

`bcut_descriptors()`로 원시 벡터를 얻을 수 있다. 입체이성질체의 구별을 위한 지표는 아니다.

### Atom Pair–Dice

원자 유형과 두 원자 사이의 최단 결합 거리로 특징을 만들고 출현 횟수를 보존한다.
`GetAtomPairGenerator(use2D=True)`와 `GetSparseCountFingerprint`를 사용하므로
3D 좌표가 필요 없고 고정 길이 bit 지문으로 접지 않는다. 기본값은 거리 1–30,
`include_chirality=True`다.
[지문 API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html)

계수는 RDKit의 `DiceSimilarity`로 계산하며, 비음수 count 벡터에 대해
`D = 2 * sum(min(x[i], y[i])) / (sum(x) + sum(y))`다.
단원자 분자처럼 지문에 원자쌍이 없는 경우도 RDKit의 `DiceSimilarity`로 계산한다.
[Atom Pair 설명](https://www.rdkit.org/docs/GettingStartedInPython.html#atom-pairs-and-topological-torsions)

동일한 bit 지문에 적용한 Dice는 `D = 2T / (1 + T)`이므로 Tanimoto와 순위가 같다.
여기서는 서로 다른 정보를 비교하도록 Morgan bit 지문 대신 Atom Pair count
지문에 Dice를 적용했다. 이 관계식은 위 bit 계수 정의에서 도출된다.

### MCS

`FindMCS`로 공통 부분구조를 찾고, 공통 원자 수 `m`과 각 분자의 원자 수 `n_A`,
`n_B`로 `S = m / (n_A + n_B - m)`를 계산한다. 이 정규화 공식도 프로젝트의
선택이다. `maximizeBonds=False`로 점수에 맞게 원자 수를 최대화한다.
원소·원자가를 비교하고 결합 차수를 정확하게 구분하며 고리 결합은 고리 결합에만
대응시킨다. 완전한 고리만 일치하도록 제한하지는 않는다.

연결된 공통 부분구조를 찾으므로 다중 성분 SMILES는 동일한 입력끼리 비교해도 점수가 1이 아닐 수 있다.
기본값 `include_chirality=False`로 입체화학을 비교하지 않는다. 원자 수는 RDKit
SMILES 파싱 후 남은 원자 기준이다. 복잡한 분자에서는 지문 비교보다 비용이 클 수
있다. 기본 제한은 10초이며, RDKit이 탐색 미완료를 보고하면 `TimeoutError`를
발생시켜 부분 해를 확정된 최대값으로 반환하지 않는다.
[MCS 설명](https://www.rdkit.org/docs/GettingStartedInPython.html#maximum-common-substructure),
[MCS API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html)

## 설치와 사용

프로젝트의 `.python-version`에 맞춰 Python 3.14 이상을 사용한다.

```bash
python -m pip install -e ".[dev]"
```

```python
from atom_pair import atom_pair_similarity
from bcut import bcut_similarity
from mcs import mcs_similarity
from tanimoto import tanimoto_similarity

smiles1 = "CCO"
smiles2 = "CCN"

print(tanimoto_similarity(smiles1, smiles2))
print(bcut_similarity(smiles1, smiles2))
print(atom_pair_similarity(smiles1, smiles2))
print(mcs_similarity(smiles1, smiles2))
```

각 모듈에서 개별 함수를 직접 import해도 된다.

```python
from bcut import BCUT_DESCRIPTOR_NAMES, bcut_descriptors
from tanimoto import tanimoto_similarity

print(dict(zip(BCUT_DESCRIPTOR_NAMES, bcut_descriptors("CCO"))))
print(tanimoto_similarity("CCO", "CCN", radius=3, fp_size=4096))
```

| 함수 | 선택 인수 |
| --- | --- |
| `tanimoto_similarity(smiles1, smiles2)` | `radius=2`, `fp_size=2048`, `include_chirality=True` |
| `bcut_similarity(smiles1, smiles2)` | `scales=None` |
| `atom_pair_similarity(smiles1, smiles2)` | `include_chirality=True` |
| `mcs_similarity(smiles1, smiles2)` | `timeout=10`, `include_chirality=False` |

유한한 유사도 점수는 0–1 범위의 실수이며 값이 클수록 유사하다.
BCUT의 수학적 범위는 유한한 거리에 대해 `(0, 1]`다. SMILES 타입·구문 오류나
원자가 없는 분자는 별도로 검사하지 않으며, RDKit이나 후속 연산의 예외를 그대로 전달한다.

`Chem.MolFromSmiles(smiles, sanitize=False)`로 직접 파싱하며, 별도 파싱 헬퍼는 사용하지
않는다. Morgan·Atom Pair·MCS에는 필요한 property cache를 `strict=False`로 초기화하고,
Morgan에는 고리 정보도 초기화한다. 입체화학 비교가 켜져 있으면 입체화학 정보를 할당한다.
`Chem.SanitizeMol()`은 호출하지 않으므로 파싱 시 원자가 검증·방향족성 인식 등 자동
sanitization을 수행하지 않는다. BCUT 계산 내부의 RDKit 처리는 그대로 사용한다.
염 제거, tautomer 통일, protonation 통일 등 추가 표준화는 호출 전에 적용해야 한다. 각 방법은 파싱 가능한
다중 성분 SMILES도 전체를 계산에 사용한다. 점수 1은 해당 표현에서 차이를 찾지 못했다는 뜻이며 동일한
분자나 생물학적 활성을 보장하지 않는다. 서로 다른 방법의 점수와 임계값은 직접
동등하게 해석하지 않는다. 3D 형상 비교는 conformer 생성·정렬이 추가로 필요하므로
이번 네 모듈의 범위에서 제외했다.

## 검증

```bash
python -m pyright
python -m ruff check src
python -m ruff format --check src
```

검증 환경: Python 3.14.7, RDKit 2026.3.6. pytest는 실행하지 않는다.
Pyright 검사 결과 파일은 보관하지 않으며, 생성한 경우 검사 후 삭제한다.
RDKit 스텁에 누락된 반환 타입은 `typing.cast()`로 보완한다. 각 모듈은 불완전한
외부 멤버 타입 진단만 제외하며, 구체적인 반환값과 프로젝트 코드는 standard Pyright로
검사한다. 별도 `rdkit-stubs` 패키지는 필요 없다.
