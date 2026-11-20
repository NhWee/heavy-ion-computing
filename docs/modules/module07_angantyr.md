# Module 7 — PYTHIA 8 Angantyr: p–Pb와 Pb–Pb를 핵자–핵자 충돌의 합으로

상태(2026-11-21): **1차 완료** — 클라우드 PYTHIA 8.318 (pp 20만, p–Pb 20만, Pb–Pb 4000 사건).
WSL2 정식 재현: `bash scripts/run_all_pythia.sh m07` (pp 100만, p–Pb 100만, Pb–Pb 2만; 16코어 ~25분).
남은 것: Pb–Pb 중심도별 dN_ch/dη 전체 (HEPData ins1410589 Table 1 — 클라우드에서 hepdata.net 접속 불가).

---

## 0. 이 모듈이 답하는 질문

> 중이온 충돌을 **집단 효과·QGP 없이** 핵자–핵자 충돌들의 합으로만 만들면 얼마나 맞고, 어디서 틀리는가?

Angantyr(Bierlich, Gustafson, Lönnblad, Shah, arXiv:1806.10820)는 그런 기준선이다:
1. 요동하는 핵자 크기·단면적을 가진 **Glauber MC**가 어떤 핵자가 어떻게(흡수적 = 비회절, 또는 회절적) 상처 입는지
   정한다 — Module 6의 black-disk Glauber와 비교할 수 있다;
2. 흡수적 하위충돌 하나하나가 보통의 PYTHIA pp 사건(MPI, 샤워)이 되고, 같은 핵자의 추가 충돌은 회절 같은
   들뜸으로 붙는다 (Fritiof의 wounded-nucleon 그림);
3. 전체를 Lund 끈으로 한꺼번에 강입자화한다.

세 질문:
- **Glauber**: Angantyr의 N_part, N_coll, σ_inel은 Module 6과 측정값에 대해 어떤가?
- **다중도**: (2/⟨N_part⟩) dN_ch/dη를 pp → p–Pb → Pb–Pb에서.
- **기묘도**: Module 2에서 PYTHIA pp의 기묘 입자 비가 다중도에 평평했다. 핵자 충돌을 쌓으면 달라지는가?

## 1. 도구

- 카드 `generators/angantyr_{pPb,PbPb}_5TeV.cmnd` — √s_NN = 5.02 TeV, `HeavyIon:mode = 1`, 기본 하위충돌 모형
  (`Angantyr:CollisionModel = 2`, 초기화 때 PYTHIA의 pp 단면적에 유전 알고리즘으로 맞춤, ~30 s),
  `HeavyIon:forceUnitWeight = on` (아래 함정 1), ALICE 일차 입자 정의 (cτ < 1 cm만 붕괴).
- `analysis/m07_angantyr.py generate --system {pp,pPb,PbPb}` — 병렬 생성, 사건마다 Glauber 정보
  (`pythia.getHeavyIonsPtr().hiInfo`: b, 상처 입은 핵자, 하위충돌 수)와 ALICE 관측 정의:
  - p–Pb는 **ALICE 실험실계**: 양성자가 −η로 가고 NN 중심계가 rapidity −0.465로 움직인다 → η_lab = −η_NN − 0.465.
  - V0A 2.8 < η_lab < 5.1, V0C −3.7 < η_lab < −1.7의 하전입자 수. Pb–Pb 중심도 = V0A + V0C (V0M),
    p–Pb = V0A (Pb 쪽), p–Pb "NSD" = V0A와 V0C 양쪽에 하전입자 (ALICE 트리거 흉내).
  - 동정·기묘 입자는 |y_lab| < 0.5 (Module 2와 같은 정의).
- `analysis/m07_angantyr.py compare` — 그림 4장과 `results/tables/m07/m07_summary.json`.

**함정 두 개 (troubleshooting #26)**
1. Angantyr는 기본적으로 b를 가우스 분포에서 뽑고 **사건 가중치**를 준다. 첫 분석은 가중치를 무시해서 같은
   중심도 클래스의 ⟨b⟩가 8 % 작게 나왔다 (σ 차이로는 1.5 %밖에 설명 안 됨). `forceUnitWeight = on`으로 고친 뒤
   ⟨b⟩는 Module 6과 대부분 클래스에서 2 % 안, 최대 4 % 차이로 같다 (30–40 %: 9.21 vs 9.20 fm; 0–1 %는 40사건뿐). 분석은 가중치가 1이 아니면 멈춘다.
2. `HIInfo.nCollTot()`은 탄성 하위충돌을, `nPartProj()`는 탄성으로만 산란된 핵자를 포함한다. Glauber 관례에 맞춰
   N_part = 흡수 + 회절로 상처 입은 핵자, N_coll = 전체 − 탄성으로 계산한다 (구성 요소의 합으로 확인).

## 2. 결과 (클라우드 8.318)

### 2.1 Glauber: Angantyr vs Module 6 (`results/figures/m07/m07_glauber_vs_m06.png`)

| | Angantyr | Module 6 (black disk) | 측정 |
|---|---|---|---|
| σ_PbPb | 7.21 ± 0.11 b | 7.64 ± 0.03 b | 7.67 ± 0.25 b (ALICE, arXiv:2204.10148) |
| σ_pPb | 2.002 ± 0.004 b | 2.10 ± 0.03 b | 2.061 ± 0.080 b (CMS, arXiv:1509.03893) |

| 중심도 (b 순서) | ⟨N_part⟩ Angantyr / Module 6 | ⟨N_coll⟩ Angantyr / Module 6 | ⟨b⟩ (fm) Angantyr / Module 6 |
|---|---|---|---|
| 4–5 % | 355 / 362 | 1768 / 1596 | 3.27 / 3.33 |
| 10–20 % | 249 / 266 | 1060 / 990 | 6.02 / 5.95 |
| 30–40 % | 117 / 132 | 362 / 342 | 9.21 / 9.20 |
| 50–60 % | 48 / 54 | 99 / 89 | 11.38 / 11.58 |
| 70–80 % | 13.4 / 15.7 | 17.9 / 16.2 | 13.29 / 13.53 |

- 같은 b에서 Angantyr는 **N_part가 2–15 % 적고 N_coll이 6–11 % 많다** (주변부로 갈수록 N_part 차이가 큼). 요동하는 핵자 단면적: 큰 단면적을 가진
  핵자는 여러 번 충돌하고(N_coll ↑), 작은 단면적을 가진 핵자는 빠진다(N_part ↓). black disk는 모든 쌍이 같은 σ.
- σ_PbPb는 Angantyr가 측정보다 6 % 낮다 (1.8σ). σ_pPb는 둘 다 측정과 맞다.
- V0M 순서와 b 순서의 N_part가 60 %까지 1 % 안에서 같다 → Pb–Pb에서 다중도 중심도는 기하를 잘 따른다.
  그 뒤로 갈라져 80–100 %에서는 10–30 % 차이 (사건마다 요동이 큰 영역).

### 2.2 다중도 (`m07_dndeta_per_participant.png`, `m07_dndeta_vs_eta.png`)

| | Angantyr | ALICE |
|---|---|---|
| p–Pb NSD, dN_ch/dη (\|η_lab\| < 0.5) | 18.40 (N_part 7.40) | 16.81 ± 0.71 (arXiv:1210.3615) → MC/data 1.09 |
| p–Pb NSD, 참여자당 | 2.49 | 2.14 ± 0.17 |
| Pb–Pb 0–5 %, dN_ch/dη (\|η\| < 0.5) | 1805 (N_part 376) | 1943 ± 54 (arXiv:1512.06104) → MC/data 0.93 |
| pp 5.02 TeV INEL>0 | 5.47 | — |

- **집단 효과 없이도 중심 Pb–Pb의 다중도를 7 % 안에서 만든다.** (2/N_part) dN_ch/dη가 N_part와 함께 2 → 9.6으로
  오르는 것은 N_coll이 N_part보다 빨리 자라서(MPI가 핵자 충돌 수를 따라감) — 참여자 스케일링이 깨지는 이유가
  "QGP"일 필요는 없다는 것을 이 기준선이 보여 준다.
- p–Pb에서는 V0A 다중도 클래스의 (2/N_part) dN_ch/dη가 N_part에 대해 단조롭지 않다 (0–5 %에서 다시 내려감).
  p–Pb에서 다중도로 고른 중심도는 N_part와 약하게만 연결된다 (클래스 0–5 %와 60–80 %의 ⟨b⟩가 3.1과 6.2 fm) —
  고다중도 클래스는 "N_part가 많은 사건"이 아니라 "핵자 충돌마다 다중도가 큰 사건"을 많이 고른다. ALICE가 p–Pb
  중심도에 별도 방법(하이브리드, 영도 열량계)을 쓰는 이유.

### 2.3 기묘도: pp → p–Pb → Pb–Pb (`m07_strangeness_vs_multiplicity.png`)

| | ⟨dN_ch/dη⟩ | K/π | Λ/π | Ξ/π | Ω/π |
|---|---|---|---|---|---|
| PYTHIA pp 5.02 (INEL>0) | 5.5 | 0.128 | 0.022 | 0.0020 | 1 × 10⁻⁴ |
| Angantyr p–Pb 0–5 % | 47 | 0.124 | 0.021 | 0.0019 | 1 × 10⁻⁴ |
| Angantyr Pb–Pb 0–5 % | 1805 | 0.128 | 0.022 | 0.0019 | 1 × 10⁻⁴ |
| ALICE pp 13 TeV 최고 다중도 | 26 | 0.130 | 0.041 | 0.0055 | 5 × 10⁻⁴ |

- **dN_ch/dη가 2에서 1800까지 늘어도 Angantyr의 기묘 입자 비는 그대로다.** Angantyr의 끈 조각화는 끈마다 독립이라,
  끈을 많이 쌓아도 s 쿼크 확률은 같다. 데이터는 pp 안에서만 이미 Λ/π가 1.5배, Ω/π가 3배 오른다 (Module 2 §14) —
  이 증가가 "큰 계의 성질"이 아니라 **다중도의 성질**이라는 ALICE의 관찰을 Angantyr는 만들지 못한다.
- 같은 이유로 무거운 입자 ⟨pT⟩가 다중도와 함께 오르지 않는다: Pb–Pb 0–5 %의 ⟨pT⟩_p 0.88 GeV/c가 pp 고다중도
  클래스(0.93)보다 낮다. radial flow가 없는 모형의 표시.

## 3. 결론과 다음

1. Angantyr는 **기하(Glauber)와 전체 다중도**의 좋은 기준선이다: σ_pPb, 중심 Pb–Pb dN_ch/dη를 10 % 안에서,
   b–중심도 관계는 Module 6과 대부분 2 % 안에서. 요동하는 핵자 단면적 때문에 N_part는 적고 N_coll은 많다 — R_AA의
   N_coll 정규화가 Glauber 모형 선택에 5–10 % 의존한다는 뜻.
2. **기묘도 증가와 radial flow는 없다** — 이것이 "핵자 충돌의 합"과 실제 중이온 충돌의 차이를 재는 지표다.
   PYTHIA 안에서 이를 더하는 시도는 rope(Module 2 §15에서 pp로 시험: 증가분의 30–60 %)와 shoving.
3. 다음: Pb–Pb 중심도별 dN_ch/dη 전체(ins1410589), Angantyr + rope/shoving의 Pb–Pb 시험, 그리고 Module 8의
   공개 데이터 재현.

## 4. 실행

```bash
python analysis/m07_angantyr.py generate --system pp   --n-events 200000
python analysis/m07_angantyr.py generate --system pPb  --n-events 200000
python analysis/m07_angantyr.py generate --system PbPb --n-events 4000
python analysis/m07_angantyr.py compare
```
