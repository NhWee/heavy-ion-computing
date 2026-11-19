# Module 4 — FastJet: 제트 재구성과 배경 차감

상태(2026-11-19): **완료** — 클라우드 컨테이너 (PYTHIA 8.318, PyPI `fastjet` 3.5.2).
WSL2 재현: `bash scripts/run_all_pythia.sh`.

---

## 0. 이 모듈이 답하는 질문

> 쿼크·글루온은 직접 보이지 않는다. 그 대신 보이는 "제트"는 어떻게 정의하고,
> 주변에 입자가 가득한 환경(중이온 충돌)에서는 그 정의가 어떻게 오염되는가?

## 1. 생성 — `generators/pythia_hardqcd_13TeV.cmnd`

Module 2 카드와 같은 튠·PDF·primary 정의·seed. 바뀐 것은 과정뿐:
`HardQCD:all = on`, `PhaseSpace:pTHatMin = 20 GeV`. σ_gen = **0.722 mb** (pT̂ > 20 GeV),
5,000 이벤트 (unweighted). pT̂ 컷 근처(< 40 GeV)의 제트는 ISR/FSR·MPI 때문에 pT̂와 어긋나므로
**편향 영역**으로 표시하고 결론에 쓰지 않는다.

```bash
python scripts/m02_run_pythia.py --card generators/pythia_hardqcd_13TeV.cmnd --out data/generated/m04_hardqcd_13TeV.root
python analysis/m04_jets.py
```

입력 입자: 중성미자를 뺀 모든 최종 상태 입자, |η| < 3. 제트: |y| < 2.

## 2. 알고리즘 — `src/hicpy/jets.py`

순차 재결합: d_ij = min(pT_i^2p, pT_j^2p) ΔR²/R², d_iB = pT_i^2p. 가장 작은 것이 d_ij면 합치고,
d_iB면 제트로 확정. 지수만 다르다:

| p | 알고리즘 | 먼저 합쳐지는 것 | 모양 |
|---|---|---|---|
| +1 | kT | 부드러운 입자 | 불규칙 |
| 0 | Cambridge/Aachen | 가까운 입자 (기하학적) | 불규칙 |
| −1 | anti-kT | 단단한 입자 | **반지름 R의 원** |

셋 다 IRC-safe(부드러운 입자 추가·공선 분할에 둔감) — 섭동 QCD로 계산 가능한 이유.

## 3. 결과 (`results/tables/m04/m04_jets_summary.json`, `results/figures/m04/`)

**(1) 제트 단면적** dσ/dpT dy (anti-kT, R = 0.2 / 0.4 / 0.6): 40 GeV 위 제트 수 216 / 454 / 698.
R이 클수록 같은 pT의 제트가 많다 — 큰 원이 out-of-cone 복사를 더 담는다.

**(2) 제트 면적 (R = 0.4, pT > 40 GeV)** — 로드맵 기준 ✅

| 알고리즘 | ⟨A⟩ | RMS | πR² ± 10 % 안 |
|---|---|---|---|
| anti-kT | 0.507 | 0.025 | **96 %** |
| C/A | 0.565 | 0.188 | 21 % |
| kT | 1.027 | 0.360 | 4 % |

πR² = 0.503. anti-kT은 원, kT은 주변의 부드러운 입자를 끌어모아 평균 2배 면적.
→ 중이온에서 anti-kT을 쓰는 이유: 배경이 들어오는 **면적이 예측 가능**하다.

**(3) dijet Δφ**: 선두 제트 > 40 GeV, 차선두 > 20 GeV인 이벤트의 54 %가 Δφ > 2.8 rad.
나머지는 3번째 제트(추가 복사) 때문에 back-to-back에서 벗어남. ⟨m_jet⟩(R = 0.4, pT > 40) = 9.8 GeV.

**(4) 배경 삽입(embedding)** — 로드맵 기준 ✅. 각 이벤트에 상관 없는 toy 배경
(dN/dη dφ = 30, ⟨pT⟩ = 0.5 GeV → ρ_true = **15 GeV / 단위 면적**)을 얹고, 선두 신호 제트
(> 40 GeV, 376개)와 ΔR < 0.2로 짝지은 제트를 비교:

| | 값 | 기대 |
|---|---|---|
| δpT 차감 전 | **7.24 ± 0.17 GeV** | ρ_true ⟨A⟩ = 7.57 GeV (−4 %) |
| 이벤트별 median ρ | 18.15 ± 1.90 | ρ_true + PYTHIA UE = 15 + 2.41 = 17.41 |
| δpT 차감 후 | **−0.71 ± 0.17 GeV** | 0 |
| 차감 후 68 % 반폭 | **2.57 GeV** | Poisson 요동 √(A · dN/dηdφ · ⟨pT²⟩) = 2.75 GeV |

해석:
- 배경은 제트 pT를 ρ·A만큼 올린다 — 4 % 안에서 예측대로.
- median ρ는 신호 이벤트 자체의 underlying event(MPI, ρ ≈ 2.4)도 잰다. 그래서 공정한 비교는
  **신호 제트에서도 자기 이벤트의 ρ·A를 뺀 것**과의 차이다 (처음 버전은 이를 빠뜨려 −2.7 GeV가
  나왔다).
- 그래도 남는 −0.7 GeV: median이 (배경 + 덩어리진 UE)의 밀도를 4 % 높게 잡는다 — 합의 median은
  median의 합이 아니다.
- **차감으로 없앨 수 없는 것은 요동이다**: 차감 후 폭이 순수 Poisson 예측과 맞는다. 실제 Pb–Pb
  (ρ ~ 100–250 GeV)에서는 이 폭이 10 GeV를 넘고, 그래서 저 pT 제트 측정이 어렵다.

### 정식 결과 (WSL2, PYTHIA 8.312, conda FastJet 3.5.1)
σ(pT̂ > 20) = 0.736 mb. 면적: anti-kT 0.504 (96 % within 10 %), C/A 0.548 (18 %), kT 1.027 (4 %).
dijet Δφ > 2.8: 59 %. embedding: δpT **7.33 ± 0.16 vs ρA 7.54 (−2.8 %)**; median ρ 18.01
(= 15 + UE 2.33 + 3.4 %); 차감 후 **−0.55 ± 0.16 GeV**, 68 % 반폭 **2.73 GeV vs Poisson 2.75** —
conda 바인딩에서도 같은 코드가 그대로 돌았다.

## 4. 한계와 다음

- toy 배경은 균일·무상관 — 실제 중이온 배경의 flow(v₂) 변조, 이벤트별 ρ 변화 없음.
  Module 8(Angantyr) 샘플로 교체 예정.
- 제트 단면적은 PYTHIA 예측일 뿐 데이터와 비교하지 않았다 (ALICE/ATLAS 13 TeV jet 표: Rivet 단계).
- FastJet 고스트 seed를 Python에서 고정할 수 없어 면적이 ~1e-3에서 흔들림 (오답노트 #20).
