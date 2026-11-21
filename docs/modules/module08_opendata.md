# Module 8 — 공개 실험 데이터 직접 분석: CMS 이중뮤온 스펙트럼 (pp 8 TeV)

상태(2026-11-21): **코드·closure 완료, 데이터 실행 대기**. 클라우드 세션은 opendata.cern.ch, eospublic.cern.ch에
접속할 수 없어서(연결 실패), 데이터 단계는 WSL2에서: `bash scripts/run_all_pythia.sh m08`.

---

## 0. 이 모듈이 답하는 질문

> 지금까지는 실험이 **이미 분석해 HEPData에 올린 숫자**와 MC를 비교했다. 실험의 **충돌 사건 자체**를 받아서
> 관측량을 처음부터 만들 수 있는가?

대상: CMS 2012 pp 8 TeV 이중뮤온. 뮤온 쌍의 불변질량
m² = (E₁ + E₂)² − |p⃗₁ + p⃗₂|²
을 쌓으면 μ⁺μ⁻로 붕괴하는 입자들(ω/φ, J/ψ, ψ(2S), Υ(1S, 2S, 3S), Z)이 연속 배경 위의 봉우리로 나타난다.

**중이온 물리와의 연결**: 쿼코늄과 Z는 중이온의 표준 탐침이다.
- Υ(2S, 3S)는 Υ(1S)보다 약하게 묶여 있어 QGP에서 먼저 녹는다 — **순차 억제**. 측정량은 Pb–Pb와 pp의
  [Υ(nS)/Υ(1S)] 이중비. 이 모듈이 그 **pp 기준**을 실제 데이터로 만든다.
- Z는 강하게 상호작용하지 않으므로 Pb–Pb에서 N_coll(Module 6)에 비례해야 한다 — R_AA 정규화의 검증. 그 pp 기준.

## 1. 데이터와 출처

| | |
|---|---|
| 레코드 | CERN Open Data [12341](https://opendata.cern.ch/record/12341) — "DoubleMuParked dataset from 2012 in NanoAOD format reduced on muons" |
| DOI | 10.7483/OPENDATA.CMS.LVG5.QT81 (라이선스 CC0) |
| 내용 | 2012 Run B, C, ~6100만 사건, 뮤온 변수만 (nMuon, Muon_pt/eta/phi/mass/charge) |
| 분석 레코드 | [12342](https://opendata.cern.ch/record/12342) (DOI 10.7483/OPENDATA.CMS.AAR1.4NZQ) — 같은 선택의 CMS 교육용 분석 |
| 파일 | `Run2012BC_DoubleMuParked_Muons.root` (~2 GB, git에 넣지 않음) |

`scripts/m08_fetch_cms_dimuon.sh`가 받고 (HTTPS 이어받기 → 안 되면 xrdcp), sha256·크기·사건 수·출처 URL·시각을
`data/raw/cms_2012_doublemuparked.json`(git에 들어감)에 기록한다. CMS와 CERN은 이 데이터로 만든 결과를 보증하지 않는다
(레코드 이용 조건).

## 2. 분석 — `analysis/m08_dimuon.py`, `src/hicpy/dimuon.py`

1. **histogram** (WSL2, 데이터 필요): uproot로 100 MB씩 읽어 `nMuon == 2`, 반대 부호 쌍의 질량을 계산.
   두 선택: 전체, 그리고 두 뮤온 모두 |η| < 1.2 (barrel, 운동량 분해능이 가장 좋음). 히스토그램만 JSON으로
   (`results/tables/m08/m08_hist.json`, git에 들어감) → 맞춤 단계는 어디서나 돈다.
2. **fit**: 세 창에서 binned Poisson likelihood.
   - J/ψ, ψ(2S) (2.6–4.2 GeV, 10 MeV bin), Υ(1S, 2S, 3S) (8.5–11.5 GeV, 20 MeV bin): Crystal Ball (가우스 핵 +
     낮은 질량 쪽 거듭제곱 꼬리 — 뮤온의 광자 방출), 한 족은 꼬리 매개변수 공유, 폭은 질량에 비례.
   - Z (70–110 GeV, 0.25 GeV bin): Breit–Wigner(Γ_Z = PDG) ⊗ 가우스 = Voigt.
   - 배경: 지수 × 2차 다항식.
   - **질량**: 족마다 공통 **질량 척도** 하나, m_i = scale × m_i(PDG). 상태 사이의 간격은 PDG로 고정 (CMS의 방식).
     scale − 1이 검출기 운동량 척도의 어긋남.
3. 출력: 질량 척도, 분해능 σ/m, 수율, S/B, 그리고 **Υ(2S)/Υ(1S), Υ(3S)/Υ(1S), ψ(2S)/J/ψ** (원시 수율 비 —
   수용·효율 보정 없음; Pb–Pb와의 이중비에서는 대부분 상쇄).

## 3. Closure 검증 (데이터 없이, 클라우드)

`python analysis/m08_dimuon.py toy` — NanoAOD와 같은 가지(nMuon, Muon_pt, …)를 가진 장난감 파일: 알려진 비율의
공명(질량 분해능 1 %)을 무작위 pT·rapidity로 만들고 등방 붕괴시킨 뮤온 쌍 + 지수 연속 배경, 40만 사건.
같은 histogram → fit 경로를 거쳐 (`--tag _toy`, 그림 `results/figures/m08/m08_fits_toy.png`):

| | 되찾은 값 (전체) | 넣은 값 |
|---|---|---|
| J/ψ 수율 | 31908 ± 211 | 32000 |
| ψ(2S) 수율 | 1224 ± 92 | 1200 |
| Υ(1S) / Υ(2S) / Υ(3S) 수율 | 8049 ± 293 / 2371 ± 107 / 1312 ± 75 | 8000 / 2400 / 1200 |
| Υ(2S)/Υ(1S), Υ(3S)/Υ(1S) | 0.295 ± 0.017, 0.163 ± 0.011 | 0.30, 0.15 |
| 질량 척도 (J/ψ, Υ, Z) | 0.99992, 0.99994, 1.00025 ± 0.00038 | 1 |
| χ²/ndf (J/ψ, Υ, Z 창) | 1.26, 0.93, 0.98 | |

`tests/test_m08_dimuon.py`: 질량 재구성 (10⁻⁹), Crystal Ball 정규화, Υ 맞춤이 수율을 4σ 안에서·척도를 10⁻³ 안에서
되찾는지.

**closure가 잡아낸 것** (그대로 두었으면 데이터에서 틀렸을 것):
1. 세 Υ 질량을 따로 풀어 두면 맞춤이 **라벨을 바꿔** 1S 봉우리에 2S를 맞췄다 → 간격을 PDG로 고정한 공통 척도.
2. Z 창에서 배경 지수가 발산(배경이 거의 없을 때) → 배경 변수를 창 폭으로 정규화하고 범위 제한.
3. 그 상태에서 수치 헤시안이 특이해져 오차가 NaN/0 → 배경 모양을 최적값에 고정하고 나머지로 오차 계산.

## 4. 데이터 결과

(WSL2 실행 후 채움: 질량 척도, σ/m (전체 vs barrel), Υ 상태 분리, 수율 비.)

## 5. 다음

- Module 9: CMS Pb–Pb 5.02 TeV의 Υ 이중비·Z R_AA (HEPData)를 이 pp 기준 + Module 6 T_AA와 맞대기.
- PYTHIA로 같은 선택의 이중뮤온 스펙트럼을 만들어 Drell–Yan 연속 배경과 Z 모양 비교.
