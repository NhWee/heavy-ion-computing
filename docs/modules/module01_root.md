# Module 1 — ROOT 파일과 첫 물리 히스토그램

첫 실행 2026-08-27 · 마지막 재생성 2026-11-09 (클라우드 컨테이너, uproot arm ·
closure · fit 방법 검증 · 체계 오차) · ROOT arm은 이 저장소 계열에서 **미실행**
(WSL2 대기, §5)

---

## 0. 이 모듈이 답하는 질문

> collider 이벤트 하나는 컴퓨터 안에서 어떤 모양의 자료구조이고,
> 거기서 pT·η·다중도·불변질량은 어떻게 나오는가?

물리는 아직 QCD가 아닙니다. **자료구조와 분석 체인**이 주제입니다.

---

## 1. 왜 CSV가 아니라 ROOT인가

collider 이벤트는 **ragged(들쭉날쭉)** 합니다. 1번 이벤트에 137개, 2번에 61개의
입자가 있습니다. 평평한 표(CSV, DataFrame)는 이걸 담지 못합니다. 담으려면

- 입자마다 한 줄씩 쓰고 `event_id`를 붙이거나 (이벤트 단위 연산이 매번 group-by),
- 최대 입자 수만큼 열을 만들어 대부분을 NaN으로 채우거나 (용량 폭발)

해야 합니다. ROOT의 **TTree**는 branch마다 가변길이 배열을 저장하고, 이벤트 수만큼
반복하면서 필요한 branch만 디스크에서 읽습니다(*column-wise, lazy*). 20,000 이벤트
119만 입자가 19 MB입니다.

이번 모듈에서는 ROOT 설치 없이 순수 Python으로 ROOT 파일을 쓰는 **uproot**를
씁니다. 같은 분석의 ROOT C++ / PyROOT 버전도 함께 두어 두 생태계를 비교합니다.

### 만들어진 파일 구조

```
data/generated/m01_toy_events.root
├── events (TTree, 20000 entries)
│   ├── nparticle      int32          이벤트당 입자 수 (counter branch)
│   ├── particle_px    float32[]      ← nparticle 길이의 가변 배열
│   ├── particle_py    float32[]
│   ├── particle_pz    float32[]
│   ├── particle_E     float32[]
│   ├── particle_pdg   int32[]        PDG Monte Carlo ID
│   └── has_z          bool           스칼라 branch
└── provenance (TObjString)           생성기 설정·시드·git 커밋 JSON
```

`nMuon / Muon_pt / Muon_eta …` 형태의 CMS NanoAOD와 같은 관례입니다.
`counter branch`가 있어야 ROOT가 각 이벤트에서 배열을 몇 개 읽어야 할지 압니다.

---

## 2. 이벤트 생성기 (toy MC) — 무엇을 가정했나

`src/hicpy/toy_generator.py`. PYTHIA가 아직 없으므로, **실험이 실제로 데이터에
fit하는 함수 형태**를 그대로 써서 그럴듯한 샘플을 만듭니다.

| 성분 | 함수 형태 | 물리적 의미 |
|---|---|---|
| 다중도 $N_{\rm ch}$ | Negative Binomial (NBD), $\bar n=60,\ k=2.2$ | UA5 이래 하드론 충돌 다중도의 표준 기술. 분산 $=\bar n + \bar n^2/k$ 로 Poisson보다 넓음 → **이벤트별 요동** |
| $dN/d\eta$ | $\|\eta\|<2$ 평탄 + Gaussian 어깨($\sigma=1.6$) | 중앙 rapidity plateau(boost-invariant 영역) + fragmentation region |
| $p_{\rm T}$ | Tsallis/Hagedorn $\left[1+\frac{m_{\rm T}-m_0}{nT}\right]^{-n}$, $T=135$ MeV, $n=6.8$ | 저-$p_{\rm T}$는 지수(열적), 고-$p_{\rm T}$는 멱함수(hard scattering). 단일 지수함수로는 tail을 절대 못 맞춤 |
| $\phi$ | 균일 | flow도 jet도 없음 (Module 8·10의 주제) |
| 입자종 | $\pi:K:p = 0.83:0.12:0.05$ | LHC pp의 대략적 하전입자 조성 |
| Z→μ⁺μ⁻ | 2 % 이벤트에 relativistic Breit-Wigner ($m_Z=91.1876$, $\Gamma_Z=2.4952$ GeV) 후 등방 2체 붕괴 + boost | 불변질량 재구성 연습용. 입력값은 LEP 시대 PDG 값이고 현재 세계 평균(RPP 2026: 91.1879, 2.4955 GeV)과 다르지만, closure의 참값은 입력값입니다 (`hicpy.pdg`) |

**이 생성기가 아닌 것**: hard scattering도, parton shower도, hadronization도,
입자간 상관도 없습니다. 이벤트는 통계적으로만 그럴듯합니다. Module 2에서 PYTHIA로
교체하면 **분석 코드는 그대로 두고** 입력만 바뀌어야 합니다 — 그게 이 구조의 목적입니다.

전체 설정은 `config/analysis.yaml`에, 실제 사용된 값·시드·git 커밋은
`data/generated/m01_toy_events.json`과 ROOT 파일 내부 `provenance`에 기록됩니다.

---

## 3. 관측량과 결과

이벤트 20,000개, 시드 20260827, |η| < 0.8, pT > 0.15 GeV/c.

### 하전입자의 정의 (2026-11-05)

"하전입자"는 `particle_charge ≠ 0`입니다 — PDG 코드 목록이 아닙니다. toy에는
π/K/p와 Z의 뮤온뿐이라 두 정의의 차이는 뮤온 두 개(이벤트의 2 %)이지만,
PYTHIA 샘플(Module 2)에는 e±, Σ±, Ξ±, Ω± 등이 있어 목록 방식은 ALICE의 "charged
primary"보다 수 % 적게 셉니다. 정의를 입력의 우연이 아니라 물리(전하)에 걸어 두면
생성기를 바꿔도 분석은 그대로입니다. 이 변경으로 $\langle N_{\rm ch}\rangle$는
11.988 → 11.999, $\langle p_{\rm T}\rangle$는 0.551 → 0.584 GeV/c (고운동량 뮤온)로
바뀌었습니다. Closure test의 다중도 비교는 NBD에서 뽑은 **하드론만** 세도록
명시적으로 남겨 두었습니다(생성기 검증이 목적이므로).

### (1) 하전입자 다중도 분포 `results/figures/m01_multiplicity.pdf`

$\langle N_{\rm ch}\rangle = 11.99 \pm 0.06$ (RMS 8.78, |η|<0.8).
RMS/평균 ≈ 0.73 으로 Poisson(≈0.29)보다 훨씬 넓습니다 — NBD의 $k$가 만든
**이벤트별 요동**이고, 중이온 물리에서 centrality를 정의할 때 바로 이 폭을 씁니다.

### (2) $dN_{\rm ch}/d\eta$ `results/figures/m01_dndeta.pdf`

|η| < 0.5 평균으로 7.50. 실제 ALICE pp √s=13 TeV INEL>0 측정은 6.5 근처이므로 toy는 조금
높습니다 — 의도적으로 맞추지 않았습니다. Module 2에서 PYTHIA로 실제 데이터와
비교하는 것이 목표입니다.

### (3) 불변 $p_{\rm T}$ 스펙트럼 `results/figures/m01_pt_spectrum.pdf`

$$\frac{1}{2\pi p_{\rm T} N_{\rm ev}}\frac{d^2N}{dp_{\rm T}\,d\eta}$$

로 그렸습니다. $2\pi p_{\rm T}$로 나누는 이유: $d^2N/(dp_x dp_y) = d^2N/(2\pi p_T dp_T)$
이므로 이 조합이 **Lorentz 불변 수율**에 해당하고, 실험 논문이 이 형태로 발표합니다.

Tsallis fit 결과 (하전입자 전체, 파이온 질량 가정, fit 범위 0.15–5.66 GeV/c = 중심이
`fits.tsallis_pt_range` [0.15, 6.0] 안에 있는 bin 25개):

```
T = 147.6 ± 0.7 (stat) ± 3.5 (syst) MeV
n = 7.496 ± 0.057 (stat) ± 0.209 (syst)          T–n 상관계수 +0.90
-2 ln λ / ndf = 25.6 / 22                         (binned Poisson likelihood)
```

**fit 방법 (2026-11-07 변경).** 물리 결과는 원시 count에 대한 binned Poisson likelihood로
얻습니다. 모형은 각 bin에서 **적분**합니다(`hicpy.fitting.fit_poisson`). 이전 방법은
불변 수율 점들을 bin **중심**에서 최소제곱으로 맞췄고, 두 가지 근사가 있었습니다.

1. 로그 bin이 넓어서(폭 15 %) bin 평균과 bin 중심값이 고-pT에서 최대 4.6 % 다릅니다.
2. 오차를 데이터에서 잡으므로 count가 적은 bin일수록 가중이 왜곡됩니다.

pseudo-experiment 300회(§7)에서 이전 방법은 T를 평균 **+0.85σ**(π) 높게 주었고,
새 방법의 평균 pull은 0.00 ± 0.06입니다. 이전 방법(χ² = 23.9/22, T = 148.2 MeV)은
세 arm 교차검증용으로만 JSON의 `tsallis` 블록에 남아 있습니다. ROOT arm 두 개가
같은 알고리즘(TF1 χ²)을 돌리기 때문입니다.

**체계 오차**는 `analysis/m01_fit_systematics.py`가 fit 범위와 binning을 하나씩 바꿔 가며
구합니다(§7). 범위 의존성이 큰 것은 이 형태가 세 질량 혼합의 **참 모형이 아니기**
때문이고, 그 자체가 물리입니다(§4).

입력값 $T=135$ MeV, $n=6.8$ 과 **맞지 않습니다.** 이건 버그가 아니라 물리입니다
(§4 참조).

### (4) 이중뮤온 불변질량 `results/figures/m01_dimuon_mass.pdf`

$p_{\rm T}^\mu > 20$ GeV/c, $|\eta_\mu| < 2.4$, 정확히 2개 뮤온. 241 후보(창 60–120 GeV 안 240).

```
m_fit = 91.057 ± 0.105 GeV/c²      (생성기 입력 91.1876; RPP 2026 91.1879 ± 0.0020)
Γ_fit = 2.28  ± 0.23  GeV          (생성기 입력 2.4952;  RPP 2026 2.4955 ± 0.0023)
-2 ln λ / ndf = 53.2 / 57           binned Poisson likelihood, 상대론적 Breit–Wigner
```

배경 항은 넣지 않습니다. toy에는 이중뮤온 배경이 **없기** 때문입니다. 0 이상으로
제한한 배경 항을 넣으면 그 항이 BW 꼬리를 흡수해 Γ를 −0.06 GeV(−0.33σ) 편향시킵니다.
Module 2(Drell–Yan, heavy-flavour 뮤온)부터는 배경 항이 필요합니다.

**해석 정정 (2026-11-07).** Session 1–3은 "중심값이 PDG보다 1.3σ 낮고 폭이 좁은 것은
$p_{\rm T}^\mu>20$ GeV 컷이 낮은 질량 꼬리를 잘라내는 acceptance sculpting"이라고
설명했습니다. 이 설명은 **틀렸습니다**. 근거는 세 가지입니다(§7,
`results/tables/m01_fit_validation.json`).

1. 같은 표본의 Z를 **컷 없이** 모두(380개) fit해도 m = 91.142 ± 0.079, Γ = 2.14 ± 0.17
   GeV가 나옵니다. 낮은 질량과 좁은 폭은 컷 이전부터 이 표본에 있던 **통계적 요동**입니다.
2. acceptance는 질량이 커질수록 **증가**합니다(+0.23 %/GeV). 그래서 컷은 피크를 **위로**
   +0.010 ± 0.002 GeV 옮깁니다. 크기는 통계 오차의 1/10, 방향은 정반대입니다.
3. 좁은 폭의 나머지 원인은 **fit 방법**이었습니다. $\sigma=\max(\sqrt N,1)$인 χ²(Neyman)는
   240개짜리 피크에서 Γ를 평균 −0.18 GeV(−0.81σ) 편향시키고 질량 오차를 12 % 과소평가합니다.
   이 방법으로 얻은 이전 값은 m = 91.06 ± 0.09, Γ = 2.06 ± 0.23이었습니다.

교훈: 1–2σ 편차에는 이유를 붙이기 전에 그것이 요동인지부터 확인합니다. 단일 표본으로는
알 수 없으므로 truth 비교나 pseudo-experiment가 필요합니다(오답노트 #15).

---

## 4. Closure test — 분석이 정답을 되찾는가

`analysis/m01_closure_test.py` → `results/tables/m01_closure_test.json`.
생성기의 진짜 파라미터를 알고 있으므로 분석이 그것을 복원해야 합니다. 참값은 현재
`config/analysis.yaml`이 아니라 **입력 파일의 provenance**에서 읽습니다. 생성기 성분
네 가지를 모두 검사하고, 모든 pull이 4σ 안이어야 통과입니다.

```
observable                      measured             truth     pull
<N_ch> in acceptance        11.988 ± 0.062         12.059     -1.1
dN/deta, |eta|<0.5           7.502 ± 0.041          7.537     -0.9
NBD k (accepted N_ch)        2.206 ± 0.032          2.200     +0.2   (bootstrap 오차)

Tsallis, Poisson likelihood, 0.15 < pT < 5.66 GeV/c (truth T = 135.0 MeV, n = 6.800)
species        N        T (MeV)      pull        n            pull   -2lnL/ndf
pi+-      198980   135.4 ± 0.7      +0.5   6.852 ± 0.052    +1.0   11.7/22
K+-        28668   135.3 ± 1.7      +0.2   6.805 ± 0.156     0.0   28.9/22
p/pbar     12109   133.7 ± 2.6      -0.5   6.695 ± 0.257    -0.4   22.9/22
incl.     239757   147.6 ± 0.7     +19.3   7.496 ± 0.057   +12.1   25.6/22  <- 닫히지 않아야 정상

Z 선형, 컷 없는 truth 380개:  m = 91.142 ± 0.079 (-0.6σ),  Γ = 2.14 ± 0.17 (-2.1σ)
```

- **입자 종류별로 fit하면 닫힙니다.** pull이 이제 0 주변 양쪽으로 흩어집니다. bin 중심 χ²로
  했을 때는 1.3/0.8/0.4σ로 모두 양수였고, 그것이 §3(3)의 방법 편향이었습니다.
- NBD $k$: binomial thinning은 NBD의 $k$를 보존하므로, 받아들인 다중도의 모멘트
  $k=\langle N\rangle^2/(\mathrm{Var}N-\langle N\rangle)$가 생성기의 2.2를 되찾아야 합니다.
- Z 폭의 −2.1σ는 이 표본 자체의 요동입니다(1000회 toy에서 이 fit은 편향이 없습니다).
- 세 질량이 섞인 inclusive 스펙트럼을 파이온 질량 하나로 fit하면 $T$가 13 MeV(19σ) 높게
  나옵니다. fit 품질(−2lnλ/ndf = 25.6/22)은 멀쩡한데 파라미터는 틀린 전형적인 경우입니다.

> **fit 품질이 좋다는 것은 함수가 데이터를 잘 따라간다는 뜻일 뿐,
> 파라미터가 물리적 의미를 갖는다는 뜻이 아니다.**

이것이 실험이 inclusive Tsallis fit 하나가 아니라 **identified-particle spectra**를
따로 발표하는 이유입니다. §7의 fit-범위 변형에서도 같은 사실이 드러납니다.
inclusive fit은 범위에 따라 T가 3.4 MeV 움직이지만, 참 모형인 파이온 fit은 움직이지 않습니다.

## 5. 두 생태계 비교

| | ROOT (C++/PyROOT) | uproot + awkward |
|---|---|---|
| 파일 | 네이티브 | 순수 Python으로 읽기/쓰기 |
| 루프 | 명시적 `for` 또는 RDataFrame 선언형 | 배열 연산(vectorized), 루프 없음 |
| 히스토그램 | `TH1F` (binning·오차 내장, `Sumw2()` 필요) | `numpy.histogram` + 직접 오차, 또는 `hist` |
| fit | `TF1` + MINUIT | `scipy.optimize` 또는 `iminuit`(=같은 MINUIT2) |
| 강점 | 실험 프레임워크 표준, 대용량 I/O, `TBrowser` | 디버깅·플로팅·ML 연동이 훨씬 편함 |

**실무 결론**: 데이터 접근과 대규모 skim은 ROOT/RDataFrame, 그 뒤의 탐색적 분석과
그림은 uproot/awkward/matplotlib. 둘 다 할 줄 알아야 합니다.
이 저장소는 같은 분석을 두 버전으로 유지합니다:

- `analysis/m01_analysis_uproot.py` — **실행·검증 완료**
- `analysis/m01_analysis_root.C`, `analysis/m01_analysis_pyroot.py` — Claude가 쓰는 클라우드
  컨테이너에는 ROOT가 없어 **이 저장소 계열에서는 미실행**입니다. 2026-11-07 검토에서
  실행하면 바로 터질 버그를 고쳤습니다. 매크로의 히스토그램이 입력 파일에 붙어 있어
  `fin->Close()` 뒤 JSON을 쓸 때 해제된 메모리를 읽었습니다(오답노트 #19).
  WSL2에서 `bash scripts/run_all.sh`를 돌리면 `m01_compare_arms.py`가 판정합니다.
  비교 전에 세 arm이 같은 입력 파일과 같은 컷을 썼는지부터 확인하고, ROOT arm이
  하나라도 빠지면 실패합니다.

---

## 6. 재현

```bash
micromamba activate hic
bash scripts/run_all.sh                  # WSL2 버전 기록까지 갱신하려면 RECORD_ENV=1
python scripts/check_results_drift.py    # 커밋된 결과 표와 같은지 확인
```

전부 `config/analysis.yaml`의 값과 고정 시드(표본 20260827, 검증 20261107,
체계 오차 bootstrap 20261108)로 결정됩니다. GitHub Actions가 매 push마다 ROOT를 뺀
같은 과정을 돌리고 결과 표가 1e-6 안에서 같은지 확인합니다.

---

## 7. fit 방법 검증과 체계 오차

**방법 검증** (`analysis/m01_fit_validation.py`, `results/figures/m01_fit_validation.pdf`).
closure test는 한 표본에 대한 질문이라 1–2σ 편차가 요동인지 편향인지 가르지 못합니다.
생성기 자신의 모형에서 실제 표본과 같은 통계로 독립 표본을 여러 번 만들고, 각 fit
방법의 pull = (fit − 참값)/오차를 모읍니다. 편향이 없으면 평균 0, 오차가 맞으면 폭 1입니다.

| fit | 파라미터 | 평균 pull (± 0.06) | pull 폭 |
|---|---|---|---|
| Tsallis, bin 중심 χ² (이전) | T: π / K / p | +0.85 / +0.74 / +0.65 | 1.02 / 1.01 / 1.06 |
| Tsallis, Poisson LL (현재) | T: π / K / p | 0.00 / +0.01 / 0.00 | 1.01 / 0.99 / 0.97 |
| Z, Neyman χ² (이전) | Γ | −0.81 | 1.08 |
| Z, Poisson LL, 신호만 (현재) | Γ | −0.10 | 0.97 |
| Z, Poisson LL + 배경 ≥ 0 | Γ | −0.33 | 0.97 |

Neyman χ²는 질량 pull 폭도 1.12로, 오차를 12 % 과소평가합니다. Poisson LL(신호만)의
Γ pull −0.10은 물리 단위 편향(−0.005 ± 0.007 GeV)이 없는데도 생깁니다. 작은 표본에서
Γ와 그 오차가 함께 움직이기 때문입니다. muon acceptance는 질량에 대해 +0.23 %/GeV로
증가하고, 그 결과 컷은 fit 피크를 +0.010 ± 0.002 GeV 옮깁니다.

**체계 오차** (`analysis/m01_fit_systematics.py`, `results/figures/m01_tsallis_systematics.pdf`).
fit 범위의 아래쪽(0.3, 0.5 GeV/c)과 위쪽(3.0, 4.5, 8.0 GeV/c), bin 수(20, 80)를 하나씩
바꿉니다. 변형으로 생긴 차이 Δ는 두 fit이 데이터를 공유하므로 **Δ 자신의 통계적 퍼짐**과
비교해야 합니다. 이 퍼짐은 입자 bootstrap 50회로 잽니다. 2배를 넘는 변형만 체계 오차로
치고, 묶음별 최대값을 제곱합합니다. 통계 오차와 체계 오차는 끝까지 따로 둡니다.

| 스펙트럼 | T (MeV) | n |
|---|---|---|
| inclusive, $m_\pi$ 형태 (물리 결과) | 147.6 ± 0.7 (stat) ± 3.5 (syst) | 7.496 ± 0.057 ± 0.209 |
| π만, 참 모형 (대조군) | 135.4 ± 0.7 (stat) ± 0.2 (syst) | 6.852 ± 0.052 ± 0.023 |

대조군의 0.2 MeV는 이 절차의 거짓 양성 하한입니다. 변형 7개를 2σ 기준으로 보면
하나쯤은 우연히 넘습니다. 흔히 쓰는 Barlow 근사 $\sqrt{|\sigma_v^2-\sigma_{\rm nom}^2|}$는
같은 데이터를 다시 binning할 때 0에 가까워져, 대조군에서 잡음을 2.7σ "유의"로 판정했습니다
(bootstrap 판정은 0.8σ, 오답노트 #18).
