# Module 1 — ROOT 파일과 첫 물리 히스토그램

실행일: 2026-08-27 · 실행 환경: Track B(클라우드 컨테이너, uproot arm)
· ROOT arm은 **미실행** (Track A 대기)

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
| Z→μ⁺μ⁻ | 2 % 이벤트에 relativistic Breit-Wigner ($m_Z=91.1876$, $\Gamma_Z=2.4952$ GeV, PDG) 후 등방 2체 붕괴 + boost | 불변질량 재구성 연습용 |

**이 생성기가 아닌 것**: hard scattering도, parton shower도, hadronization도,
입자간 상관도 없습니다. 이벤트는 통계적으로만 그럴듯합니다. Module 2에서 PYTHIA로
교체하면 **분석 코드는 그대로 두고** 입력만 바뀌어야 합니다 — 그게 이 구조의 목적입니다.

전체 설정은 `config/analysis.yaml`에, 실제 사용된 값·시드·git 커밋은
`data/generated/m01_toy_events.json`과 ROOT 파일 내부 `provenance`에 기록됩니다.

---

## 3. 관측량과 결과

이벤트 20,000개, 시드 20260827, |η| < 0.8, pT > 0.15 GeV/c.

### (1) 하전입자 다중도 분포 `results/figures/m01_multiplicity.pdf`

$\langle N_{\rm ch}\rangle = 11.99 \pm 0.06$ (RMS 8.78, |η|<0.8).
RMS/평균 ≈ 0.73 으로 Poisson(≈0.29)보다 훨씬 넓습니다 — NBD의 $k$가 만든
**이벤트별 요동**이고, 중이온 물리에서 centrality를 정의할 때 바로 이 폭을 씁니다.

### (2) $dN_{\rm ch}/d\eta$ `results/figures/m01_dndeta.pdf`

η=0에서 7.54. 실제 ALICE pp √s=13 TeV INEL>0 측정은 6.5 근처이므로 toy는 조금
높습니다 — 의도적으로 맞추지 않았습니다. Module 2에서 PYTHIA로 실제 데이터와
비교하는 것이 목표입니다.

### (3) 불변 $p_{\rm T}$ 스펙트럼 `results/figures/m01_pt_spectrum.pdf`

$$\frac{1}{2\pi p_{\rm T} N_{\rm ev}}\frac{d^2N}{dp_{\rm T}\,d\eta}$$

로 그렸습니다. $2\pi p_{\rm T}$로 나누는 이유: $d^2N/(dp_x dp_y) = d^2N/(2\pi p_T dp_T)$
이므로 이 조합이 **Lorentz 불변 수율**에 해당하고, 실험 논문이 이 형태로 발표합니다.

Tsallis fit 결과 (하전입자 전체, 파이온 질량 가정):

```
T = 147.9 ± 0.6 MeV,   n = 7.475 ± 0.055,   chi2/ndf = 1.09
```

입력값 $T=135$ MeV, $n=6.8$ 과 **맞지 않습니다.** 이건 버그가 아니라 물리입니다
(§4 참조).

### (4) 이중뮤온 불변질량 `results/figures/m01_dimuon_mass.pdf`

$p_{\rm T}^\mu > 20$ GeV/c, $|\eta_\mu| < 2.4$, 정확히 2개 뮤온. 241 후보.

```
m_fit = 91.06 ± 0.09 GeV/c²      (PDG: 91.1876)
Γ_fit = 2.06  ± 0.23 GeV         (PDG: 2.4952)
chi2/ndf = 0.43
```

중심값은 PDG보다 1.3σ 낮고 폭은 좁게 나옵니다. 원인은 **acceptance sculpting**:
$p_{\rm T}^\mu>20$ GeV 컷이 Breit-Wigner의 낮은 질량 꼬리를 비대칭적으로 잘라냅니다.
실험이 Z 질량을 측정할 때 acceptance·efficiency 보정을 반드시 하는 이유입니다.

---

## 4. Closure test — 분석이 정답을 되찾는가

`analysis/m01_closure_test.py`. 생성기의 진짜 파라미터를 알고 있으므로, 분석이
그것을 복원해야 합니다.

```
generator truth:  T = 135.0 MeV,  n = 6.800

species          N          T (MeV)   pull               n   pull  chi2/ndf
pi+-        198980    136.1 +-  0.6    1.6    6.873 +- 0.049    1.5      0.91
K+-          28668    136.0 +-  1.6    0.6    6.829 +- 0.149    0.2      1.18
p/pbar       12109    135.7 +-  2.6    0.3    6.910 +- 0.279    0.4      0.91
incl.       239757    147.9 +-  0.6   20.1    7.475 +- 0.055   12.2      1.09
```

**입자 종류별로 fit하면 닫힙니다** (pull < 2σ). 반면 세 질량이 섞인 inclusive
스펙트럼을 파이온 질량 하나로 fit하면 $T$가 13 MeV(20σ) 높게 나옵니다.
χ²/ndf = 1.09 로 "fit은 잘 됐는데" 파라미터는 틀린 전형적인 경우입니다.

> **χ²/ndf가 1에 가깝다는 것은 함수가 데이터를 잘 따라간다는 뜻일 뿐,
> 파라미터가 물리적 의미를 갖는다는 뜻이 아니다.**

이것이 실험이 inclusive Tsallis fit 하나가 아니라 **identified-particle spectra**를
따로 발표하는 이유입니다.

---

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
- `analysis/m01_analysis_root.C`, `analysis/m01_analysis_pyroot.py` — **미실행**
  (Track B에 ROOT 없음). Track A에서 실행 후 결과를 `docs/research_log.md`에 기록할 것.

---

## 6. 재현

```bash
micromamba activate hic
bash scripts/run_all.sh
```

전부 `config/analysis.yaml`의 값과 고정 시드(20260827)로 결정됩니다.
