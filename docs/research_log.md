# 연구 노트

세션마다 Question / Physics / Tool / Procedure / Result / Interpretation /
Problems / Next step 형식으로 기록합니다.

---

## Session 1 — 2026-08-27 · Phase 0 + Module 1

### Question
collider/heavy-ion 계산 환경을 어디에 어떻게 세울 것인가, 그리고
"ROOT 파일 만들고 읽어서 첫 pT/η 히스토그램 그리기"가 실제로 되는가.

### Physics
아직 QCD 없음. 다루는 개념은 (a) 상대론적 운동학 — $p_{\rm T}$, $\eta$, $y$,
불변질량, (b) 다중도 요동(NBD), (c) $p_{\rm T}$ 스펙트럼의 열적/멱함수 이중 구조
(Tsallis/Hagedorn), (d) Breit-Wigner 공명.

### Tool
uproot 5.7.6 / awkward 2.13 / numpy 2.4.6 / scipy 1.17.1 / matplotlib 3.11.1.
FastJet 3.5.1.4와 pyhepmc 2.16.1(HepMC3)도 설치·검증만 완료(사용은 Module 3–4).
ROOT·PYTHIA는 이번 세션 환경에서 설치 불가 → WSL2로 이관.

### Procedure
1. 세 환경(Windows 호스트 / Cowork 로컬 VM / 클라우드 컨테이너) 조사 →
   `docs/environment_report.md`
2. 저장소 뼈대 + `config/analysis.yaml` + `src/hicpy/{style,kinematics,toy_generator}.py`
3. `scripts/m01_make_toy_events.py` → TTree 20,000 이벤트 (19 MB)
4. `analysis/m01_analysis_uproot.py` → 관측량 4종 + Tsallis/BW fit
5. `analysis/m01_closure_test.py` → 생성기 입력 복원 검증
6. `environment/{environment.yml,setup_wsl2.sh,verify_environment.py}` 작성,
   verify는 컨테이너에서 실행(17 PASS / 0 FAIL / 6 SKIP)

### Result
- kinematics 단위 검증 4종 통과 (round-trip 오차 ~1e-16, back-to-back 45.6 GeV
  뮤온 → m = 91.200 GeV)
- $\langle N_{\rm ch}\rangle = 11.99 \pm 0.06$ (|η|<0.8), RMS 8.78
- $dN_{\rm ch}/d\eta$ (|η| < 0.5 평균) = 7.50
- Tsallis(inclusive): $T = 148.2 \pm 0.7$ MeV, $n = 7.484 \pm 0.057$, χ²/ndf 0.98
- Tsallis(종류별): π 136.0±0.7, K 136.3±1.6, p 136.1±2.6 MeV → **입력 135 MeV 복원**
- Z 피크: $m = 91.06 \pm 0.09$ GeV/c² (PDG 91.1876), $\Gamma = 2.06 \pm 0.23$ GeV

### Interpretation
1. **Closure test가 핵심 결과.** 종류별로는 닫히고 inclusive로는 20σ 어긋난다는
   사실이, χ²/ndf가 좋아도 파라미터가 틀릴 수 있음을 정량적으로 보여준다.
   실험이 identified-particle spectra를 따로 내는 이유와 직결된다.
2. Z 질량이 1.3σ 낮고 폭이 좁은 것은 $p_{\rm T}^\mu>20$ GeV 컷의 acceptance
   sculpting. 보정 없는 "측정"이 어떻게 편향되는지의 최소 예제.
3. 다중도 RMS/평균 ≈ 0.73 ≫ Poisson 0.29 — centrality 정의가 왜 가능한지의 씨앗.

### Problems
- 클라우드 컨테이너에서 conda-forge·root.cern·hepdata·opendata.cern.ch가 모두
  차단 → ROOT/PYTHIA/실제 데이터 접근 불가. Module 2부터 WSL2 필수.
- uproot 5.7의 RNTuple 기본 동작, jagged counter branch 공유 실패,
  awkward에서 `np.clip`/`np.isin` 사용 불가, numpythia 빌드 실패,
  마운트 폴더에서 git lock 삭제 불가 → 전부 `docs/troubleshooting.md`에 기록
- ROOT C++ / PyROOT 버전은 작성만 하고 **실행하지 못함**. 검증 미완.

### Next step
1. WSL2에서 `bash environment/setup_wsl2.sh` 실행 → `verify_environment.py` 출력 공유
2. `analysis/m01_analysis_root.C` / `m01_analysis_pyroot.py` 실행, uproot 결과와
   히스토그램 대조 (같은 파일 → 같은 숫자가 나와야 함)
3. Module 2: PYTHIA 8 minimum-bias pp 13 TeV, toy 생성기를 대체하고
   `dN_{\rm ch}/d\eta`를 ALICE 측정과 비교

---

## Session 2 — 2026-08-27 · WSL2 환경 실체화 + 두 arm 정합성 사전 정비

### Question
WSL2 `hic` 환경이 실제로 동작하는가. 그리고 실행해 보지 못한 ROOT arm이
uproot arm과 같은 숫자를 낼 수 있는 상태인가.

### Tool
micromamba 2.9.0 / ROOT 6.40.02 / PYTHIA 8 / FastJet 3.5.1.4 / iminuit 2.32.0,
Python 3.11.16, WSL2 kernel 6.6.114.1, **16 논리코어**.

### Result
`verify_environment.py`: **23 PASS / 0 FAIL / 2 SKIP**.
- ROOT: TFile/TTree/TH1 100 entries 왕복, RDataFrame 합계 일치
- PYTHIA 8: pp 13 TeV 10 이벤트 → 최종상태 입자 1473개
- FastJet: anti-kT R=0.4 → 예상대로 2 jets
- SKIP은 `pyhepmc`(패키지 이름 누락, 추가함)와 `lhapdf`(Module 6 예정)

### Problems (실행 전 코드 검토에서 발견)
ROOT arm을 돌리기 전에 두 arm이 **어긋날 수밖에 없는** 불일치 두 건을 찾아 고침:
1. 불변 수율의 $1/p_{\rm T}$를 bin 중심 vs 입자별로 적용 (오답노트 #7)
2. "η=0에서의 dN/dη"가 bin edge 위치에 따라 다른 bin을 고름 (오답노트 #8)

둘 다 입자별 가중치 / binning-무관 정의로 통일. 부수적으로 파이온 closure fit의
χ²/ndf가 0.91 → 0.50으로 개선되어 이전 방식에 편향이 있었음이 확인됨.

### Next step
1. WSL2에서 `bash scripts/run_all.sh` → ROOT arm 실행 → `m01_compare_arms.py`가
   두 arm 일치를 자동 판정
2. Module 2: PYTHIA 8 카드 작성 및 minimum-bias 생성

---

## Session 3 — 2026-11-05 · 점검 후 검증 장치 보강 (클라우드 컨테이너, ROOT 없음)

### Question
10주 만에 다시 연 저장소가 여전히 재현되는가. Module 2(PYTHIA)로 넘어가기 전에
"분석 버그인지 생성기 차이인지"를 가를 장치가 충분한가.

### Tool
uproot 5.7.7 / awkward 2.14 / numpy 2.5.3 / scipy 1.18.1 / matplotlib 3.11.2 /
pytest. ROOT 없음 → ROOT arm은 코드만 수정, 실행은 다음 세션(WSL2).

### Procedure
1. 전체 재독 + uproot arm 처음부터 재실행 → 커밋된 `m01_fit_results.json`과 대조
2. 점검 결과를 `docs/next_steps.md`로 정리
3. ROOT 없이 가능한 보완 항목 반영:
   - `m01_compare_arms.py`: 세 arm(uproot / PyROOT / C++ macro)을 스칼라 →
     히스토그램 bin별(내용·오차) → Tsallis (T, n) σ 단위의 3단계로 비교
   - C++ 매크로가 같은 스키마의 JSON 요약을 쓰도록 수정 (이전에는 대조 불가)
   - 세 arm의 track 선택 통일 (`pt < pt_max` 누락, 오답노트 #11)
   - `m01_closure_test.py`: 생성기 설정에서 해석적으로 계산한 ⟨N_ch⟩와
     dN/dη|_{|η|<0.5} 기대값 추가 (`--infile` 옵션도)
   - `hicpy.style.save_figure`: 그림 메타데이터 제거 → 바이트 단위 재현 (#10)
   - `tests/` 22건 (kinematics 7, toy_generator 9, pipeline 6) + GitHub Actions
   - `run_all.sh`: `|| true` 제거(FAIL은 중단), pytest 단계 추가
   - 로드맵 번호를 README의 Module 번호로 통일

### Result
- 재실행 수치 = 커밋 수치 (Z fit 파라미터만 상대 ~1e-9 최적화기 잡음)
- 새 closure: ⟨N_ch⟩ 11.988 ± 0.062 vs 기대 12.059 (−1.1σ),
  dN/dη 7.502 ± 0.041 vs 7.537 (−0.9σ). Tsallis 종류별 pull 1.5/0.8/0.4σ 불변
- 그림 8개: 두 번 그려 `cmp` 동일. 테스트 22/22 통과, `run_all.sh` 전 구간 exit 0
- `pt_max` 통일 후에도 uproot arm 수치는 변화 없음 (50 < pT < 60 GeV/c 입자가
  |η|<0.8 안에 없었음) — ROOT arm에서는 달라졌을 것

### Interpretation
1. 두 다중도 closure가 모두 약간 음(−1σ)인 것은 생성기의 η 샘플러가 CDF를
   `cumsum`(사각형 적분)으로, 기대값 계산은 사다리꼴로 하기 때문일 가능성이
   있다 (상대 ~2.5e-4). 통계 오차(0.5 %)보다 작아 지금은 무시하지만, 이벤트 수를
   100배 늘리면 보일 것이다 — 그때 샘플러를 사다리꼴 CDF로 바꾼다.
2. `invariant_mass` 테스트에서 pT = 50, η = 5(p ≈ 3.7 TeV)인 파이온의 질량이
   8자리까지만 맞는다. $E^2 - p^2$의 상쇄 때문이며 오답노트 #6의 정량적 예다.
3. "스칼라 5개 일치"로는 arm 정합성을 주장할 수 없었다. 코드를 읽는 것만으로
   cut 하나(#11)가 빠진 것을 찾았고, 이는 스칼라 비교에서도 보였겠지만
   히스토그램의 어느 bin에서 어긋나는지는 bin별 비교가 있어야 보인다.

### Problems
- ROOT arm은 여전히 실행 못 함. 새로 쓴 TF1 fit, JSON 출력, `compare_arms`의
  fit tolerance(0.05σ)는 **실행 검증 전**이다.
- `run_all.sh`가 컨테이너에서 `environment/versions_installed.txt`를 덮어써
  되돌렸다 (next_steps C3).

### Module 2 착수 (같은 세션, 실행 없이 작성만)
- `generators/pythia_mb_13TeV.cmnd`: 물리에 영향을 주는 설정은 기본값이라도 전부
  명시 (Monash 2013, NNPDF2.3 LO, tau0Max = 10 mm, MPI/ISR/FSR/CR on). 이벤트
  클래스(INEL>0, NSD)는 카드가 아니라 분석에서 고르도록 설계.
- `scripts/m02_run_pythia.py`: Module 1과 같은 TTree 브랜치 + `particle_charge`,
  `process_code`, `n_ch_eta1`, `weight`. provenance에 카드 전문과 sigmaGen 기록.
- 열린 결정: Module 1 분석의 "charged hadron = π/K/p"는 toy의 우연. PYTHIA
  샘플에서는 charge ≠ 0으로 바꿔야 ALICE primary 정의와 맞음 (module02 §5).
- `tests/test_m02_card.py` 13건: 카드가 필수 설정을 빠짐없이, 중복 없이 담는지.

### 하전입자 정의 변경 (같은 세션) → float32 버그 발견
- 결정 (a): `charge` 브랜치(int8) 추가, 세 arm 모두 `particle_charge ≠ 0`.
  ⟨N_ch⟩ 11.988 → 11.999, ⟨pT⟩ 0.551 → 0.584 GeV/c (Z 뮤온 ~270개가 |η|<0.8에).
- 뮤온이 pT 20–50 GeV/c에 들어왔는데 Tsallis fit 범위가 여전히 6.1 GeV/c에서
  끝나는 게 이상해서 파 보니: float32 가중치 + numpy 누적합 → 꼬리 bin의 Σw²가
  정확히 0 → 오차 0인 점이 fit에서 **조용히 제외**되고 있었다 (오답노트 #12).
  커밋돼 있던 Module 1 결과도 이 영향을 받았다 (고 pT bin 내용이 1/128 단위로 양자화).
- 수정: 읽자마자 float64, 히스토그램 헬퍼에서 캐스팅, PyROOT 헬퍼 double 반환,
  fit 범위를 config에 명시(`fits.tsallis_pt_range = [0.15, 6.0]`).
- 수정 후: inclusive n 7.484 → 7.479, χ²/ndf 0.98 → 1.08; 종류별 closure 불변
  (pull 1.3/0.8/0.4σ). 다중도 closure는 하드론만 세도록 명시 (NBD 검증이 목적).

### `m02_compare_data.py` (같은 세션)
- 이벤트 클래스를 저장된 플래그가 아니라 입자(INEL>0)와 process code(NSD)에서
  다시 계산한다; 저장된 `n_ch_eta1`과 다르면 경고 (생성 스크립트와 분석의
  "charged"/"eta" 정의 불일치 검출).
- 발표 dN/dη는 pT → 0 외삽값이므로 비교에는 pT cut을 걸지 않는다 — Module 1의
  `pT > 0.15` 히스토그램을 재사용하지 않는 이유.
- 테스트: MC 자신을 HEPData YAML + provenance JSON으로 둔갑시켜 넣으면 ratio가
  정확히 1이어야 한다(6건 통과). 실데이터는 provenance를 채운 뒤에만.

### Next step
`docs/next_steps.md` Step 1: WSL2에서 `run_all.sh` → `ARMS AGREE` 확인 →
`results/` 재커밋. 이어서 Step 2: `m02_run_pythia.py --n-events 10`부터,
그 다음 HEPData 표 provenance 채우기 → `m02_compare_data.py`.

---

## Session 4 — 2026-11-19 · ALICE dN/dη 표 확보

- HEPData가 자동 요청(urllib·curl 모두)을 Cloudflare 403으로 막아, 사용자가 브라우저로
  `HEPData-ins1395253-v1-Table_1.yaml`을 받음 → `data/raw/alice_pp13TeV_dndeta_table1.yaml`
  (payload는 gitignore) → `m02_fetch_hepdata.py --local`로 sha256 기록·`verified: true`.
- 표 내용: 열 2개(INEL, INEL>0), |η| < 1.8에서 0.2 폭 18 bins. INEL>0 열(index 1) 자동 선택.
  오차는 점마다 **라벨 없는 대칭 오차 하나** → 비교 스크립트는 계통 오차로 취급(stat = 0).
- INEL>0, |η| < 0.4 (중앙 4 bins) 평균: **6.44 ± 0.19**.
- 발견: 판정용 "|η| < 0.5" 중앙값이 데이터 쪽에서는 bin 경계 때문에 실제로 |η| < 0.4였고,
  MC는 |η| < 0.5로 평균 → 범위 불일치. MC도 데이터 bin이 덮는 범위로 맞추도록 수정
  (`verdict.eta_range`에 기록).
- toy 파일로 체인 동작 확인: MC/data = 1.18. **물리적 의미 없음** (toy는 데이터에 맞춘 적이 없음).

### Next step
WSL2에서 PYTHIA 8 생성 → `python analysis/m02_compare_data.py` → 판정 (±15 %).

### PYTHIA 첫 실행 (같은 날)
- WSL2 (PYTHIA 8.312), 10 이벤트: 정상 종료. sigma_gen = 78.05 mb (ND 56.4, SD 6.4 × 2,
  DD 8.8 mb). INEL>0 5/10 — 10개 중 6개가 회절이라 통계 요동.
- 첫 이벤트가 stdout에 출력됨 → 카드 줄 순서 문제 발견·수정 (오답노트 #14).
- 클라우드 컨테이너에서 PyPI `pythia8mc` (**8.318**, WSL2와 버전 다름)로 2,000 이벤트
  예비 실행: INEL>0 85.5 %, process mix ND 72 / SD 17 / DD 11 %.
  ALICE INEL>0 비교 (|η| < 0.4): MC 7.08 ± 0.17 vs data 6.44 ± 0.19 → **MC/data = 1.10**
  (±15 % 기준 안). η 전 구간에서 MC가 1–13 % 높음.
  **예비값**: 버전·통계가 정식 조건이 아님. 결과는 커밋하지 않음.
- 10 % 높은 원인 후보 (확인 전): primary 정의 차이 (tau0Max = 10 mm에서 Σ±, Ξ±, Ω± 등
  cτ > 1 cm 하전 입자가 그대로 남아 세어짐 — ALICE 정의와 일치하는지 논문 확인 필요),
  Monash 튠 자체의 과대 예측, 통계.

### Next step
WSL2에서 `git pull` 후 20,000 이벤트 → `python analysis/m02_compare_data.py` →
정식 결과 커밋. 이어서 primary 정의를 논문과 대조.

### Module 1 종료 — 세 arm 일치 (같은 날, WSL2, ROOT 6.40.02)
- 처음으로 PyROOT, C++ 매크로를 실행해 uproot arm과 3단계 대조 (`m01_compare_arms.py`).
- 첫 실행에서 나온 문제 세 개를 고침:
  #15 MINUIT의 C 범위(10³⁰)가 정규화를 얼림 → PyROOT fit T = 56 MeV, χ²/ndf ≈ 10⁴;
  #16 히스토그램이 입력 TFile 소유 → `fin->Close()` 후 segfault, 요약 JSON 0 바이트;
  #17 C++에서 `Float_t` 곱 → 전방 η가 float32 정밀도, dN/dη 19 bins 어긋남.
- 최종: 세 arm이 스칼라 5개, 히스토그램 220 bins(dN/dη 한 bin만 경계 입자 1개), Tsallis
  (T, n) 모두 일치 → **ARMS AGREE**. 결과물 WSL2에서 재생성해 커밋 (`eabd61f`).

### Module 2 첫 결과 — PYTHIA 8 vs ALICE (같은 날, WSL2, PYTHIA 8.312)
- 20,000 이벤트, seed 20260827, 35 s. σ_gen = 78.05 mb; INEL>0 86.4 %;
  ND / SD / DD = 72.7 / 16.3 / 11.0 %. PYTHIA 내부 오류 0.2 % 수준, 모두 재시도로 처리.
- **ALICE INEL>0 dN/dη, |η| < 0.4: MC 6.84 ± 0.05 vs data 6.44 ± 0.19 → MC/data = 1.062**.
  η 모양은 재현, 크기는 |η| < 1.8 전 구간에서 3–7 % 높음 (중앙에서 가장 큼).
  ±15 % 기준 통과, 차이는 약 2σ (데이터 오차가 지배).
- 해석(가설): 약붕괴 strange baryon(Σ±, Ξ, Ω)은 |η| < 0.4 하전 트랙의 ~1 % (2k 예비 샘플)
  → primary 정의 차이로 설명되는 몫은 작음. Monash 튠 자체의 과대 예측이 유력.
- 결과: `results/tables/m02/m02_compare_summary.json`, `results/figures/m02/`.

### Next step
1. ALICE 논문의 primary 정의 확인 → tau0Max 설정과 대조.
2. 튠·colour reconnection을 바꿔 dN/dη 민감도 확인.
3. 같은 ALICE 측정의 pT 스펙트럼·다중도 분포로 관측량 확대.
4. Module 1 분석 코드(`m01_analysis_uproot.py`)를 PYTHIA 파일에 수정 없이 적용.

---

## Session 5 — 2026-11-19 · 6 % 원인, 관측량 확대, Module 3

### 6 % 원인
- 문헌: ALICE primary 정의 = strange 약붕괴 산물만 제외 → 카드 tau0Max = 10 mm와 일치.
- 민감도 (클라우드, 8.318, 각 20k): CR off +37 %, 4C +7 %, pT0Ref +5 % → −3.4 %, 모든 붕괴 +9.5 %.
- INEL 열 대조: 데이터 f(INEL>0) = 0.821, PYTHIA 0.860 → INEL에서만 +5 % 추가. INEL>0의 6 %와는 별개.
- 결론: 정의·이벤트 클래스는 아님. MPI 튠(13 TeV 외삽)이 가장 유력. 상세: module02_pythia.md §7.

### 관측량
- `analysis/m02_compare_pt.py`: ALICE Table 2 (pT 불변 수율, |η| < 0.8, INEL>0) 비교 스크립트.
  저 pT / 고 pT 비를 따로 보고 → 6 %가 soft 성분인지 판별용. 표는 HEPData 수동 다운로드 대기.
- Module 1 분석이 PYTHIA 파일에 수정 없이 동작 확인. 빈 Z fit 버그(#19), 범례 문구 수정.

### Module 3 (HepMC)
- `src/hicpy/eventrecord.py`: PYTHIA → HepMC3 변환기 (PyPI PYTHIA에 Pythia8ToHepMC3 없음).
  두 번 실패 후 완성 (#18): 엄마 집합 묶기 ✗, PYTHIA 번호순 쓰기 ✗ (ISR backward evolution으로
  엄마가 뒤 번호), 위상 정렬 ✓.
- 100 이벤트: 전역 4-운동량 1.5e-8, 전하 100/100, 개수 100/100. 꼭짓점별 보존은 붕괴 100 %,
  강입자화 99.7 %, FSR 1 %, remnant 0.3 % → 원래 검증 기준이 PYTHIA에 맞지 않았음을 확인,
  로드맵 기준 수정.
- 테스트 51개 (Module 3 5개, pT 비교 2개 추가). CI에 pyhepmc 추가.

### Next step
1. ALICE Table 2 다운로드 → `m02_compare_pt.py` → 6 %가 soft 성분인지 확인.
2. WSL2에서 Module 3 재현, `scripts/m02_run_variations.sh`로 민감도 재현 (8.312).
3. Module 4 FastJet.

### pT 스펙트럼 (같은 날)
- ALICE Table 2 (사용자 브라우저 다운로드 → `--local` 검증): |η| < 0.8, INEL>0, 0.15–20 GeV/c, stat/sys 분리.
- Monash (8.318, 20k): MC/data 0.86 (0.15–0.2) → 1.24 (0.5–0.55) → ~1.08 (1–2) → ~1.2 (2–8).
  적분 1.09 (dN/dη 결과 1.076과 일치). ⟨pT⟩ 0.667 vs 0.647 (+3 %).
- 변형: pT0Ref는 크기만 (⟨pT⟩ −0.8 %), CR off는 모양 (⟨pT⟩ −14 %). 4C도 모양 개선 없음.
- **결론: 6 % dN/dη 초과는 0.3–1 GeV/c 벌크 영역의 초과 + 스펙트럼 모양 차이. 손잡이 하나(pT0Ref)로는
  dN/dη만 맞고 pT 모양은 못 맞춘다 → MPI와 CR을 함께 다시 맞춰야 함.**

---

## Session 6 — 2026-11-19 · 2차원 스캔, Module 4 (FastJet)

### 2차원 스캔 (pT0Ref × CR range, 20점 × 10k, 클라우드 8.318)
- 계통 상관: 논문에서 확인되는 것은 INEL>0 정규화 1.2 %뿐 → uncorr / corr 두 극단으로 χ².
- dN/dη 최적: 두 모형 모두 pT0Ref 2.52, CR range 1.4 (MC/data 0.999). Monash는 χ² 61/18.
- pT 모양: 어디서도 안 맞음 (corr에서 χ² ≥ 738/47). → 다음 손잡이는 StringPT:sigma 등.

### Module 4
- 카드: HardQCD, pT̂ > 20 GeV, Monash; σ = 0.722 mb; 5k 이벤트.
- 면적 (R = 0.4): anti-kT 0.507 (96 % within 10 % of πR²), C/A 0.565, kT 1.027.
- embedding (ρ_true = 15): δpT 7.24 ± 0.17 vs ρA 7.57; median ρ = 18.15 (= 15 + PYTHIA UE 2.4 + 4 %);
  차감 후 −0.71 ± 0.17, 68 % 폭 2.57 vs Poisson 2.75.
- 버그/함정: 첫 비교에서 신호 제트의 UE 차감을 빠뜨려 −2.7 GeV; FastJet ClusterSequence 수명,
  고스트 seed (#20).
- 테스트 56개. CI에 fastjet.
- WSL2 재현 스크립트: `scripts/run_all_pythia.sh`.

### Next step
1. WSL2: `bash scripts/run_all_pythia.sh` → git diff로 8.318 vs 8.312 비교 후 커밋.
2. pT 모양: StringPT:sigma 1차원 스캔.
3. 로드맵 Module 5 (LHAPDF, nPDF).

### WSL2 정식 재생성 (같은 날, `1cc643e`)
- `run_all_pythia.sh` 첫 실행이 Table 2 원본 부재로 중단 + 스캔 결과 삭제 → preflight·단계 선택 추가 (#21).
- 스캔 χ²에서 빈 고 pT MC bin이 χ²를 ~100씩 부풀림 → pT χ²를 < 5 GeV/c로 제한.
- 정식 결론: **pT0Ref 2.28 → 2.40 (CR 그대로)이면 dN/dη MC/data 1.003–1.004** (변형·스캔 독립 일치).
  pT 모양은 미해결. Module 3·4 수치는 클라우드와 일치 (module03/04 문서 "정식 결과").
- pyhepmc를 WSL2 hic 환경에 설치 (`micromamba install -n hic -c conda-forge pyhepmc`).

---

## Session 7 — 2026-11-19 · StringPT 스캔, Module 5 (PDF)

### StringPT:sigma (pT0Ref 2.40, 클라우드 8.318, 6 × 20k)
- σ 0.335 → 0.31: ⟨pT⟩ MC/data 1.022 → 0.995, 최저 pT 0.88 → 0.92, χ² 152 → 126 / 39.
- 0.4–0.8 GeV/c의 +12–14 %는 σ·pT0Ref·CR 어느 것으로도 안 움직임 → identified 스펙트럼 필요.

### Module 5
- LHAPDF 불가 (PyPI 없음, 세트 서버 차단) → PYTHIA 내장 격자 (NNPDF2.3 LO, EPS09 LO Pb).
- 합 규칙: u_v 2.003, d_v 1.001, s_v 0.000, 운동량 1.000 (Q² = 10), 0.997 (10⁴). 글루온 몫 0.47 → 0.50.
- LHC x: 저 pT 하드론 1.5e-4 (pp 13) / 4e-4 (Pb-Pb 5.02); Z 7e-3 / 1.8e-2; 100 GeV 제트 1.5e-2 / 4e-2.
- EPS09 R_g^Pb: 하드론 영역 0.67 ± 0.28 (90 % CL), Z 1.06 ± 0.03, 제트 1.09 ± 0.06.
  → 저 pT 억제의 일부는 냉핵물질(PDF)만으로도 가능, 불확실성 큼; 제트 억제는 PDF로 설명 불가.
- 테스트 5개 추가 (합 규칙, 동위스핀, shadowing 부호). CI에 pythia8mc.
