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

### Next step
`docs/next_steps.md` Step 1: WSL2에서 `run_all.sh` → `ARMS AGREE` 확인 →
`results/` 재커밋. 이어서 Step 2: `m02_run_pythia.py --n-events 10`부터.
