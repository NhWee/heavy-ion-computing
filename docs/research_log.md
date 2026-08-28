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
