# heavy-ion-computing

중이온 충돌(Heavy-Ion Collision) 및 collider 현상론을 위한 **계산 워크플로 구축 프로젝트**.

목표 체인:

```
이론/모델 → 이벤트 생성 → 입자 수준 이벤트 → 실험 데이터 포맷
        → ROOT 분석 → 물리 관측량 → 실제 실험 데이터와 비교
```

이 저장소는 튜토리얼 모음이 아니라, **재현 가능한 연구 환경**을 단계적으로 쌓아
올리는 것을 목적으로 합니다. 6개월 뒤에 다시 열어도 `bash scripts/run_all.sh`
한 줄로 모든 그림과 숫자가 재생성되어야 합니다.

---

## 1. 현재 상태 (2026-11-09)

| 모듈 | 내용 | 상태 |
|---|---|---|
| Module 0 | HEP 계산 환경 구축 (WSL2 + micromamba) | **완료** — WSL2 `hic` 환경, verify 23 PASS / 0 FAIL / 2 SKIP |
| Module 1 | ROOT 파일 I/O, 첫 물리 히스토그램 | uproot arm 완료: Poisson-likelihood fit, closure(생성기 성분 4종), pseudo-experiment 방법 검증, 체계 오차, 테스트 65개, CI 실가동 / **ROOT arm은 이 계열에서 실행 대기** |
| Module 2 | PYTHIA 8 minimum bias pp 13 TeV | 카드·생성 스크립트·설정 작성 완료, **WSL2 실행 대기** |
| Module 3 | HepMC · PDG 이벤트 레코드 | 예정 |
| Module 4 | FastJet | 예정 |

자세한 로드맵: [`docs/roadmap.md`](docs/roadmap.md) ·
보완 사항과 바로 다음 할 일: [`docs/next_steps.md`](docs/next_steps.md)

> **주의 (2026-11-09).** 사용자 로컬 저장소에는 아직 push되지 않은 Session 3·4 작업이
> 있습니다(Notion 작업 로그: ROOT arm 실행 `d3379f8`, 2026-09-26 독립 검증).
> GitHub의 Claude 브랜치들은 그 작업을 모른 채 `main`에서 갈라져 나왔습니다.
> 합치기 전에 [`docs/next_steps.md`](docs/next_steps.md) §3을 보세요.

---

## 2. 환경 (Module 0)

### 왜 WSL2 + micromamba 인가

Windows 네이티브에서 ROOT/PYTHIA/FastJet을 쓰는 것은 가능하지만 사실상 모든
HEP 소프트웨어의 빌드 시스템·문서·커뮤니티가 Linux를 전제로 합니다. 세 가지
후보 중:

- **WSL2 + micromamba (채택)** — Windows 파일에 그대로 접근 가능, GUI(ROOT
  TBrowser)도 WSLg로 동작, conda-forge 바이너리라 컴파일 불필요, 호스트 CPU
  8코어를 전부 사용 가능.
- Docker — 재현성은 최고지만 GUI·파일 공유·디스크 사용량 부담.
- Linux VM — 리소스 격리가 심하고 파일 왕래가 번거로움.

패키지 관리자는 `conda`가 아니라 **micromamba**를 씁니다. 단일 정적 바이너리이고
base 환경을 만들지 않으며 solver가 훨씬 빠릅니다.

### 설치

```bash
# WSL2 Ubuntu 셸에서
git clone <this-repo> && cd heavy-ion-computing
bash environment/setup_wsl2.sh          # 시스템 라이브러리 → micromamba → hic 환경 → 검증
micromamba activate hic
```

설치되는 것: `environment/environment.yml` 참조 (ROOT, PYTHIA 8, HepMC3,
FastJet, uproot/awkward/hist/vector, iminuit, JupyterLab).
**MadGraph, Delphes, LHAPDF, Rivet, JETSCAPE 등은 아직 설치하지 않습니다** —
해당 모듈에 도달할 때 추가합니다.

### 검증 (설치 성공 ≠ 동작)

```bash
python environment/verify_environment.py --write environment/versions_installed.txt
```

각 패키지를 import만 하는 것이 아니라 실제로 **사용**합니다: ROOT 파일 생성·재읽기,
PYTHIA 10 이벤트 생성 후 최종상태 입자 수 확인, FastJet anti-kT 클러스터링,
LHAPDF PDF 값 평가. 결과는 `environment/versions_installed.txt`에 기록되며
이 파일이 재현성의 근거입니다.

---

## 3. 재현 방법

```bash
micromamba activate hic
bash scripts/run_all.sh                  # WSL2에서 환경 기록까지 갱신: RECORD_ENV=1 bash scripts/run_all.sh
python scripts/check_results_drift.py    # 커밋된 결과 표가 그대로 재현됐는지
```

수행 순서 (전체 출력은 `results/logs/run_all_<시각>.log`에도 남습니다):

1. `environment/verify_environment.py` — 환경 검증 (FAIL이면 중단). 버전 기록
   `environment/versions_installed.txt`는 `RECORD_ENV=1`일 때만 다시 씁니다
2. `pytest` — 단위·파이프라인·교차검증 판정 테스트 (`tests/`, 65건; `RUN_TESTS=0`이면 생략)
3. `scripts/m01_make_toy_events.py` — toy 이벤트 20,000개 생성 → ROOT TTree
4. `analysis/m01_analysis_uproot.py` — 관측량 4종 + Poisson-likelihood fit → `results/`
5. `analysis/m01_closure_test.py` — **closure test**: 다중도 평균·NBD k, dN/dη, 종류별
   Tsallis, truth 수준 Z 선형을 파일 provenance의 생성기 입력과 대조
6. `analysis/m01_fit_validation.py` — pseudo-experiment로 fit **방법**의 편향과 오차 포함률 측정
7. `analysis/m01_fit_systematics.py` — fit 범위·binning 변형으로 Tsallis 체계 오차
   (통계 오차와 분리, 참 모형 대조군 포함)
8. ROOT가 있으면 `m01_analysis_pyroot.py`와 `m01_analysis_root.C`(컷 값은 YAML에서 전달)를
   실행하고, `analysis/m01_compare_arms.py`가 같은 입력·같은 컷인지 확인한 뒤 세 arm을
   스칼라·bin별·fit으로 대조 (ROOT arm이 하나라도 빠지면 실패)

GitHub Actions(`.github/workflows/ci.yml`)가 push마다 ROOT를 뺀 전 과정을 돌리고,
커밋된 `results/tables/*.json`이 재현되는지 `scripts/check_results_drift.py`로 확인합니다
(중심값 1e-6, 기계마다 흔들리는 수치 Hessian 오차는 1e-2). 그림은 메타데이터를
제거해 저장하므로 같은 기계·같은 seed면 **바이트까지 같습니다**.

모든 물리 파라미터는 코드가 아니라 [`config/analysis.yaml`](config/analysis.yaml)에
있습니다. 이 파일이 결과의 provenance입니다.

---

## 4. 디렉토리 구조

```
heavy-ion-computing/
├── config/analysis.yaml        # 모든 cut·binning·생성기 파라미터 (숫자는 여기에만)
├── environment/                # environment.yml, setup_wsl2.sh, verify_environment.py
├── src/hicpy/                  # 재사용 라이브러리: kinematics, events, histograms, fitting,
│                               #   toy_generator, pdg(기준값), provenance, style
├── scripts/                    # 이벤트 생성 / run_all.sh / check_results_drift.py
├── analysis/                   # 관측량 계산 (uproot arm, ROOT arm, closure, 방법 검증, 체계 오차)
├── generators/                 # PYTHIA 카드 (pythia_mb_13TeV.cmnd)
├── data/{raw,generated,processed}/
├── results/{figures,tables,logs}/
├── docs/                       # README 외 모든 문서
│   ├── environment_report.md   # Phase 0 환경 조사 결과
│   ├── roadmap.md              # 마일스톤
│   ├── research_log.md         # 연구 노트 (세션별)
│   ├── troubleshooting.md      # 오답노트 (오류·실패 기록)
│   └── modules/module0N_*.md   # 모듈별 코드·결과 설명
├── exercises/                  # 모듈별 연습문제
├── tests/                      # pytest: kinematics, 생성기, fit, 교차검증 판정, 소규모 end-to-end
└── notebooks/
```

**규칙**: 노트북은 오케스트레이션만, 계산은 `src/`. 두 번 쓰이거나 15줄이 넘는
함수는 `src/hicpy/`로 옮깁니다.

---

## 5. 데이터 출처

현재 커밋에는 **실험 데이터가 없습니다**. Module 1의 이벤트는 toy Monte Carlo이며
`data/generated/*.json`에 생성기 설정·시드·git 커밋이 함께 기록됩니다.

실제 데이터는 Module 9에서 HEPData / CERN Open Data / ALICE·CMS·ATLAS 공개
데이터에서 받아오며, 받을 때마다 DOI·실험·충돌계·√s_NN·관측량·선택 조건·불확실성을
`data/raw/<dataset>.json`에 기록합니다. 데이터를 지어내지 않습니다.

---

## 6. GitHub / Notion

- 원격: `https://github.com/NhWee/heavy-ion-computing.git` — 설정과 첫 push 절차는
  [`docs/github_setup.md`](docs/github_setup.md)
- Notion 허브: [HIC Computing — 연구 허브](https://app.notion.com/p/3c96974402a681d5bbfdd050e0af6ebe)
  - 작업 로그 (Work Log) — 세션별 기록, 커밋 해시·seed 포함
  - 실패 보고 (Failure Report) — 오답노트, 심각도·영역·가드 위치 포함

**"무엇을 했다"는 GitHub, "왜 그렇게 했고 무엇이 안 됐다"는 Notion.**

## 7. 문서 규약

- [`docs/research_log.md`](docs/research_log.md) — 세션마다 Question / Physics /
  Tool / Procedure / Result / Interpretation / Problems / Next step
- [`docs/troubleshooting.md`](docs/troubleshooting.md) — **오답노트**. 실제로 만난
  오류, 원인, 해결, 그리고 "왜 그런 오류가 나는지"의 물리적/소프트웨어적 이유
- `docs/modules/moduleNN_*.md` — 각 모듈의 코드와 결과 설명
- `docs/howto/*.md` — 따라 하기 가이드 (명령 → 무엇을 왜 하는지 → 화면에 보일 것).
  첫 번째: [`fetch_hepdata.md`](docs/howto/fetch_hepdata.md) 실험 데이터 받아 MC와 비교하기
