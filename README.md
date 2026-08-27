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

## 1. 현재 상태 (2026-08-27)

| 모듈 | 내용 | 상태 |
|---|---|---|
| Module 0 | HEP 계산 환경 구축 (WSL2 + micromamba) | 스크립트 작성 완료, **사용자 실행 대기** |
| Module 1 | ROOT 파일 I/O, 첫 물리 히스토그램 | **완료** (uproot arm 실행·검증 완료 / ROOT arm 미실행) |
| Module 2 | PYTHIA 8 | 예정 |
| Module 3 | HepMC · PDG 이벤트 레코드 | 예정 |
| Module 4 | FastJet | 예정 |

자세한 로드맵: [`docs/roadmap.md`](docs/roadmap.md)

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
bash scripts/run_all.sh
```

수행 순서:

1. `environment/verify_environment.py` — 환경 검증 및 버전 기록
2. `scripts/m01_make_toy_events.py` — toy 이벤트 20,000개 생성 → ROOT TTree
3. `analysis/m01_analysis_uproot.py` — 관측량 4종 + fit → `results/`
4. `analysis/m01_closure_test.py` — **closure test** (분석이 생성기 입력값을 되찾는지)
5. ROOT가 설치되어 있으면 `m01_analysis_pyroot.py`와 `m01_analysis_root.C`도 실행

모든 물리 파라미터는 코드가 아니라 [`config/analysis.yaml`](config/analysis.yaml)에
있습니다. 이 파일이 결과의 provenance입니다.

---

## 4. 디렉토리 구조

```
heavy-ion-computing/
├── config/analysis.yaml        # 모든 cut·binning·생성기 파라미터 (숫자는 여기에만)
├── environment/                # environment.yml, setup_wsl2.sh, verify_environment.py
├── src/hicpy/                  # 재사용 라이브러리: style, kinematics, toy_generator
├── scripts/                    # 이벤트 생성 / run_all.sh
├── analysis/                   # 관측량 계산 (uproot arm, ROOT arm, closure test)
├── generators/                 # (예정) PYTHIA·Angantyr 등 생성기 설정 카드
├── data/{raw,generated,processed}/
├── results/{figures,tables,logs}/
├── docs/                       # README 외 모든 문서
│   ├── environment_report.md   # Phase 0 환경 조사 결과
│   ├── roadmap.md              # 마일스톤
│   ├── research_log.md         # 연구 노트 (세션별)
│   ├── troubleshooting.md      # 오답노트 (오류·실패 기록)
│   └── modules/module01_root.md
├── exercises/                  # 모듈별 연습문제
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

## 6. 문서 규약

- [`docs/research_log.md`](docs/research_log.md) — 세션마다 Question / Physics /
  Tool / Procedure / Result / Interpretation / Problems / Next step
- [`docs/troubleshooting.md`](docs/troubleshooting.md) — **오답노트**. 실제로 만난
  오류, 원인, 해결, 그리고 "왜 그런 오류가 나는지"의 물리적/소프트웨어적 이유
- `docs/modules/moduleNN_*.md` — 각 모듈의 코드와 결과 설명
