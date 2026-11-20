# 보완 사항과 Next step — 2026-11-05 점검 (같은 날 1차 반영)

마지막 작업(Session 2, 2026-08-27) 이후 약 10주 만의 점검.
저장소 전체(코드·config·결과·문서)를 읽고, 클라우드 컨테이너에서 uproot arm을
처음부터 다시 실행해 확인한 뒤, ROOT가 필요 없는 항목은 **Session 3에서 바로
반영**했다 (`docs/research_log.md` Session 3).

---

## 0. 점검 결과 요약

| 항목 | 점검 시 | 반영 후 |
|---|---|---|
| uproot arm 재현성 | ✅ 같은 seed → 커밋된 수치와 일치 | ✅ CI가 매 push마다 커밋된 JSON과 대조 |
| Closure test | ✅ Tsallis (T, n)만 | ✅ 다중도·dN/dη 해석적 기대값 추가 |
| ROOT arm (PyROOT, C++ macro) | ❌ 한 번도 실행 안 됨 | ❌ **여전히 WSL2 실행 대기** (A1) |
| `m01_compare_arms.py` | ⚠️ 스칼라 5개, PyROOT만 | ✅ 세 arm, 히스토그램 bin별 + fit 파라미터 |
| 커밋된 결과의 provenance | ⚠️ `git_commit: "unknown"` | ✅ 커밋 해시 기록됨 (platform은 아직 컨테이너) |
| README 상태표 | ⚠️ Module 0 "실행 대기" | ✅ 갱신 |
| 자동 테스트 / CI | ❌ 없음 | ✅ `tests/` 22건, `.github/workflows/ci.yml` |
| 그림 파일 재현성 | ⚠️ 매번 다른 바이트 | ✅ 바이트 동일 (오답노트 #10) |
| 문서 번호 체계 | ⚠️ README/로드맵 어긋남 | ✅ Module 번호로 통일 |

---

## 1. 남은 보완 사항

### A1. ROOT arm 실행과 세 arm 대조 — **WSL2 필요, 최우선**
Session 2에서 정합성 문제 두 건(오답노트 #7, #8), Session 3에서 한 건(#11,
`pt_max` 누락)을 실행 전에 고쳤지만, 고친 것이 맞는지는 실행해야만 안다.
세 arm 모두 같은 스키마의 JSON을 쓰도록 바꿨고 `compare_arms`가 bin별로
비교하므로, 실행만 하면 판정이 나온다. **처음 실행이므로 API drift나 fit
tolerance 조정이 필요할 수 있다** — 그것도 오답노트에 기록할 내용.

### A2. 정식 환경에서 결과 재생성 후 커밋 — WSL2 필요
현재 `results/`는 클라우드 컨테이너(Python 3.11, numpy 2.5)에서 만들어졌다.
숫자는 환경에 무관하게 같음을 확인했지만, provenance의 `platform`이 WSL2
`hic` 환경을 가리키도록 A1과 함께 재생성해 커밋한다.

### B4. 남은 SKIP 두 건 정리 — WSL2 필요
`micromamba env update -f environment/environment.yml`로 `pyhepmc` 설치 후
verify 재실행 → `environment/versions_installed.txt` 갱신. LHAPDF는 Module 5.

### C2. 체계 오차 추정 틀 (Module 2 전후)
결과는 통계 오차뿐이다. Module 2(실데이터 비교)부터는 fit range·binning 변화에
따른 **체계 오차**를 기록하는 스크립트가 필요하다 (예: fit range 상·하한을
흔들어 T의 변화 폭을 표로).

### C3. `run_all.sh`가 `versions_installed.txt`를 덮어쓴다
어느 환경에서 돌리든 WSL2의 검증 기록을 덮어쓴다. 기록 파일은 WSL2에서만 쓰고,
다른 환경에서는 `--write` 없이 검증만 하도록 분리하는 것이 좋다.
(이번 세션에서는 덮어쓴 파일을 되돌렸다.)

---

## 2. Next step

### Step 1 — Module 1 닫기 (WSL2, 예상 1 세션)
```bash
micromamba activate hic
micromamba env update -f environment/environment.yml   # pyhepmc 반영
bash scripts/run_all.sh
```
**완료 기준**
- [ ] `verify_environment.py` FAIL 0, SKIP은 LHAPDF 1건만
- [ ] `pytest` 22건 통과 (WSL2에서도)
- [x] PyROOT arm, C++ macro 실행 성공 → `results/tables/m01_{pyroot,root}_summary.json`
- [x] `m01_compare_arms.py` → `ARMS AGREE` (2026-11-19, WSL2, ROOT 6.40.02)
- [x] 불일치 원인 기록: 오답노트 #15 (MINUIT 범위), #16 (TH1 소유권 segfault), #17 (C++ float32 곱)
- [x] `results/` 재생성분 커밋 (`eabd61f`, provenance platform = WSL2)
- [x] `docs/research_log.md`에 Session 4 기록, README 상태표에서 Module 1을 "완료"로

### Step 2 — Module 2: PYTHIA 8 minimum bias pp 13 TeV
`docs/modules/module02_pythia.md`에 상세. 진행:
- [x] `generators/pythia_mb_13TeV.cmnd` (빔, √s, `SoftQCD:inelastic`, Monash, PDF, tau0Max, seed, 이벤트 수)
- [x] `scripts/m02_run_pythia.py` — Module 1과 같은 TTree 스키마 + charge/process_code/n_ch_eta1
- [x] `config/analysis_m02_pythia.yaml`, `data/raw/README.md`(provenance 템플릿), `tests/test_m02_card.py`
- [x] WSL2에서 10 → 20,000 이벤트 생성 (PYTHIA 8.312; 카드 줄 순서 버그 → 오답노트 #14)
- [ ] Module 1 분석 코드를 **수정 없이** 실행 — 안 되면 인터페이스 설계 문제로 기록
- [x] "하전입자" 정의 → `particle_charge ≠ 0` (module02_pythia.md §5; 부수적으로 오답노트 #12 발견)
- [x] HEPData ALICE pp 13 TeV dN_ch/dη 표 확보 (ins1395253 v1 Table 1, 브라우저 수동 다운로드 →
      `--local` 검증, 2026-11-19). INEL>0 열 자동 선택, |η| < 1.8 18 bins
- [x] PYTHIA 샘플로 `m02_compare_data.py` 실행 → MC/data = 1.062 (|η| < 0.4), ±15 % 기준 통과, 약 2σ
- [x] 6 % 초과의 원인 (2026-11-19, `results/tables/m02/m02_sensitivity.json`):
      primary 정의는 ALICE와 일치(틀렸다면 +9.5 %), 이벤트 클래스도 일치. MPI pT0Ref +5 % → −3.4 %,
      CR 끔 → +37 %, 4C → +7 %. **MPI 튠(13 TeV로의 pT0 외삽)이 유력.** 별개로 INEL 정규화에는
      회절 모델 차이(f(INEL>0) MC 0.86 vs data 0.82)가 약 5 % 더해짐
- [x] Module 1 분석 코드를 PYTHIA 파일에 수정 없이 실행 — 돌아감. 빈 히스토그램 Z fit "수렴" 버그 수정 (#19)
- [x] pT 스펙트럼 비교 (ALICE Table 2, 2026-11-19): ⟨pT⟩ MC/data 1.03, MC/data 0.86 (0.15 GeV/c) →
      1.24 (0.5 GeV/c). pT0Ref는 크기만, CR은 모양을 움직임 → MPI·CR 동시 튜닝 필요 (module02 §8)
- [x] pT0Ref × CR range 2차원 스캔 (module02 §9): dN/dη 최적 (2.52, 1.4), 두 상관 모형 모두.
      pT 모양은 격자 밖 → 다음 후보 `StringPT:sigma`
- [x] WSL2에서 `bash scripts/run_all_pythia.sh` — Module 2–4 전부 PYTHIA 8.312로 재생성 (`1cc643e`)
- [x] pT 모양: `StringPT:sigma` 스캔 (module02 §11) — σ ≈ 0.31이 ⟨pT⟩를 맞추지만 0.4–0.8 GeV/c 초과는 남음
- [ ] ALICE pp 13 TeV identified π/K/p 스펙트럼과 비교 (0.4–0.8 GeV/c 초과의 원인 판별)
- [x] `analysis/m02_compare_data.py`: 이벤트 클래스 선택 + 데이터 binning + ratio panel + 판정 (toy로 self-comparison 테스트, 실데이터는 아직)

**완료 기준**: dN_ch/dη|₀가 ALICE 측정의 ±15% 이내, 비교 그림에 ratio panel.

### Step 3 — Module 3 HepMC ✅ (2026-11-19)
`docs/modules/module03_hepmc.md`. 남은 것: WSL2에서 같은 두 명령으로 재현 (pyhepmc 설치 확인).

### Step 4 — Module 4 FastJet ✅ (2026-11-19)
`docs/modules/module04_fastjet.md`.

### 그 다음
로드맵: Module 5 LHAPDF·nPDF → 6 Glauber MC → 7 Angantyr p–Pb/Pb–Pb (Module 4 배경 차감을
실제 중이온 배경으로).

### Step 5 — Module 5 PDF ✅ (2026-11-19)
`docs/modules/module05_pdfs.md`. 남은 것: WSL2 재현 (`bash scripts/run_all_pythia.sh stringpt m05`),
LHAPDF로 EPPS21 교체.

### Step 6 — Module 6 Glauber ✅ (2026-11-19)
`docs/modules/module06_glauber.md`. σ_PbPb 7.64 b, ALICE 중심도 표와 일치.
남은 것: 변형당 40k로 계통 재계산, 다중도(NBD) 기반 중심도.

### Step 7 — identified π/K/p ✅ (2026-11-19)
module02 §12. 0.4–0.8 GeV/c 초과의 95 %가 π. K/π 0.92, K·p는 0.8 GeV/c 위에서 15–25 % 부족.
WSL2 재현 ✅. π 출처 분해 ✅ (module02 §13): 조성(벡터 중간자 비)·Bose–Einstein은 모양을 못 바꾸고,
강입자 재산란이 모양 비를 1.38 → 1.23 (조정과 합쳐 1.19)으로.
남은 것: WSL2 재현 (`bash scripts/run_all_pythia.sh origin`), rope/shoving 시험, 다중도별 K/π·⟨pT⟩
(Table 17, 20; V0M 수용 범위로 MC 클래스 정의) — 재산란은 고다중도에서 커지므로 다중도 의존성이 판별 시험.

### 그 다음
Module 7: PYTHIA Angantyr p–Pb / Pb–Pb — Module 6의 N_part 분포와 비교, Module 4 배경 차감을
실제 중이온 배경으로.
