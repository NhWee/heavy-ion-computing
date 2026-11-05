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
- [ ] PyROOT arm, C++ macro 실행 성공 → `results/tables/m01_{pyroot,root}_summary.json`
- [ ] `m01_compare_arms.py` → `ARMS AGREE`
- [ ] 불일치가 나오면 원인을 `docs/troubleshooting.md`에 #12부터 기록
- [ ] `results/` 재생성분 커밋 (provenance `platform`이 WSL2인지 확인)
- [ ] `docs/research_log.md`에 Session 4 기록, README 상태표에서 Module 1을 "완료"로

### Step 2 — Module 2: PYTHIA 8 minimum bias pp 13 TeV
`docs/modules/module02_pythia.md`에 상세. 진행:
- [x] `generators/pythia_mb_13TeV.cmnd` (빔, √s, `SoftQCD:inelastic`, Monash, PDF, tau0Max, seed, 이벤트 수)
- [x] `scripts/m02_run_pythia.py` — Module 1과 같은 TTree 스키마 + charge/process_code/n_ch_eta1
- [x] `config/analysis_m02_pythia.yaml`, `data/raw/README.md`(provenance 템플릿), `tests/test_m02_card.py`
- [ ] WSL2에서 10 → 1,000 → 20,000 이벤트 생성 (스크립트는 **미실행**, API drift 가능)
- [ ] Module 1 분석 코드를 **수정 없이** 실행 — 안 되면 인터페이스 설계 문제로 기록
- [x] "하전입자" 정의 → `particle_charge ≠ 0` (module02_pythia.md §5; 부수적으로 오답노트 #12 발견)
- [ ] HEPData에서 ALICE pp 13 TeV dN_ch/dη (INEL>0) 표를 받아 비교
      (DOI·표 번호·선택 조건을 `data/raw/<dataset>.json`에 먼저 기록)
- [ ] `analysis/m02_compare_data.py`: INEL>0 선택 + ratio panel

**완료 기준**: dN_ch/dη|₀가 ALICE 측정의 ±15% 이내, 비교 그림에 ratio panel.

### 그 다음
로드맵 순서대로 Module 3 HepMC 이벤트 레코드 → Module 4 FastJet.
