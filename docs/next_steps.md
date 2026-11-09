# 보완 사항과 Next step — 2026-11-09 갱신

- 2026-11-05: 10주 만의 점검과 같은 날 1차 반영(브랜치 `claude/vigilant-sagan-aatxhs`)
- 2026-11-07~09: 2차 반영(브랜치 `claude/confident-bohr-ie66ff`, `docs/research_log.md`
  Session 2026-11-09)

> **먼저 §3을 읽으세요.** 사용자 로컬 저장소에는 아직 push되지 않은 Session 3·4
> 작업이 있고, GitHub의 두 Claude 브랜치는 그 작업을 모른 채 `main`에서 갈라져
> 나왔습니다. 같은 버그를 양쪽에서 따로 고친 곳이 여럿입니다.

---

## 0. 상태 요약

| 항목 | 11-05 점검 | 11-05 1차 반영 | 11-09 2차 반영 |
|---|---|---|---|
| uproot arm 재현성 | ✅ | ✅ | ✅ 같은 기계에서 표 4개·그림 전부 동일, CI에서 중심값 1e-6 이내 |
| fit 방법 | 검사 안 함 | 검사 안 함 | ✅ Poisson likelihood로 교체. pseudo-experiment로 편향·포함률 측정 (#15, #16) |
| Z 질량 해석 | ❌ 틀림 (acceptance sculpting) | ❌ 그대로 | ✅ 철회·정정 (#15) |
| Closure test | Tsallis만 | + 다중도, dN/dη | + NBD k, truth Z 선형, 참값은 파일 provenance, JSON 기록 |
| ROOT arm | ❌ 미실행 | ❌ 미실행 (코드 수정) | ❌ 미실행. 첫 실행에서 터졌을 use-after-free 등 수정 (#19) |
| `m01_compare_arms.py` | 스칼라 5개 | 세 arm, bin별, fit | + 같은 입력·같은 컷 확인, arm 누락 = 실패, 판정 테스트 6건 |
| CI | 없음 | ⚠️ 추가했지만 **11회 전부 설정 단계에서 실패** | ✅ 실제로 테스트 65건과 전체 재현이 돈다 (#14) |
| 체계 오차 틀 (C2) | 없음 | 없음 | ✅ `m01_fit_systematics.py`: bootstrap 유의성, 참 모형 대조군 |
| 버전 기록 덮어쓰기 (C3) | 열림 | 열림 | ✅ `RECORD_ENV=1`일 때만 기록 |
| provenance | `git_commit: unknown` | 해시 기록 | ✅ 코드·설정이 커밋과 다르면 `-dirty` |
| Z 값 출처 표기 | "PDG 2024" (틀림) | 그대로 | ✅ LEP 시대 값과 RPP 2025/2026 구분 (`hicpy.pdg`) |

---

## 1. 남은 보완 사항

### A1. ROOT arm 실행과 세 arm 대조 — **WSL2 필요, 최우선**
ROOT C++ 매크로와 PyROOT arm은 이 계열에서 한 번도 실행되지 않았습니다. 2026-11-07 수정분도
`g++ -fsyntax-only`(stub 헤더)로 문법만 확인했습니다. 사용자 로컬 Session 3에서는 실행해
`ARMS AGREE`를 얻었으므로, §3대로 그쪽 ROOT 코드와 비교하는 것이 가장 빠릅니다.
`m01_compare_arms.py`는 이제 ROOT arm 둘 다 있어야 통과합니다.

### A2. 정식 환경에서 결과 재생성 후 커밋 — WSL2 필요
`RECORD_ENV=1 bash scripts/run_all.sh` → `python scripts/check_results_drift.py` → 커밋.
현재 `results/`의 provenance는 클라우드 컨테이너(Python 3.13)를 가리킵니다.

### B4. 남은 SKIP 정리 — WSL2 필요
`micromamba env update -f environment/environment.yml`로 `pyhepmc` 설치 후 verify 재실행.

### B5. 생성기 η 샘플러의 경계 결함 (#17)
CDF가 0에서 시작하지 않아 입자 약 20개가 정확히 η = −6.0에 놓입니다. 고치면 모든 η가
미세하게 움직이므로, 다음에 표본을 의도적으로 다시 만들 때(Module 2 착수 시점) 함께 고칩니다.

### B6. Module 2의 Z 분석에는 배경 항이 필요
toy에는 이중뮤온 배경이 없어 신호만으로 맞춥니다. PYTHIA 표본(Drell–Yan, heavy-flavour
뮤온)에서는 배경 모형을 넣고, `m01_fit_validation.py`의 toy에도 같은 배경을 넣어
편향을 다시 재야 합니다. 0 이상 제한이 Γ를 −0.33σ 편향시킨다는 점에 유의하세요.

### C4. CI 결과를 실제로 본다
push 뒤 Actions 실행 결과를 확인하는 습관이 필요합니다(#14). GitHub 알림 설정에서 실패 메일을
켜 두는 것을 권합니다.

### C5. (선택) 오차의 기계 의존성
MINUIT HESSE 오차는 기계마다 상대 1e-5–2e-3 다릅니다. 물리적으로는 무의미하지만, 논문용 숫자는
MINOS 비대칭 오차나 `strategy=2`로 확인하는 것이 좋습니다.

---

## 2. Next step

### Step 1 — 로컬 작업과 이 브랜치 합치기 (§3)

### Step 2 — Module 1 닫기 (WSL2)
```bash
micromamba activate hic
micromamba env update -f environment/environment.yml   # pyhepmc 반영
RECORD_ENV=1 bash scripts/run_all.sh
python scripts/check_results_drift.py
```
**완료 기준**
- [ ] `verify_environment.py` FAIL 0, SKIP은 LHAPDF 1건만
- [ ] `pytest` 65건 통과 (WSL2에서도)
- [ ] PyROOT arm, C++ macro 실행 성공 → `results/tables/m01_{pyroot,root}_summary.json`
- [ ] `m01_compare_arms.py` → `ARMS AGREE` (같은 입력·같은 컷 확인 포함)
- [ ] 불일치가 나오면 원인을 `docs/troubleshooting.md`에 다음 번호로 기록
- [ ] `results/` 재생성분 커밋 (provenance `platform`이 WSL2인지 확인)

### Step 3 — Module 2: PYTHIA 8 minimum bias pp 13 TeV
`docs/modules/module02_pythia.md`에 상세합니다. B5(η 샘플러)와 B6(Z 배경)을 함께 처리합니다.
**완료 기준**: dN_ch/dη|₀가 ALICE 측정의 ±15 % 이내, 비교 그림에 ratio panel.

---

## 3. 로컬 작업과 GitHub 브랜치의 관계 — 합치기 전에 꼭 읽기

### 세 계열

| 계열 | 위치 | 마지막 | 내용 |
|---|---|---|---|
| `main` | GitHub | `7d84a31` (08-28) | Session 2까지 |
| 로컬 Session 3–4 | 사용자 PC, Notion 작업 로그 | `d3379f8` (08-31) + 09-26 수정분 (push 전) | ROOT arm 실행·`ARMS AGREE`, float32 버그, Z Poisson fit, 독립 검증 `docs/audit_2026-09-26.md` |
| `claude/vigilant-sagan-aatxhs` | GitHub | `5ed1a17` (11-05) | 11-05 점검, 세 arm 비교 틀, closure 확장, 테스트, CI(실패), Module 2 착수 |
| `claude/confident-bohr-ie66ff` | GitHub | vigilant-sagan + 2차 반영 | fit 방법 교정, 방법 검증, 체계 오차, ROOT arm 검토, CI 실가동 |

Notion의 2026-09-26 기록에는 다음 단계가 "WSL2에서 수정본 재실행 → GitHub push"로 남아 있습니다.
GitHub에서는 그 작업을 볼 수 없었으므로, 아래 대조표는 **Notion 기록만 보고** 만든 것입니다.

### 양쪽에서 따로 고친 것

| 문제 | 로컬 (Notion 기록) | GitHub 브랜치 | 합칠 때 |
|---|---|---|---|
| float32 가중 히스토그램 | #11, `d3379f8` | #12 (vigilant-sagan) | 같은 해결(float64). 하나로 |
| 저통계 Z fit + 해석 오류 | #13: Poisson LL, BW + 상수 배경, toy 300회 | #15: Poisson LL, BW만, toy 1000회 | **결과가 같다**(91.057 ± 0.105, 2.28 ± 0.23). toy에는 배경이 없으므로 BW만이 맞는 모형. 0 이상 배경 항은 Γ를 −0.33σ 편향시킨다(이 브랜치의 측정) |
| 세 arm 대조의 허점 | #14: 입력 없으면 실패, 히스토그램 4종 × 2 arm + fit, double + asinh, TH1D, PyROOT 리터럴, float32 η 테스트 | vigilant-sagan 틀 + #19 (use-after-free 등) | 대부분 겹침. **ROOT 쪽 코드는 실제로 실행된 로컬 판을 기준으로** 하고, 이 브랜치의 provenance·컷 확인과 판정 테스트를 옮긴다 |
| 설치 스크립트 셸 hook | #10, `d19b5bb` (`d3379f8`이 되돌렸다가 Session 4에서 복원) | 없음 (vigilant-sagan #13은 수동 해결만 기록) | 로컬 판을 쓴다 |
| Z 값 출처 | "PDG 2022" 정정, 2024 = 91.1880/2.4955 | `hicpy.pdg`: RPP 2025 = 91.1880, RPP 2026 = 91.1879 | 일치 |
| ALICE dN/dη 인용 | 6.46 ± 0.19, 1차 출처 확인 | `data/raw/alice_pp13TeV_dndeta.json` (unverified) | 로컬 확인값으로 JSON을 채운다 |
| 연구 노트 Session 1 | 원래 값 복원, 철회 표시 | 철회 표시만 추가 | 문장 단위로 합친다 |
| `git describe --dirty` | 있음 | `hicpy.provenance` (코드·설정만 dirty 판정) | 하나로 |

### 번호가 겹친다
- 오답노트: 로컬 #10–#14와 GitHub #10–#19가 **같은 번호로 다른 내용**을 가리킵니다.
- 세션: Notion의 Session 3(08-31)·4(09-26)와 이 저장소 `research_log.md`의 Session 3(11-05)이 다릅니다.
- 합친 뒤 한 번호 체계로 다시 매기고, Notion 실패 보고의 "가드 위치"도 갱신하세요.

### 권하는 절차 (오답노트 #12 "사용자 폴더가 원본" 규칙을 지키며)
1. PC에서 Session 4 수정분을 커밋하고(아직이라면) **새 브랜치로만** push합니다.
   아직 아무것도 합치지 않습니다.
   ```bash
   git push origin HEAD:local/session4
   ```
2. 두 계열을 파일 단위로 비교합니다.
   ```bash
   git fetch origin
   git diff --stat origin/local/session4 origin/claude/confident-bohr-ie66ff
   ```
3. ROOT 쪽 코드(`m01_analysis_root.C`, `m01_analysis_pyroot.py`, `setup_wsl2.sh`)는 실행된 로컬 판에서
   시작하고, 이 브랜치에서 `src/hicpy/`(fitting, histograms, events, pdg, provenance), `m01_fit_validation.py`,
   `m01_fit_systematics.py`, closure 확장, CI 수정, `check_results_drift.py`, `tests/`를 가져옵니다.
4. WSL2에서 `RECORD_ENV=1 bash scripts/run_all.sh` → `ARMS AGREE` → 결과 커밋 → `main`으로 PR.
5. 오답노트·세션 번호를 하나로 정리하고 Notion을 갱신합니다.
