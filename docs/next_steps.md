# 보완 사항과 Next step — 2026-11-05 점검

마지막 작업(Session 2, 2026-08-27) 이후 약 10주 만의 점검.
저장소 전체(코드·config·결과·문서)를 읽고, 클라우드 컨테이너에서 uproot arm을
처음부터 다시 실행해 확인했다.

---

## 0. 점검 결과 요약

| 항목 | 상태 |
|---|---|
| uproot arm 재현성 | ✅ 같은 seed로 재생성 시 모든 수치가 커밋된 `m01_fit_results.json`과 일치 (Z fit 파라미터만 상대 ~1e-9 수준의 최적화기 잡음) |
| Closure test | ✅ 종류별 pull 1.5 / 0.8 / 0.4σ, inclusive는 의도대로 20σ 불일치 |
| ROOT arm (PyROOT, C++ macro) | ❌ **아직 한 번도 실행되지 않음** — Module 1의 유일한 미완 항목 |
| `m01_compare_arms.py` | ⚠️ 작성됨, 실행된 적 없음 |
| 커밋된 결과의 provenance | ⚠️ `git_commit: "unknown"`, 생성 환경이 WSL2 `hic` 환경이 아닌 클라우드 컨테이너 (python 3.11.15) |
| README 상태표 | ⚠️ Module 0을 "사용자 실행 대기"로 표시 — 실제로는 Session 2에서 23 PASS 완료 |

---

## 1. 보완할 부분

우선순위 순. **A = Module 2 진입 전에 반드시**, B = 곧, C = 여유 있을 때.

### A1. ROOT arm 실행과 두 arm 대조 (Module 1 마무리)
Session 2에서 정합성 문제 두 건(오답노트 #7, #8)을 실행 전에 고쳤지만, 고친 것이
맞는지는 실행해야만 안다. "작성했지만 돌려보지 않은 분석 코드"는 검증되지 않은 것.

### A2. 정식 환경에서 결과 재생성 후 다시 커밋
현재 `results/`는 클라우드 컨테이너에서 git 커밋이 기록되지 않은 채 만들어졌다.
WSL2 `hic` 환경에서 `run_all.sh`를 돌려, `git_commit`이 채워진 결과로 교체해야
"이 그림은 이 커밋·이 환경에서 나왔다"가 성립한다.

### A3. `compare_arms`가 대조하는 범위가 좁다
- 지금은 **스칼라 5개**(이벤트 수, ⟨N_ch⟩, dN/dη, 이중뮤온 후보 수, ⟨m⟩)만 비교.
  히스토그램 **bin별 내용과 오차**, **Tsallis fit 파라미터**는 비교하지 않는다.
  스칼라가 일치해도 binning이나 정규화 버그는 그대로 남을 수 있다.
- C++ macro(`m01_analysis_root.C`)는 콘솔 출력과 `.root`만 남기고 JSON 요약을
  쓰지 않으므로 **대조 대상에서 빠져 있다**. 실제로는 세 개의 arm인데 두 개만 확인.
- 보완: 세 arm 모두 같은 스키마의 JSON(bin 내용 포함)을 쓰고, `compare_arms`가
  bin별 상대 차이와 fit 파라미터의 σ 단위 차이를 함께 보고.

### A4. `run_all.sh`가 환경 검증 실패를 삼킨다
`verify_environment.py ... || true` 때문에 FAIL이 나도 파이프라인이 계속 진행된다.
오답노트의 원칙("exit 0이었으니 됐다를 믿지 않는다")과 반대 방향.
`verify_environment.py`는 이미 FAIL이 있으면 종료 코드 1을, SKIP만 있으면 0을
반환하므로 `|| true`만 지우면 된다 (FAIL은 중단, SKIP은 허용).

### B1. 자동 테스트와 CI가 없다
`kinematics` 단위 검증(round-trip, back-to-back 뮤온 → 91.2 GeV)은 Session 1에서
손으로 확인했지만 `tests/`로 남아 있지 않다. Module 2부터 코드가 늘어나면 회귀를
잡을 장치가 필요하다.
- `tests/test_kinematics.py`, `tests/test_toy_generator.py` (작은 seed 고정 샘플)
- closure test를 소규모(예: 2,000 이벤트)로 pytest에 연결
- GitHub Actions: uproot arm + 테스트만 (ROOT 없이도 돌아가는 부분)

### B2. 그림 파일이 바이트 단위로 재현되지 않는다
같은 seed·같은 숫자인데도 PDF/PNG가 매번 다른 바이너리로 나온다(생성 시각 등의
메타데이터, 라이브러리 버전). 그래서 `git diff`에 의미 없는 변경이 생긴다.
보완: `savefig(..., metadata={"CreationDate": None})` 및 `SOURCE_DATE_EPOCH` 고정.

### B3. 문서 번호 체계가 둘이다
README는 "Module 2 = PYTHIA", 로드맵은 "M2 = 환경, M3 = PYTHIA".
하나로 통일해야 커밋 메시지·Notion 로그가 헷갈리지 않는다.
(제안: Module 번호를 기준으로 하고, 로드맵의 M 번호는 Module 번호와 같게.)

### B4. 남은 SKIP 두 건 정리
`pyhepmc`는 `environment.yml`에 추가했지만 아직 설치 확인 전 (SKIP 상태로 기록됨).
Module 4 전에 `micromamba env update` 후 verify 재실행. LHAPDF는 Module 6까지 보류.

### C1. Closure test 범위 확대
현재는 Tsallis (T, n)만 확인. 생성기 입력에서 해석적으로 예측 가능한
- 다중도: NBD (n̄, k) → acceptance 안의 기댓값
- dN/dη: plateau 모양 → |η|<0.5 기대값
- Z: acceptance 보정 후 질량/폭
도 closure에 넣으면 "분석이 생성기를 되찾는다"는 주장이 관측량 4종 전체로 확장된다.

### C2. 불확실성 표기
결과는 통계 오차뿐이다. toy MC라 당연하지만, Module 3(실데이터 비교)부터는 fit range·
binning 변화에 따른 **체계 오차 추정 틀**을 미리 만들어 두는 것이 좋다
(예: fit range를 흔들어 T의 변화를 기록하는 스크립트).

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
- [ ] PyROOT arm, C++ macro 실행 성공
- [ ] `m01_compare_arms.py` → `ARMS AGREE`
- [ ] 불일치가 나오면 원인을 `docs/troubleshooting.md`에 #10으로 기록
- [ ] `results/` 재생성분 커밋 (`git_commit` 필드가 채워졌는지 확인)
- [ ] `docs/research_log.md`에 Session 3 기록, README 상태표 갱신

### Step 2 — 검증 장치 보강 (A3, A4, B1, B2; 1 세션)
Module 2에서 생성기를 PYTHIA로 바꾸기 **전에** 해야 한다. 생성기를 바꾼 뒤에 숫자가
달라지면, 분석 버그인지 생성기 차이인지 구별할 기준이 지금의 toy closure뿐이다.

### Step 3 — Module 2: PYTHIA 8 minimum bias pp 13 TeV
`docs/roadmap.md`의 PYTHIA 항목 그대로. 시작점:
1. `generators/pythia_mb_13TeV.cmnd` 작성 (빔, √s, `SoftQCD:inelastic`, tune, seed, 이벤트 수)
2. 10 → 1,000 → 20,000 이벤트로 확대, Module 1과 **같은 TTree 스키마**로 저장
3. Module 1 분석 코드를 **수정 없이** 실행 — 수정이 필요하면 그것이 인터페이스 설계 문제
4. HEPData에서 ALICE pp 13 TeV dN_ch/dη (INEL>0) 표를 받아 비교
   (DOI·표 번호·선택 조건을 `data/raw/<dataset>.json`에 먼저 기록)

**완료 기준**: dN_ch/dη|₀가 ALICE 측정의 ±15% 이내, 비교 그림에 ratio panel.

### 그 다음
로드맵 순서대로 HepMC 이벤트 레코드 → FastJet.
