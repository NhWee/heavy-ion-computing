# 로드맵 — 다음 마일스톤

각 마일스톤은 "돌아가는 코드 + 검증 + 문서"까지를 하나로 봅니다.
설치는 필요해질 때만 합니다. 번호는 README의 **Module 번호와 같습니다**
(2026-11-05에 통일; 이전에는 이 문서의 M 번호가 하나씩 어긋나 있었습니다).
당장의 할 일과 완료 기준은 [`next_steps.md`](next_steps.md)에 있습니다.

---

## Module 0 ✅ 환경 — WSL2 `hic` 환경 (완료, 2026-08-27)

micromamba 2.9.0 / ROOT 6.40.02 / PYTHIA 8 / FastJet 3.5.1.4, 16 논리코어.
`verify_environment.py`: 23 PASS / 0 FAIL / 2 SKIP (`pyhepmc` 추가 후 재검증 필요,
`lhapdf`는 Module 5까지 보류).

---

## Module 1 — ROOT 파일 I/O와 첫 히스토그램 (uproot arm 완료, ROOT arm 대기)

- ✅ toy 이벤트 20,000개 → TTree → P(N_ch) / dN/dη / 불변 pT 수율 / m(μμ) + fit
- ✅ closure test: 다중도·dN/dη(해석적 기대값)·종류별 Tsallis (T, n) 복원
- ✅ `tests/` 22건 + GitHub Actions CI, 그림 바이트 단위 재현
- ⬜ **WSL2에서 ROOT arm 두 개 실행** → `m01_compare_arms.py`가 세 arm의
  스칼라·히스토그램 bin·fit 파라미터 일치를 판정 → `ARMS AGREE`
- ⬜ 재생성한 `results/`를 WSL2 환경에서 커밋 (provenance의 platform 갱신)

**검증 기준**: 세 arm의 `<N_ch>`, `dN/dη`, 이중뮤온 수가 상대 1e-6 이내,
히스토그램 bin이 float round-off 또는 경계 입자 1개 이내, Tsallis (T, n)이
0.05σ 이내.

---

## Module 2 — PYTHIA 8: pp minimum bias와 실제 데이터 비교

**할 일**
1. `generators/pythia_mb_13TeV.cmnd` — 빔·에너지·process·tune·시드·이벤트 수 전부 명시
2. 10 → 1,000 → 20,000 이벤트 단계적 확대, Module 1과 **같은 TTree 스키마**로 저장
3. Module 1과 **같은 분석 코드**로 `dN_ch/dη`, pT 스펙트럼, 다중도
4. HEPData에서 ALICE pp 13 TeV `dN_ch/dη` (INEL>0) 표를 내려받아 비교 플롯 + ratio panel
   (받기 전에 DOI·표 번호·선택 조건을 `data/raw/<dataset>.json`에 기록)

**검증 기준**
- toy 생성기를 PYTHIA 출력으로 바꿨을 때 분석 스크립트가 **수정 없이** 동작
  (CI의 `tests/test_m01_pipeline.py`가 그대로 통과)
- `dN_ch/dη|₀`가 ALICE 측정의 ±15 % 이내 (tune 차이를 넘어서면 설정 오류 의심)

**배우는 것**: hard scattering → ISR → MPI → FSR → hadronization → decay의 각
단계를 켜고 끄며 관측량이 어떻게 변하는지. PYTHIA는 블랙박스가 아님.

---

## Module 3 — HepMC + PDG: 이벤트 레코드 해부

**할 일**
1. PYTHIA 출력을 HepMC3로 저장 (`pyhepmc`로 읽기)
2. 이벤트 하나를 골라 parton → shower → hadron → decay 사슬을 **직접 추적**
3. status code, mother/daughter, vertex 구조를 그림으로 정리
4. 보존량 검사: 각 vertex에서 4-운동량 보존, 전체 전하 보존

**검증 기준**
- 임의 이벤트 100개에서 vertex별 4-운동량 보존이 수치 오차(1e-6 상대) 이내
- 최종상태 입자 수와 status==1 개수 일치

---

## Module 4 — FastJet: 제트 재구성 첫 단계

**할 일**
1. PYTHIA hard QCD 이벤트(pTHat > 20 GeV) 생성
2. anti-kT, kT, Cambridge/Aachen을 R = 0.2/0.4/0.6으로 클러스터링
3. 제트 pT 스펙트럼, leading jet, dijet $\Delta\phi$, jet mass
4. 같은 이벤트에 toy underlying event를 얹어 제트 pT가 어떻게 오염되는지 →
   jet area 기반 배경 차감의 필요성 확인

**검증 기준**
- anti-kT 제트가 "원형"이고 kT/CA는 불규칙 — jet area 분포로 정량 확인
- 배경을 넣기 전/후 leading jet pT 차이가 $\rho \cdot A$ 예측과 일치

**배우는 것**: sequential recombination의 거리 척도 $d_{ij}$, IRC safety,
중이온 환경에서 제트 정의가 왜 까다로운지 (Module 10의 전제).

---

## 그 다음 (개요)

Module 5 LHAPDF·nuclear PDF → 6 Glauber MC 직접 구현 → 7 Angantyr/HIJING/AMPT로
pPb·PbPb → 8 CERN Open Data / HEPData 실측 재현 → 9 jet quenching(R_AA,
dijet asymmetry) → 10 JETSCAPE/JEWEL/MUSIC.

> 이전 문서·코드 주석에 남아 있는 "Module 6 LHAPDF", "Module 8·10 flow/jet",
> "Module 9 open data"는 이 개요 번호 기준으로 각각 5, 7·9, 8에 해당합니다.
> 해당 모듈에 도착할 때 주석을 고칩니다.

최종 목표: **공개 데이터와 MC로 LHC 중이온 관측량 하나를 재현하고 모델 파라미터를
바꿔가며 기술이 개선되는지 평가하기.**
