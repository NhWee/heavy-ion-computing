# 로드맵 — 다음 5개 마일스톤

각 마일스톤은 "돌아가는 코드 + 검증 + 문서"까지를 하나로 봅니다.
설치는 필요해질 때만 합니다.

---

## M1 ✅ ROOT 파일 I/O와 첫 히스토그램 — **완료(uproot arm)**

- toy 이벤트 20,000개 → TTree → pT/η/다중도/불변질량 + fit + closure test
- 남은 것: **WSL2에서 ROOT arm 실행**하여 uproot 결과와 숫자 일치 확인
- 산출물: `results/figures/m01_*.pdf`, `docs/modules/module01_root.md`

---

## M2 — 환경 실체화: WSL2 `hic` 환경

**할 일**
1. `bash environment/setup_wsl2.sh`
2. `python environment/verify_environment.py --write environment/versions_installed.txt`
3. `bash scripts/run_all.sh` — ROOT arm까지 전부 실행

**검증 기준**
- ROOT: TFile 생성·재읽기 100 entries, RDataFrame 합계 일치
- PYTHIA: pp 13 TeV 10 이벤트, 최종상태 입자 > 100개
- FastJet: anti-kT R=0.4 → 예상 jet 수
- ROOT arm의 `<N_ch>`, `dN/dη|₀`, Z 질량이 uproot arm과 **소수점 둘째 자리까지 일치**

**배우는 것**: conda 환경 관리, `root-config`, `LD_LIBRARY_PATH`,
`PYTHONPATH`가 왜 HEP에서 늘 문제인지.

---

## M3 — PYTHIA 8: pp minimum bias와 실제 데이터 비교

**할 일**
1. `generators/pythia_mb_13TeV.cmnd` — 빔·에너지·process·tune·시드·이벤트 수 전부 명시
2. 10 → 1,000 → 20,000 이벤트 단계적 확대
3. Module 1과 **같은 분석 코드**로 `dN_ch/dη`, pT 스펙트럼, 다중도
4. HEPData에서 ALICE pp 13 TeV `dN_ch/dη` (INEL>0) 표를 내려받아 비교 플롯 + ratio panel

**검증 기준**
- toy 생성기를 PYTHIA 출력으로 바꿨을 때 분석 스크립트가 **수정 없이** 동작
- `dN_ch/dη|₀`가 ALICE 측정의 ±15 % 이내 (tune 차이를 넘어서면 설정 오류 의심)

**배우는 것**: hard scattering → ISR → MPI → FSR → hadronization → decay의 각
단계를 켜고 끄며 관측량이 어떻게 변하는지. PYTHIA는 블랙박스가 아님.

---

## M4 — HepMC + PDG: 이벤트 레코드 해부

**할 일**
1. PYTHIA 출력을 HepMC3로 저장 (`pyhepmc`로 읽기)
2. 이벤트 하나를 골라 parton → shower → hadron → decay 사슬을 **직접 추적**
3. status code, mother/daughter, vertex 구조를 그림으로 정리
4. 보존량 검사: 각 vertex에서 4-운동량 보존, 전체 전하 보존

**검증 기준**
- 임의 이벤트 100개에서 vertex별 4-운동량 보존이 수치 오차(1e-6 상대) 이내
- 최종상태 입자 수와 status==1 개수 일치

---

## M5 — FastJet: 제트 재구성 첫 단계

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

M6 LHAPDF·nuclear PDF → M7 Glauber MC 직접 구현 → M8 Angantyr/HIJING/AMPT로
pPb·PbPb → M9 CERN Open Data / HEPData 실측 재현 → M10 jet quenching(R_AA,
dijet asymmetry) → M11 JETSCAPE/JEWEL/MUSIC.

최종 목표: **공개 데이터와 MC로 LHC 중이온 관측량 하나를 재현하고 모델 파라미터를
바꿔가며 기술이 개선되는지 평가하기.**
