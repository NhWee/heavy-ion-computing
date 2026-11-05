# Module 2 — PYTHIA 8: pp minimum bias와 실제 데이터 비교

상태(2026-11-05): **카드·생성 스크립트·설정 작성 완료, 실행 대기** (WSL2 `hic` 환경).

---

## 0. 이 모듈이 답하는 질문

> Module 1의 toy 생성기를 "진짜" 생성기로 바꾸면 분석 체인은 그대로 돌아가는가,
> 그리고 그 결과는 실측 $dN_{\rm ch}/d\eta$와 얼마나 맞는가?

처음으로 QCD가 들어옵니다. 다만 이 모듈의 PYTHIA는 **파라미터를 건드리는 대상**이지
블랙박스가 아닙니다 — 단계를 하나씩 끄며 관측량이 어떻게 변하는지 봅니다.

---

## 1. 카드 — `generators/pythia_mb_13TeV.cmnd`

| 설정 | 값 | 왜 |
|---|---|---|
| `Beams:idA/idB`, `eCM` | 2212 / 2212 / 13000 | LHC Run 2 pp; ALICE·CMS·ATLAS 13 TeV 데이터와 비교 |
| `SoftQCD:inelastic` | on | ND + SD + DD + CD = 실험의 INEL 클래스. **INEL>0·NSD는 분석에서** 고른다 |
| `Tune:pp` | 14 (Monash 2013) | PYTHIA 8.3 기본 tune. PDF·MPI·string·CR을 한 세트로 맞춘 것이라 PDF만 바꾸면 tune이 깨짐 |
| `PDF:pSet` | 13 (NNPDF2.3 LO) | Monash가 맞춰진 PDF. 기본값이지만 **보이게** 적음 |
| `ParticleDecays:tau0Max` | 10 mm/c | $c\tau < 1$ cm인 입자는 붕괴시킴 → ALICE "primary" 정의와 같은 입자 목록 |
| `PartonLevel:MPI/ISR/FSR`, `HadronLevel:*`, `ColourReconnection` | on | 전부 기본값. Module 2 연습은 이것들을 하나씩 **끄는** 것 |
| `Random:seed` | 20260827 | Module 1과 같은 규약. 파일 안에도 기록 |

카드에 없는 것: 이벤트 클래스 선택, 운동학 cut. 둘 다 분석의 몫입니다.
생성 단계에서 자르면 한 샘플로 INEL / INEL>0 / NSD를 비교할 수 없습니다.

---

## 2. 생성 스크립트 — `scripts/m02_run_pythia.py`

Module 1과 **같은 TTree 브랜치**(`nparticle`, `particle_px/py/pz/E/pdg`, `has_z`)를
쓰고, 거기에 Module 2가 필요한 브랜치를 **더합니다**:

| 브랜치 | 용도 |
|---|---|
| `particle_charge` | 하전입자 선택 (`≠ 0`). Module 1 분석도 이 브랜치를 쓴다 (§5) |
| `process_code` | 101 ND, 103/104 SD, 105 DD, 106 CD → NSD 선택 |
| `n_ch_eta1` | \|η\|<1 하전입자 수 → INEL>0 선택 |
| `weight` | 이벤트 가중치 (MB는 1) |

provenance에는 카드 **전문**, PYTHIA 버전, seed, `sigmaGen ± err` (mb), 시도/수락
이벤트 수, git 커밋이 들어갑니다. 단면적이 있어야 "이벤트 수"를 "mb"로 바꿔
데이터와 정규화를 맞출 수 있습니다.

`--hepmc`를 주면 HepMC3 파일도 씁니다 (Module 3용; 이 pythia8 빌드에 HepMC
인터페이스가 없으면 경고만 내고 무시).

---

## 3. 실행 순서 (WSL2)

```bash
micromamba activate hic
# 1. 작동 확인: 10 이벤트
python scripts/m02_run_pythia.py --n-events 10 --out data/generated/m02_test10.root
# 2. 통계 확보: 1,000 -> 20,000
python scripts/m02_run_pythia.py --n-events 1000 --out data/generated/m02_test1k.root
python scripts/m02_run_pythia.py                      # 카드의 20000, 기본 출력 경로
# 3. Module 1 분석을 수정 없이 적용
python analysis/m01_analysis_uproot.py --config config/analysis_m02_pythia.yaml \
       --figdir results/figures/m02 --tabdir results/tables/m02
```

**확인할 것**
- [ ] 10 이벤트: `sigma_gen` ≈ 78–80 mb 수준인지 (13 TeV INEL 단면적 자릿수;
      정확한 값은 PYTHIA 로그의 cross-section 표를 기록)
- [ ] 1,000 이벤트: INEL>0 비율, 이벤트당 최종상태 입자 수가 toy(~60)와 같은 자릿수인지
- [ ] 20,000 이벤트: Module 1 분석 스크립트가 **수정 없이** 끝까지 도는지.
      안 되면 그것이 Module 1 인터페이스의 결함 → 오답노트에 기록 후 수정
- [ ] `tests/test_m01_pipeline.py`를 PYTHIA 파일로도 돌려 볼 것 (closure test는
      toy 전용이므로 제외)

---

## 4. 데이터 비교 (다음 단계)

1. `data/raw/README.md`의 템플릿대로 ALICE pp 13 TeV $dN_{\rm ch}/d\eta$ (INEL>0)
   provenance JSON을 **먼저** 채우고 HEPData에서 표를 받는다.
2. `analysis/m02_compare_data.py` (작성 예정): PYTHIA INEL>0 선택 →
   $dN_{\rm ch}/d\eta$ → ALICE 표와 같은 binning으로 비교 + ratio panel.
3. 판정: $dN_{\rm ch}/d\eta|_{|\eta|<0.5}$가 측정의 ±15 % 이내. 벗어나면 tune 차이보다
   설정 오류(primary 정의, 이벤트 클래스, 하전입자 정의)를 먼저 의심.

---

## 5. 결정됨 — "하전입자"의 정의 → (a) `particle_charge ≠ 0` (2026-11-05)

Module 1 분석의 `CHARGED_HADRONS = {211, 321, 2212}`는 toy 생성기의 산물입니다.
PYTHIA 샘플에서 이대로 세면 ALICE의 "charged primary"보다 수 % 적게 셉니다
(e±, μ±, Σ±, Ξ±, Ω±가 빠짐). 선택지:

- (a) 분석을 `particle_charge != 0`로 바꾸고 toy 파일에도 `charge` 브랜치를 추가
  (Module 1 결과 재생성; 값은 동일해야 함 — 그것 자체가 검증)
- (b) Module 1 분석은 그대로 두고 Module 2 비교 스크립트만 charge를 사용

**(a)로 결정·반영 완료**: toy 생성기와 writer에 `charge` 브랜치를 추가하고 세 arm의
분석을 `particle_charge ≠ 0`으로 바꿨다. "같은 분석 코드가 수정 없이 돈다"는 목표의
참뜻은 "입력이 바뀌어도 정의가 흔들리지 않는다"이고, 그러려면 정의가 입력의
우연(toy에 π/K/p뿐)이 아니라 물리(전하)에 걸려 있어야 한다. Module 1 결과는
뮤온 두 개만큼 바뀌었고(`module01_root.md`), 그 과정에서 float32 히스토그램 버그
(오답노트 #12)를 찾았다. 따라서 `m02_run_pythia.py`가 쓰는 `particle_charge`(int8)는
Module 1 파일과 같은 브랜치·dtype이다.

---

## 6. 연습 — 단계별로 끄기

같은 seed로 아래를 하나씩 바꿔 1,000 이벤트씩 생성하고
$\langle N_{\rm ch}\rangle$, $dN/d\eta|_0$, $\langle p_{\rm T}\rangle$을 표로:

| 변경 | 예측(먼저 적을 것) |
|---|---|
| `PartonLevel:MPI = off` | 다중도 급감 (LHC 다중도의 주된 원천) |
| `PartonLevel:ISR = off` | |
| `PartonLevel:FSR = off` | |
| `ColourReconnection:reconnect = off` | $\langle p_{\rm T}\rangle(N_{\rm ch})$ 기울기 변화 |
| `HadronLevel:Decay = off` | tau0Max 안의 붕괴가 사라져 다중도 감소, 입자 종류 변화 |
| `SoftQCD:nonDiffractive = on`만 | SD/DD가 빠져 저다중도 꼬리 감소 |

실행 전 예측을 `docs/research_log.md`에 적고, 실행 후 대조합니다.
