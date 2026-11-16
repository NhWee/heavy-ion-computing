# generators/

이벤트 생성기 설정 카드를 둡니다. **설정을 숨기지 않는다**가 원칙입니다.

PYTHIA 카드(`*.cmnd`)에는 최소한 다음이 명시되어야 합니다:
빔 입자, 충돌 에너지, 켠 process, tune, PDF set, 운동학 cut,
난수 시드, 이벤트 수. 기본값이라도 물리에 영향을 주는 것은 **적어 둡니다** —
카드만 읽고 무엇을 돌렸는지 알 수 있어야 합니다.

| 카드 | 용도 | 실행 |
|---|---|---|
| `pythia_mb_13TeV.cmnd` | Module 2: pp 13 TeV, `SoftQCD:inelastic` (ND+SD+DD+CD), Monash 2013 | `python scripts/m02_run_pythia.py --card generators/pythia_mb_13TeV.cmnd --n-events 10` |

이벤트 클래스(INEL, INEL>0, NSD)는 **카드가 아니라 분석에서** 고릅니다.
`process_code`와 `n_ch_eta1` 브랜치가 그 선택을 가능하게 합니다.
