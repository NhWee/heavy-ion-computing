# How-to: 실험 데이터(HEPData) 받아서 MC와 비교하기

누구든 이 저장소를 clone해서 아래를 그대로 치면 같은 결과가 나와야 합니다.
각 단계마다 **무엇을 하는지, 왜 하는지, 화면에 무엇이 보여야 하는지**를 적었습니다.
예제는 Module 2의 ALICE pp 13 TeV $dN_{\rm ch}/d\eta$ 표입니다.

필요한 것: 인터넷(`www.hepdata.net`), Python 환경(`hic` 또는
`pip install uproot awkward scipy pyyaml matplotlib`). PYTHIA는 4단계부터만 필요.

## 시작 전에 — 두 가지 전제 (여기서 가장 많이 막힙니다)

**(1) 모든 명령은 저장소 루트에서 칩니다.** 스크립트 경로가 `scripts/...`,
`data/raw/...`처럼 상대 경로이기 때문입니다. 홈 디렉토리에서 치면
`can't open file '.../scripts/m02_fetch_hepdata.py'`가 납니다.

```bash
cd "/mnt/c/Users/Administrator/Claude/Projects/HIC analysis/heavy-ion-computing"   # WSL2의 저장소 위치
git fetch origin && git checkout claude/vigilant-sagan-aatxhs && git pull
```
(저장소가 없으면 `git clone https://github.com/NhWee/heavy-ion-computing.git` 후 그 폴더로.)

**(2) 새 터미널에서는 micromamba 셸 초기화가 먼저입니다.** 안 하면
`'micromamba' is running as a subprocess and can't modify the parent shell`이 납니다.

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate hic
```
매번 치기 싫으면 한 번만 `micromamba shell init --shell bash --root-prefix=~/.local/share/mamba`
(다음 터미널부터 자동). 프롬프트 앞에 `(hic)`가 보이면 준비 끝입니다.
`(base)`는 conda 기본 환경이라 uproot/pyyaml이 없을 수 있습니다.

---

## 0. 원칙 — 왜 이렇게 복잡하게 하나

"HEPData에서 csv 받아서 플롯"은 10분이면 됩니다. 그런데 6개월 뒤에
"이 점들은 어느 논문의 어느 표였지? INEL이었나 INEL>0이었나? 계통 오차는 넣었나?"에
답할 수 없으면 그 그림은 못 씁니다. 그래서 **출처를 먼저 적고, 그 다음에 받습니다.**
데이터 파일(payload)은 커밋하지 않고(`.gitignore`), 출처 JSON만 커밋합니다 —
payload는 JSON만 있으면 언제든 다시 받을 수 있습니다.

---

## 1. 출처 JSON 확인 — `data/raw/alice_pp13TeV_dndeta.json`

```bash
cat data/raw/alice_pp13TeV_dndeta.json
```

**무엇**: 논문(arXiv 번호, 저널), HEPData 레코드 번호(`ins1395253`), 표 이름(`Table 1`),
관측량, 이벤트 클래스(`INEL>0`), 받을 파일 이름이 적혀 있습니다.
`"verified": false`는 "이 식별자들을 아직 HEPData 자체와 대조하지 않았다"는 뜻입니다.

**왜 verified가 false인가**: 이 JSON은 웹 검색 요약을 보고 썼습니다. 검색 요약은
1차 출처가 아닙니다. 다음 단계의 스크립트가 HEPData 레코드와 대조해 맞을 때만
`true`로 바꿉니다. `analysis/m02_compare_data.py`는 `false`인 JSON을 **거부**합니다.

새 데이터셋을 추가할 때는 `data/raw/README.md`의 템플릿을 복사해 식별자만 채웁니다.
숫자(데이터 값)는 절대 손으로 적지 않습니다.

---

## 2. 표 받기 — `scripts/m02_fetch_hepdata.py`

```bash
python scripts/m02_fetch_hepdata.py data/raw/alice_pp13TeV_dndeta.json
```

스크립트는 네 단계를 하고, 각 단계를 `=== step N/4` 줄로 알립니다.

| 단계 | 하는 일 | 왜 |
|---|---|---|
| 1/4 | `https://www.hepdata.net/record/ins1395253?format=json`을 받아 레코드의 arXiv·INSPIRE 번호를 JSON과 비교 | 레코드 번호를 잘못 적었으면 **다른 논문의 표**를 받게 됨. 불일치면 여기서 멈추고 아무것도 쓰지 않음 |
| 2/4 | 레코드의 표 목록을 출력하고 `Table 1`이 있는지 확인 | 표 이름이 바뀌었거나 번호를 잘못 알았을 때 목록을 보고 고를 수 있게 |
| 3/4 | `Table 1`을 YAML로 받아 `data/raw/alice_pp13TeV_dndeta_table1.yaml`에 저장, sha256 계산 | YAML이 HEPData의 정식 포맷(CSV는 export). sha256은 "그때 받은 파일과 지금 파일이 같다"의 증거 |
| 4/4 | 표의 헤더(qualifier)와 η 범위 출력 | **사람이 눈으로** 이벤트 클래스(INEL / INEL>0)와 η 범위를 확인하는 단계. 스크립트는 이것을 판단하지 않음 |

끝나면 JSON에 `verified: true`, `downloaded_utc`, `sha256`, 표 DOI가 채워집니다.

**화면에 보여야 하는 것** (예시 — 실제 문구는 HEPData가 주는 대로):
```
=== step 1/4: fetch the HEPData RECORD ins1395253 and check it is the paper named in ...
[get] https://www.hepdata.net/record/ins1395253?format=json
      title   : Pseudorapidity and transverse-momentum distributions of charged particles ...
      arXiv   : 1509.08734    inspire: 1395253    doi: 10.17182/hepdata.70847.v1
=== step 2/4: list the record's tables and pick 'Table 1'
[rec] N table(s):
      - Table 1: Average pseudorapidity density of charged particles ...
      ...
=== step 3/4: download 'Table 1' as YAML ...
[out] data/raw/alice_pp13TeV_dndeta_table1.yaml (NNNN bytes, sha256 ...)
=== step 4/4: print the table header -- CONFIRM BY EYE ...
      dependent: {'name': 'DN/DETARAP'}  qualifiers: [... 'INEL>0' ...]
      independent: {'name': 'ETARAP'}  range [-1.8, 1.8]  (N bins)
NOTE: 2 dependent variables in this table (e.g. INEL and INEL>0) ...
[out] data/raw/alice_pp13TeV_dndeta.json updated (verified = true)
```

**4/4에서 할 일**: 표에 INEL과 INEL>0 두 열이 있으면 어느 열이 INEL>0인지 보고
JSON에 `"dependent_index": 1`처럼 적습니다(0부터 셈). η 범위를 `selection.eta_range`에
적고 "(to be confirmed)"를 지웁니다. 이 한 줄이 사람이 책임지는 부분입니다.

**실패하면**:
- `ABORT -- the record does not match` → JSON의 arXiv/INSPIRE 번호가 틀림. HEPData
  페이지에서 번호를 확인하고 JSON을 고친 뒤 재실행.
- `ENOTFOUND` / `CONNECT ... 403` → 네트워크(프록시)가 `www.hepdata.net`을 막고 있음.
  클라우드 세션이면 환경 설정의 Network access에 추가 후 **새 세션**에서.
- `HTTP Error 403: Forbidden` (서버까지 갔는데 거부) → HEPData 앞의 Cloudflare가
  자동 다운로드를 막은 것. 스크립트는 브라우저형 헤더 → `curl` 순으로 재시도한다.
  그래도 안 되면 **브라우저로 받는다**: 레코드 페이지에서 Table 1 → "Download" →
  YAML을 `data/raw/alice_pp13TeV_dndeta_table1.yaml`로 저장한 뒤
  ```bash
  python scripts/m02_fetch_hepdata.py data/raw/alice_pp13TeV_dndeta.json --local
  ```
  `--local`은 다운로드 없이 sha256·헤더 출력·`verified: true`만 한다. 이때 식별자
  대조(1/4)는 사람이 페이지에서 한 것이므로 JSON의 `_status`에 "by hand"로 남는다.

---

## 3. 데이터 없이 MC만 먼저 보기 (선택)

```bash
python analysis/m02_compare_data.py --config config/analysis_m02_pythia.yaml --no-data \
       --mc data/generated/m01_toy_events.root
```
PYTHIA 파일이 아직 없어도 toy 파일로 스크립트가 도는지 볼 수 있습니다.
`[data] running MC only`가 뜨고 판정(verdict)은 null입니다.

---

## 4. PYTHIA 생성 (WSL2 `hic` 환경)

```bash
python scripts/m02_run_pythia.py --n-events 10      # 작동 확인
python scripts/m02_run_pythia.py                    # 카드의 20,000 이벤트
```
자세한 것은 `docs/modules/module02_pythia.md` §3.

---

## 5. 비교

```bash
python analysis/m02_compare_data.py --config config/analysis_m02_pythia.yaml
```

| 출력 줄 | 뜻 |
|---|---|
| `[data] ALICE dN_ch/deta INEL>0, N bins` | 출처 JSON을 통해 표를 읽음 |
| `[mc] ... N pass INEL>0 (xx %)` | 이벤트 클래스를 **입자에서 다시 계산**해 적용 (저장된 플래그를 믿지 않음) |
| `[mc] process mix: ND .. SD .. DD ..` | PYTHIA process code 분포 — SD/DD 비율이 이상하면 카드 의심 |
| `[mc] dN_ch/deta (|eta|<0.5, pT > 0) = ...` | pT cut 없음: 발표값은 pT→0 외삽값이므로 |
| `[cmp] data ... MC/data = 0.9xx -> WITHIN +-15 %` | 판정 |
| `[out] results/figures/m02/m02_dndeta_vs_data.pdf, results/tables/m02/m02_compare_summary.json` | 그림(위: 분포, 아래: ratio, 통계·계통 분리)과 JSON(양쪽 출처 포함) |

판정이 ±15 %를 벗어나면 tune이 아니라 **설정**을 먼저 의심합니다: 이벤트 클래스,
primary 정의(`tau0Max`), 하전입자 정의, pT cut 여부.

---

## 6. 커밋할 것 / 하지 않을 것

- 커밋: `data/raw/*.json`(verified = true가 된 것), `results/tables/m02/*.json`,
  `results/figures/m02/*`, 연구 노트
- 커밋 안 함: `data/raw/*.yaml|csv`(payload), `data/generated/*.root`
