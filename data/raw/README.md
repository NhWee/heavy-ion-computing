# data/raw/ — 실험 데이터 provenance

실험 데이터 payload(csv/yaml/root)는 커밋하지 않습니다(`.gitignore`).
대신 데이터셋마다 `<dataset>.json`에 **출처를 먼저** 기록하고, 그 다음에 받습니다.
숫자를 지어내지 않습니다: 아래 템플릿의 값은 HEPData에서 받은 뒤에만 채웁니다.

## 템플릿 (`alice_pp13TeV_dndeta.json`, Module 2)

```json
{
  "experiment": "ALICE",
  "title": "<논문 제목 그대로>",
  "arxiv": "<arXiv ID -- 확인 후 기입>",
  "doi_paper": "<10.xxxx/...>",
  "hepdata": {"record": "<ins...>", "table": "<Table N>", "doi": "<10.17182/hepdata....>"},
  "system": "pp", "sqrt_s_TeV": 13.0,
  "observable": "dN_ch/deta",
  "event_class": "INEL>0",
  "selection": {"eta_range": "<from the table header>", "pt_cut": "<if any>",
                "primary_definition": "<from the paper>"},
  "uncertainties": {"stat": "<column>", "syst": "<column>", "correlated": "<yes/no, how>"},
  "downloaded_utc": "<YYYY-MM-DDThh:mm:ssZ>",
  "file": "alice_pp13TeV_dndeta.csv",
  "sha256": "<of the downloaded file>"
}
```

## 채우는 절차 (`alice_pp13TeV_dndeta.json`)

1. JSON에 **식별자만** 적는다 (`verified: false`). 지금 커밋된 상태가 이것이다:
   ALICE, "Pseudorapidity and transverse-momentum distributions of charged
   particles in proton–proton collisions at √s = 13 TeV", arXiv:1509.08734,
   PLB 753 (2016) 319, HEPData `ins1395253` Table 1 — **웹 검색 요약에서 얻은
   값이라 미검증**. 이 세션의 클라우드 환경은 hepdata.net·inspirehep.net·arxiv.org에
   접근할 수 없었다.
2. 네트워크가 되는 곳(WSL2)에서:
   ```bash
   python scripts/m02_fetch_hepdata.py data/raw/alice_pp13TeV_dndeta.json
   ```
   레코드 메타데이터와 JSON의 arXiv/INSPIRE 번호가 **일치할 때만** 표(YAML)를 받고
   sha256·시각·DOI를 채워 `verified: true`로 바꾼다. 불일치면 아무것도 쓰지 않는다.
3. 스크립트가 출력하는 표 헤더(qualifier)로 이벤트 클래스와 η 범위를 **눈으로 확인**하고
   `selection`의 "(to be confirmed)" 문구를 지운다. 표에 INEL과 INEL>0이 같이 있으면
   `"dependent_index"`로 열을 고른다.
4. `analysis/m02_compare_data.py`는 `verified: false`인 JSON을 거부한다.
