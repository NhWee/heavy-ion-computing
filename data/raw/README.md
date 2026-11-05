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

후보 논문(**미검증** — HEPData/INSPIRE에서 확인하고 채울 것):
ALICE, "Pseudorapidity and transverse-momentum distributions of charged
particles in proton–proton collisions at √s = 13 TeV", arXiv:1509.08734.
