# Figure manifest

| figure | produced by | input | config |
|---|---|---|---|
| `m01_multiplicity.pdf/png` | `analysis/m01_analysis_uproot.py` | `data/generated/m01_toy_events.root` | `config/analysis.yaml` (seed 20260827) |
| `m01_dndeta.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |
| `m01_pt_spectrum.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |
| `m01_dimuon_mass.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |

fit 결과 전문(파라미터, 오차, 공분산행렬, χ²/ndf, fit range)은
`results/tables/m01_fit_results.json`.

모든 그림은 코드로만 생성되며 수작업 편집이 없습니다.
재생성: `bash scripts/run_all.sh`
