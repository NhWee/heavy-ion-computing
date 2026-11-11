# Figure manifest

| figure | produced by | input | config |
|---|---|---|---|
| `m01_multiplicity.pdf/png` | `analysis/m01_analysis_uproot.py` | `data/generated/m01_toy_events.root` | `config/analysis.yaml` (seed 20260827) |
| `m01_dndeta.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |
| `m01_pt_spectrum.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |
| `m01_dimuon_mass.pdf/png` | `analysis/m01_analysis_uproot.py` | 〃 | 〃 |
| `m01_fit_validation.pdf/png` | `analysis/m01_fit_validation.py` | pseudo-experiments from `hicpy.toy_generator` (no input file) | `config/analysis.yaml` (`validation`, seed 20261107) |
| `m01_tsallis_systematics.pdf/png` | `analysis/m01_fit_systematics.py` | `data/generated/m01_toy_events.root` | `config/analysis.yaml` (`systematics`, bootstrap seed 20261108) |

fit 결과 전문(파라미터, 오차, 공분산·상관행렬, −2lnλ 또는 χ², ndf, fit range, multi-start
점검)은 `results/tables/m01_fit_results.json`. closure, 방법 검증, 체계 오차는 각각
`m01_closure_test.json`, `m01_fit_validation.json`, `m01_tsallis_systematics.json`.

`m01_root_overview.pdf`(C++ 매크로)는 ROOT의 PDF에 생성 시각이 들어가 바이트 재현이 안 되므로
커밋하지 않습니다(`.gitignore`). 정식 그림은 위 matplotlib 그림입니다.

모든 그림은 코드로만 생성되며 수작업 편집이 없습니다. 생성 시각·라이브러리
버전 메타데이터를 넣지 않으므로 같은 코드·설정·seed면 바이트까지 같습니다
(`hicpy.style.save_figure`, `docs/troubleshooting.md` #10).
재생성: `bash scripts/run_all.sh`
