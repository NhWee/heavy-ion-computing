# Module 1 연습문제

각 문제는 재현 → 파라미터 수정 → 예측 → 실행 → 해석 순서로 되어 있습니다.
**실행 전에 예측을 먼저 적어 두는 것**이 핵심입니다.

---

## E1-1 (재현) — 파일 구조 직접 확인하기

`data/generated/m01_toy_events.root`를 열어

```python
import uproot
f = uproot.open("data/generated/m01_toy_events.root")
print(f.classnames())
f["events"].show()
print(str(f["provenance"])[:400])
```

**답할 것**
1. `nparticle` branch가 없으면 무슨 일이 생기는가? (왜 counter branch가 필요한가)
2. `particle_px`의 typename이 `float[]`인 이유와, `double[]`로 바꾸면 파일 크기가
   어떻게 될지 예측한 뒤 `--format`을 바꾸지 말고 코드에서 `np.float64`로 저장해
   실제로 확인하시오.

---

## E1-2 (파라미터 수정 + 예측) — Tsallis 파라미터 흔들기

`config/analysis.yaml`에서 `tsallis_n`을 6.8 → 4.5로 바꿉니다.

**실행 전에 예측하시오**
1. $p_{\rm T}$ 스펙트럼의 고-$p_{\rm T}$ 꼬리는 더 단단해지는가(hard) 부드러워지는가?
2. $\langle p_{\rm T}\rangle$는 커지는가 작아지는가?
3. closure test의 $T$ 복원은 여전히 성공하는가?

그 다음 `bash scripts/run_all.sh`로 확인하고, 예측이 틀렸다면 **왜 틀렸는지**
`docs/research_log.md`에 적으시오.

---

## E1-3 (해석) — inclusive fit이 왜 틀리는가

closure test에서 inclusive Tsallis fit이 $T$를 13 MeV 높게 준다는 것을 확인했습니다.

1. $\pi$, $K$, $p$의 질량만 다르고 $(T, n)$은 같은데 왜 합치면 $T$가 올라가는가?
   힌트: 같은 $p_{\rm T}$에서 $m_{\rm T} = \sqrt{p_{\rm T}^2+m_0^2}$가 무거운 입자일수록
   크므로, 저-$p_{\rm T}$에서 무거운 입자가 상대적으로 억제된다 ($m_{\rm T}$ scaling).
2. `frac_K`와 `frac_p`를 0으로 두고 파이온만 생성하면 inclusive fit이 닫히는지
   확인하시오. 닫힌다면 §1의 설명이 맞다는 증거가 된다.
3. 실험 논문에서 "inclusive charged-particle spectrum"에 대한 Tsallis fit 파라미터를
   인용할 때 어떤 주의가 필요한가?

---

## E1-4 (설계) — Z 질량의 편차: 요동인가, 편향인가

Session 1의 χ² fit은 $m_Z = 91.06 \pm 0.09$ GeV로 생성기 입력 91.1876보다 1.4σ 낮았습니다.
이를 "pT 컷의 acceptance sculpting"으로 설명했다가 철회했습니다(오답노트 #15).
이 문제는 그 판단을 직접 다시 해 보는 것입니다.

1. `results/tables/m01_closure_test.json`의 truth 수준 Z(컷 없음)와
   `m01_fit_results.json`의 `z_fit`(컷 있음)을 비교하시오. 편차는 컷 이전부터 있었는가?
2. **실행 전에 예측하시오.** 뮤온 $p_{\rm T}$ 컷을 20 → 30 GeV/c로 올리면 acceptance의
   질량 기울기(`m01_fit_validation.json`의 `relative_slope_per_GeV`)는 커지는가,
   작아지는가? 힌트: Z 정지계에서 $p^* = m/2$이고, 두 뮤온이 모두 컷을 넘으려면
   $\sin\theta^* > p_{\rm T}^{\rm cut}/p^*$이다. `config/analysis.yaml`의
   `muon_cuts.pt_min`을 바꾸고 `python analysis/m01_fit_validation.py`로 확인하시오.
   피크 이동(`peak_shift_from_cuts`)은 통계 오차와 비교해 언제 무시할 수 없게 되는가?
3. 같은 파일에서 Neyman χ²와 Poisson likelihood의 Γ pull을 비교하시오.
   `z_fraction`을 0.02 → 0.2로 올리면(후보 약 10배) 차이가 어떻게 될지 예측한 뒤 확인하시오.
   χ²가 likelihood에 가까워지는 조건을 한 문장으로 쓰시오.
4. 실험은 acceptance를 어떻게 다루는가? (검색어: acceptance correction, template fit,
   unfolding) 한 문단으로 정리하시오.

---

## E1-5 (도전) — ROOT arm 대조

WSL2 환경이 준비되면:

```bash
root -l -b -q 'analysis/m01_analysis_root.C("data/generated/m01_toy_events.root")'
python analysis/m01_analysis_pyroot.py
```

**같은 입력 파일**에서 나온 `<N_ch>`, `dN/dη|₀`, Z 질량이 uproot arm과 일치하는지
확인하시오. 불일치가 있다면 그 원인은 거의 항상 (a) bin edge 정의, (b) 정규화
(`Scale(1/N, "width")` vs 직접 나누기), (c) 컷 순서 중 하나입니다. 어느 것이었는지
`docs/troubleshooting.md`에 기록하시오.
