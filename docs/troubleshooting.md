# 오답노트 — 실제로 만난 오류와 그 원인

각 항목: **증상 → 원인 → 해결 → 배운 점**.
"명령이 exit 0이었으니 됐다"를 신뢰하지 않기 위한 기록입니다.

---

## #1 uproot 5.7이 TTree 대신 RNTuple을 써버림

**증상**
```python
with uproot.recreate(out) as f:
    f["events"] = arrays          # awkward 배열 dict
```
로 만든 파일을 다시 열어 `f.classnames()`를 찍으면
```
{'events;1': 'ROOT::RNTuple', 'provenance;1': 'TObjString'}
```
TTree가 아니라 **RNTuple**이 나온다.

**원인**
uproot 5.7부터 `directory[name] = dict-of-arrays` 형태의 대입은 기본적으로
RNTuple로 기록된다. RNTuple은 TTree의 차세대 대체 컨테이너지만 ROOT 6.34 이상에서만
읽을 수 있고, 대부분의 실험 프레임워크·기존 매크로는 아직 TTree를 전제한다.

**해결**
`mktree()`로 명시적으로 TTree를 만들고 `extend()`로 채운다.
```python
tree = f.mktree("events", branch_types, counter_name=lambda n: "n" + n)
tree.extend(arrays)
```
`scripts/m01_make_toy_events.py --format rntuple`로 RNTuple을 선택할 수도 있게 남겨둠.

**배운 점**
파일을 쓰고 나면 **반드시 `classnames()`로 실제 컨테이너 타입을 확인**한다.
포맷 호환성은 "읽는 쪽의 ROOT 버전"이 결정한다.

---

## #2 여러 jagged branch가 counter branch를 공유하지 못함

**증상**
```python
f.mktree("events", types, counter_name=lambda _: "nparticle")
# struct.error: required argument is not an integer
```

**원인**
uproot 5.7.6에서 서로 다른 가변길이 branch들이 **같은 이름의 counter branch**를
공유하도록 강제하면 내부 직렬화가 깨진다. 기본값(`lambda c: "n"+c`)을 쓰면
`npx, npy, npz, nE, npdg` 처럼 counter가 5개나 생겨 파일이 지저분해진다.

**해결**
입자 변수들을 하나의 **record branch**로 묶는다.
```python
particle = ak.zip({"px":…, "py":…, "pz":…, "E":…, "pdg":…})
arrays = {"particle": particle, "has_z": …}
```
결과 branch 구조:
```
nparticle, particle_px, particle_py, particle_pz, particle_E, particle_pdg, has_z
```
이는 CMS NanoAOD의 `nMuon / Muon_pt / Muon_eta …` 와 정확히 같은 관례다.

**배운 점**
"보기 좋게" 만들려던 우회가 오히려 실제 실험 파일 구조와 일치했다.
HEP 파일 레이아웃은 대체로 이유가 있어서 그 모양이다.

---

## #3 numpy의 `clip`/`isin`은 awkward 배열에서 터진다

**증상**
```python
eta = np.arctanh(np.clip(pz / p, -1+1e-12, 1-1e-12))
# ValueError: cannot convert to RegularArray because subarray lengths are not regular
```

**원인**
awkward 배열은 numpy **ufunc**(`sqrt`, `hypot`, `arctanh`, `+`, `>` …)에 대해서만
element-wise 브로드캐스트를 지원한다. `np.clip`, `np.isin`은 ufunc이 아니라
`__array_function__` 디스패치이고, awkward는 이를 처리하려고 배열을 직사각형으로
만들려다 실패한다 — 이벤트마다 입자 수가 다르니 당연히 불가능하다.

**해결**
- `np.clip(x,-1,1)` → `x * (1 - 1e-12)` (|pz|/|p| ≤ 1 이 이미 보장됨)
- `np.isin(apdg, [211,321,2212])` → `(apdg==211)|(apdg==321)|(apdg==2212)`
- 꼭 필요하면 `ak.where`, 또는 `ak.flatten` 후 numpy로.

**배운 점**
jagged 데이터에서는 **ufunc만 쓴다**가 안전한 기본 규칙.

---

## #4 numpythia(PyPI PYTHIA 바인딩) 빌드 실패

**증상**
```
Error compiling Cython file: numpythia/src/_libnumpythia.pyx:25:0:
'cpython/cobject.pxd' not found
error: failed-wheel-build-for-install
```

**원인**
`numpythia` 1.3.1은 Python 2 시절의 `PyCObject` API에 의존한다. 이 API는 Python 3.2
에서 제거되었고 Cython도 해당 `.pxd`를 더 이상 제공하지 않는다. 패키지가 사실상
유지보수되지 않는 상태.

**해결**
클라우드 컨테이너에서는 PYTHIA를 포기한다. Module 2부터는 WSL2에서
`micromamba install pythia8` (conda-forge 공식 빌드, Python 바인딩 포함)을 쓴다.

**배운 점**
"PyPI에 이름이 있다"는 것이 "쓸 수 있다"를 뜻하지 않는다. HEP 소프트웨어는
**conda-forge가 사실상의 표준 배포 경로**다.

---

## #5 mounted 폴더에서 git이 `index.lock`을 지우지 못함

**증상**
```
warning: unable to unlink '.../.git/index.lock': Operation not permitted
```
이후 git 명령이 stale lock 때문에 막힘.

**원인**
Cowork 로컬 VM에서 Windows 폴더 마운트는 기본적으로 **삭제 권한이 없다**.
git은 거의 모든 명령에서 lock 파일을 만들고 지운다.

**해결**
세션에서 해당 폴더에 대한 삭제 권한을 요청해 승인받음. 이후 정상 동작.

**배운 점**
버전 관리를 마운트된 네트워크/가상 파일시스템 위에서 할 때는 삭제·파일잠금 권한을
먼저 확인한다.

---

## #6 float32로 저장한 4-momentum에서 질량을 되계산하면 안 된다

**증상**
저장된 `(px,py,pz,E)`에서 `m = sqrt(E²-p²)`를 다시 계산하면 파이온이
0.099–0.150 GeV로 흩어지고, 고-pT 입자에서는 음수 m²가 나와 0으로 clip된다.

**원인**
float32는 유효숫자 ~7자리. E ≈ |p| 인 상대론적 입자에서 `E²-p²`는 거의 완전한
상쇄(catastrophic cancellation)이고, 상대오차가 m²에 증폭된다.

**해결**
질량은 **되계산하지 말고 PDG ID로 조회**한다 (`particle` 패키지 또는
`hicpy.toy_generator.MASS`). 불변질량처럼 여러 입자를 합할 때는 합이 충분히
무거워(m ≫ 개별 오차) 문제가 없다 — Z→μμ의 91 GeV 피크는 float32로도 멀쩡하다.

**배운 점**
실험 파일이 pT/eta/phi/mass를 따로 저장하는 이유가 여기에 있다.

---

## #7 두 arm이 어긋나는 진짜 원인: 1/pT Jacobian을 어디서 적용하는가

**증상**
uproot arm과 ROOT arm이 같은 파일에서 불변 pT 수율을 몇 %씩 다르게 냈다
(실행 전 코드 검토에서 발견).

**원인**
불변 수율 $\frac{1}{2\pi p_{\rm T}}\frac{d^2N}{dp_{\rm T}d\eta}$ 에서 $1/p_{\rm T}$를
- Python 초판: **bin 중심** $p_{\rm T}$로 나눔 (히스토그램을 채운 뒤 일괄 처리)
- ROOT 매크로: `Fill(pt, 1/(2π·pt·Δη))` 로 **입자마다** 자기 $p_{\rm T}$ 사용

로그 bin은 폭이 넓어서(고-pT에서 한 bin이 수 GeV) 두 방식이 수 % 차이가 난다.

**해결**
둘 다 **입자별 가중치**로 통일. 물리적으로도 이쪽이 옳다 — bin 중심 근사는
bin 안에서 스펙트럼이 평평하다고 가정하는데, 급격히 떨어지는 스펙트럼에서는 거짓이다.
가중 히스토그램이므로 오차도 $\sqrt{N}$이 아니라 $\sqrt{\sum w_i^2}$ (= ROOT의 `Sumw2`).
`hist_weighted()` 헬퍼로 분리.

**부수 효과**
파이온 closure fit의 χ²/ndf가 0.91 → 0.50으로 개선. 편향이 있었다는 방증.

**배운 점**
같은 관측량의 두 구현이 어긋날 때 원인은 거의 항상 순서대로:
**bin edge → 정규화 → 컷 순서 → float 정밀도**. 이번엔 정규화였다.

---

## #8 "η = 0에서의 값"은 정의가 모호하다

**증상**
`dN/dη|_{η≈0}`을 uproot arm은 `argmin(|bin center|)`, ROOT arm은 `FindBin(0.0)`으로
잡았는데, bin edge가 0에 정확히 놓이면 **서로 다른 bin**을 고른다.

**해결**
"η=0의 bin 값" 대신 **|η| < 0.5 구간의 평균 밀도**를 세어서 정의
(`n_central / (N_ev · Δη)`, 히스토그램을 아예 거치지 않음). 두 arm 모두 동일.

**배운 점**
비교 가능한 숫자는 **binning에 의존하지 않게** 정의해야 한다.
실험 논문이 "dN/dη at midrapidity"를 쓸 때 항상 |η| 범위를 명시하는 이유.

---

## #9 WSL2 첫 실행 시 systemd 경고 / conda 환경의 SKIP 2건

**증상 1**
```
wsl: Failed to start the systemd user session for 'uno'. See journalctl for more details.
```
**원인/해결**: WSL2에서 systemd user session이 뜨지 않는 것은 흔한 경고이며
셸·컴파일·conda 사용에는 영향이 없다. 실제로 그 뒤 모든 설치와 검증이 정상 동작했다.
systemd 서비스(예: docker daemon)를 WSL 안에서 돌릴 때만 문제가 된다.

**증상 2**
`verify_environment.py`에서 `pyhepmc`와 `lhapdf`가 SKIP.

**원인/해결**
- `pyhepmc`: `environment.yml`에 `hepmc3`만 넣었는데 이건 **C++ 라이브러리와 헤더**만
  준다. Python 바인딩은 별도 패키지 `pyhepmc`. → environment.yml에 추가함.
- `lhapdf`: Module 6에서 설치 예정이라 의도적으로 SKIP. PDF set은 수 GB가 될 수
  있어 필요할 때 받는다.

**배운 점**
conda-forge에서 C++ 라이브러리 이름과 Python 바인딩 이름이 다른 경우가 흔하다
(`hepmc3` vs `pyhepmc`, `fastjet` vs `fastjet-cxx`). 검증 스크립트가 없었으면
Module 3에 가서야 알았을 것이다.

---

## #10 같은 seed, 같은 숫자인데 그림 파일이 매번 달라진다

**증상**
`run_all.sh`를 다시 돌리면 `m01_fit_results.json`의 숫자는 그대로인데
`git status`에 `results/figures/*.pdf|png` 8개가 전부 수정으로 뜬다.
`cmp`로 보면 바이트가 다르다.

**원인**
matplotlib이 PDF에는 `/CreationDate`, `/ModDate`, `/Producer`를, PNG에는
`Software` 텍스트 청크를 기본으로 심는다. 생성 **시각**과 matplotlib **버전**이
파일 안에 들어가므로, 같은 데이터로 그려도 같은 파일이 되지 않는다.

**해결**
`hicpy.style.save_figure()`가 `metadata={"CreationDate": None, ...}`로 그 필드를
비운다. 모든 분석 스크립트는 `fig.savefig()` 대신 이 함수를 쓴다.
`tests/test_m01_pipeline.py::test_figures_are_byte_reproducible`가 두 번 그려
바이트 비교로 이를 지킨다.

**배운 점**
"재현 가능"의 정의를 숫자에서 **파일**까지 넓히면 `git diff`가 곧 회귀 검출기가
된다. 그림이 바뀌었다는 것은 코드·설정·seed 중 하나가 바뀌었다는 뜻이어야 한다.

---

## #11 arm마다 cut 집합이 달랐다 (`pt < pt_max`)

**증상**
실행 전 코드 대조에서 발견. uproot arm의 $N_{\rm ch}$는
`pt_min < pT < pt_max`를, PyROOT arm과 C++ 매크로는 `pT > pt_min`만 적용하고
있었다. 반대로 uproot arm의 $p_{\rm T}$ 스펙트럼은 pT cut 없이 |η|만 걸고 있었다.

**원인**
cut을 세 곳에 따로 썼기 때문. `config/analysis.yaml`에는 `pt_max: 50.0`이
있지만 두 ROOT arm은 그 키를 읽지 않았다. 생성 범위가 60 GeV/c라
50–60 GeV/c 사이 입자가 실제로 존재하므로 숫자가 달라졌을 것이다.

**해결**
세 arm 모두 track 선택을 "하전 하드론 ∧ |η| < 0.8 ∧ pt_min < pT < pt_max"로
통일하고 $N_{\rm ch}$와 $p_{\rm T}$ 스펙트럼이 **같은 입자 집합**을 쓰게 했다.
$dN/d\eta$는 세 arm 모두 pT 상한 없이 `pT > pt_min`만 적용한다.
`m01_compare_arms.py`가 이제 네 히스토그램을 bin별로 비교하므로, 이런 차이는
다음부터 자동으로 잡힌다.

**배운 점**
"같은 분석"을 두 언어로 쓰면 cut 하나가 빠지는 건 시간문제다. 스칼라 요약 몇 개가
아니라 **히스토그램 전체**를 대조해야 하고, C++ 매크로처럼 설정 파일을 못 읽는
코드는 기본 인자 옆에 "config와 같아야 함"을 적는 것으로는 부족하다 — 대조
스크립트가 cut 값 자체를 비교해야 한다.

---

## #12 float32 가중치로 히스토그램을 채우면 희소한 bin의 오차가 0이 된다

**증상**
불변 $p_{\rm T}$ 수율의 6 GeV/c 위 bin들이 내용은 있는데 오차가 **정확히 0**.
`do_fit`은 오차 0인 점을 제외하므로 Tsallis fit 범위가 "우연히" 6 GeV/c에서
끝났고, 아무도 눈치채지 못했다. 내용 자체도 1/128의 배수로 양자화되어 있었다.

**원인**
ROOT 파일의 4-운동량은 float32로 저장된다(크기 때문; 맞는 선택). 그런데
`pt = np.hypot(px, py)`, `w = 1/(2π pt Δη)`가 전부 **float32로 전파**됐다.
`np.histogram(values, weights=w)`는 가중치를 정렬 후 **누적합**해서 bin 경계의
차이로 bin 내용을 구한다. 누적합이 float32로 전체 합(큰 수)에 도달하면 꼬리 bin의
차이는 그 수의 ulp 아래로 떨어진다. $w^2$는 더 작아서 꼬리에서 정확히 0.

**해결**
읽자마자 float64로 올린다 (`ak.values_astype(arr, np.float64)`), 그리고
`hist_weighted` / `hist_with_errors`는 인자의 dtype을 믿지 않고 float64로 캐스팅한다.
PyROOT arm의 `calcPt/calcEta`도 `RVec<double>`을 반환하도록 바꿨다 (C++ 매크로는
처음부터 double이었다). fit 범위는 이제 `config/analysis.yaml`의
`fits.tsallis_pt_range = [0.15, 6.0]`에 **명시**한다 — 범위가 "빈 bin이 나올 때까지"로
암묵적으로 정해지면 이런 버그가 숨는다.

**배운 점**
저장 정밀도와 연산 정밀도는 다른 문제다. float32 저장은 괜찮지만 그 값으로 합·누적합·
차이를 만들면 안 된다. 그리고 "fit이 잘 됐다"(χ²/ndf ≈ 1)는 fit에 **들어간 점들**에
대한 말이다 — 어떤 점이 왜 빠졌는지 보지 않으면 좋은 χ²는 아무것도 보증하지
않는다. 이 버그는 하전입자 정의를 바꾸다가(#11 다음 단계) 발견했다: 뮤온이 꼬리에
들어왔는데도 fit 범위가 늘어나지 않는 것이 이상했다.

---

## #13 `can't open file 'scripts/...'` / `micromamba ... can't modify the parent shell`

**증상 1**: `python scripts/m02_fetch_hepdata.py ...` → `No such file or directory`.
**원인/해결**: 홈 디렉토리에서 실행함. 이 저장소의 모든 명령은 **저장소 루트**에서
친다(상대 경로). WSL2의 저장소는 `/mnt/c/Users/Administrator/Claude/Projects/HIC analysis/heavy-ion-computing`.

**증상 2**: `micromamba activate hic` → `'micromamba' is running as a subprocess and
can't modify the parent shell`.
**원인/해결**: 새 터미널에서 셸 hook이 안 읽힌 상태. `eval "$(micromamba shell hook
--shell bash)"` 후 activate. 영구히 하려면 `micromamba shell init --shell bash
--root-prefix=~/.local/share/mamba` 한 번.

**배운 점**: 가이드에 "어디서", "어떤 환경에서"가 빠지면 첫 줄에서 막힌다.
`docs/howto/*.md`는 이제 "시작 전에" 절로 시작한다.

**증상 3**: `micromamba shell init --shell bash --root-prefix=~/.local/share/mamba` 실행 후
새 터미널에서 `hic` 환경이 안 보임.
**원인/해결**: root prefix는 **환경이 실제로 있는 곳**이어야 한다. `hic`는
`~/micromamba/envs/hic`에 있으므로 `micromamba shell init --shell bash --root-prefix=~/micromamba`.
(micromamba의 안내 메시지가 보여주는 기본값을 그대로 쓰면 안 된다 —
`micromamba env list`로 위치를 먼저 확인.)

---

## #14 PYTHIA 카드의 `Next:numberShowEvent = 0`이 무시되고 첫 이벤트가 통째로 출력됨

**증상**: WSL2에서 `m02_run_pythia.py --n-events 10` 실행 시 첫 이벤트의 Info / hard
process / complete event listing이 stdout에 찍힘. 시작 시 "changes only" 설정표에도
`Next:*` 항목이 없음 (0으로 바꿨으면 기본값 1과 달라 표에 나와야 함).

**원인**: `Print:quiet`은 **상위 스위치**다. `on`이든 `off`든 설정하는 순간 모든
`Init:show*`, `Next:show*`를 자기 기본값으로 되돌린다. 카드에서 `Print:quiet = off`가
`Next:*` 뒤에 있어서 방금 0으로 바꾼 값을 1로 되돌렸다. 한 줄씩 `readString`으로 넣어
값을 추적해 확인함 (PYTHIA 8.318: 75번 줄에서 0, 77번 줄 `Print:quiet` 뒤에 다시 1).

**해결**: `Print:quiet`을 `Next:*` **앞**으로 옮김. 카드에 이유를 주석으로 남김.

**배운 점**: PYTHIA 설정은 순서가 의미를 갖는 경우가 있다 (`Print:quiet`, `Tune:pp`처럼
다른 설정을 한꺼번에 바꾸는 스위치). 시작 시 출력되는 "changes only" 표가 카드의
의도와 맞는지 보는 것이 가장 싼 검사다.

---

## #15 PyROOT arm의 Tsallis fit: T = 56 MeV, χ²/ndf ≈ 10⁴ (히스토그램은 완전히 일치)

**증상**: WSL2에서 처음 실행한 3-arm 대조. PyROOT arm의 스칼라 5개와 히스토그램 4종
220 bins가 uproot arm과 일치(최대 상대 차이 1e-10; dN/dη 한 bin만 경계 입자 1개)하는데
Tsallis fit만 T = 0.0559 GeV (uproot 0.1482), n = 10.45, χ²/ndf = 10876.

**원인 (가설, 재실행으로 확인 예정)**: `SetParLimits(0, 0, 1e30)`. MINUIT은 범위가 있는
파라미터를 내부적으로 sin 변환해 다루는데, 범위가 10³⁰이면 시작값(첫 bin × 2)이 변환의
평평한 끝에 놓여 기울기가 ~0이 된다. C가 시작값(참값의 약 2배)에 얼어붙고 T와 n이 그
어긋남을 떠안음. scipy의 bounds는 이렇게 스텝을 왜곡하지 않아서 uproot arm은 정상이었다.
같은 데이터·같은 함수·같은 시작값에서 결과가 갈렸다는 것이 근거.

**해결**: C에는 범위를 두지 않음 (PyROOT, C++ 매크로 둘 다). T, n의 물리적 범위는 유지.

**배운 점**: 히스토그램 대조를 fit 대조와 분리해 둔 덕에 "계산은 같고 fit만 다르다"가 한눈에
보였다. MINUIT에 넓은 범위를 "안전장치"로 주면 오히려 fit을 망가뜨릴 수 있다 — 범위는
물리적으로 필요한 파라미터에만.

(같은 실행에서 C++ 매크로의 `m01_root_summary.json`은 첫 글자부터 JSON이 아니었음 —
원인 조사 중. `m01_compare_arms.py`가 그 지점에서 traceback으로 죽던 것은 FAIL로 보고하도록 수정.)

**#15 확인**: C 범위 제한을 뺀 뒤 C++ 매크로 fit = T 148.193 ± 0.66 MeV, n 7.4786, χ²/ndf 1.0847
→ uproot arm과 일치. 진단이 맞았음.

---

## #16 C++ 매크로가 끝에서 segfault, `m01_root_summary.json`이 0 바이트

**증상**: `[out] histograms -> ...` 직후 `*** Break *** segmentation violation`. 요약 JSON은
생성은 됐지만 비어 있음 → `m01_compare_arms.py`에서 "Expecting value: line 1 column 1".

**원인**: ROOT 객체 소유권. 입력 `TFile`을 연 **뒤에** 만든 `TH1`은 그 파일(gDirectory)에
소속된다. `fin->Close()`가 소속 히스토그램을 지우고, 그 뒤의 JSON 작성 코드가 해제된
메모리를 읽다가 죽었다. `std::ofstream`의 버퍼는 flush 전에 프로세스가 죽어 0 바이트.

**해결**: 히스토그램 생성 전에 `TH1::AddDirectory(kFALSE)` — 히스토그램을 파일에 소속시키지
않는다(쓰기는 `fout`이 열린 상태에서 `Write()`로 그대로 됨). 덤으로 `Fit` 옵션에 `0`을
넣어 batch에서 불필요한 기본 캔버스 `c1`이 생기지 않게 함.

**배운 점**: PyROOT/C++에서 "객체를 만든 순간 어느 디렉터리에 붙는가"는 Python 쪽
(uproot)에는 없는 함정이다. 요약 파일이 비어 있으면 쓰기 도중에 죽었다는 신호다.

**#16 확인**: 수정 후 매크로 정상 종료, `m01_root_summary.json` 생성. PyROOT arm은 스칼라·히스토그램·
fit 모두 uproot와 일치 (fit T 0.148193 vs 0.148193).

---

## #17 C++ 매크로의 dN/dη만 19 bins에서 어긋남 (최대 0.44 %)

**증상**: 3-arm 대조에서 C++ 매크로가 스칼라, P(N_ch), pT, m(μμ), Tsallis fit은 모두 일치하는데
dN/dη 히스토그램만 60 bins 중 19개가 다름(각각 0.5σ 이내). `compare_arms`가 "경계 입자
하나로 설명하기엔 너무 많다"로 FAIL 처리. 중앙 |η| < 0.5 스칼라는 정확히 일치.

**원인**: `px[i] * px[i]` — 분기 버퍼가 `Float_t`라 곱 자체가 float32로 계산된 뒤 double로
승격됨. 전방 영역에서 η = ½ ln((p+pz)/(p−pz))의 p − pz가 상쇄되어 float32 오차가 증폭 →
입자가 옆 bin으로 이동. 중앙에서는 상쇄가 없어서 스칼라는 맞았다. 오답노트 #12(numpy
float32)와 같은 계열의 실수가 C++ 쪽에 남아 있던 것.

**해결**: 곱하기 전에 `Double_t`로 승격 (`const Double_t dpx = px[i]` …). PyROOT arm은 #12 때
이미 double로 바꿨기 때문에 이 문제가 없었다.

**배운 점**: "같은 정밀도 규칙"은 언어마다 따로 지켜야 한다. bin별 대조가 없었다면 중앙
스칼라만 보고 일치한다고 결론냈을 것이다.

---

## #18 직접 만든 HepMC3 파일이 다시 읽히지 않음 / pyhepmc segfault

**증상 1**: `pyhepmc`로 쓴 파일을 읽으면 `ReaderAscii: too few vertices were parsed
(152 vs 159 expected)`, 이벤트 0개. 그런데 검사 스크립트는 "GLOBAL CHECKS PASSED"를 출력.
**증상 2**: `f.read().particles[0]`의 속성을 읽는 순간 segfault.

**원인 1**: PYTHIA 레코드 번호는 **시간 순서가 아니다**. ISR은 하드 산란에서 빔 쪽으로
**역방향 진화**하므로, 하드 과정 파톤(3번)의 엄마가 7번에 있다. HepMC3 ASCII reader는
입자가 생성 꼭짓점의 입력 입자들보다 뒤에 나와야 하는데, PYTHIA 순서대로 쓰면 앞선
참조(forward reference)가 생겨 파싱이 실패한다. (그 전 시도 — 엄마 *집합*으로 꼭짓점을
묶는 방식 — 은 HepMC의 "입자당 끝 꼭짓점 하나" 규칙과도 충돌했다.)
**원인 2**: pyhepmc에서 입자·꼭짓점은 GenEvent가 소유한다. 임시 이벤트가 사라지면
입자는 해제된 메모리를 가리킨다 (#16의 ROOT TH1 소유권과 같은 계열).
**원인 3 (검사의 결함)**: 이벤트 0개면 모든 "이벤트별" 검사가 공허하게 참이 된다.

**해결**: `hicpy.eventrecord.topological_order`로 엄마가 항상 먼저 오게 정렬해서 쓰고,
꼭짓점 규칙은 Pythia8ToHepMC3와 같은 "첫 엄마의 끝 꼭짓점"으로. 읽을 때는 이벤트 객체를
변수에 붙잡아 둔다. 검사 스크립트는 이벤트 0개면 FAIL.

**배운 점**: 레코드의 "번호"와 물리적 "순서"는 다르다 — 그게 backward evolution이다.
그리고 "모든 X에 대해 참"인 검사는 X가 비어 있으면 아무것도 말하지 않는다.

---

## #19 Module 1 분석을 PYTHIA 파일에 돌리면 Z fit이 빈 히스토그램에 "수렴"

**증상**: `m01_analysis_uproot.py --config config/analysis_m02_pythia.yaml`이 수정 없이 끝까지
돌지만 `Z peak m = 91.200 +- 1.3e10 GeV, chi2/ndf = 0.00` — converged = true. 그림 범례는
"toy MC".
**원인**: 최소 바이어스 PYTHIA에는 Z가 없다 (후보 0). curve_fit은 0인 bin들에 시작값
그대로 "수렴"하고, 공분산이 발산한다. 범례·스탬프 문구는 toy 전제로 하드코딩돼 있었다.
**해결**: 후보 20개 미만이면 fit을 시도하지 않고 사유를 기록 (`converged: false`,
`error: "not attempted: 0 candidates < 20"`). 범례와 스탬프는 config의 생성기 이름에서.
**배운 점**: "코드가 수정 없이 돈다"는 기준은 통과했지만, "돌았다"와 "의미 있는 결과를
냈다"는 다르다. 실패해야 할 fit은 실패를 보고해야 한다.

---

## #20 FastJet Python: "ClusterSequence has gone out of scope" / 고스트 seed를 못 넘김

**증상 1**: `cluster(...)[0][0].area()` → `FastJetError: you requested information about the
internal structure of a jet, but its associated ClusterSequence has gone out of scope`.
**원인**: 제트의 면적·구성 입자 정보는 ClusterSequence가 들고 있다. 결과 튜플에서 제트만
꺼내고 시퀀스를 버리면 Python이 시퀀스를 해제한다. #16(ROOT TH1), #18(pyhepmc)과 같은 계열 —
**C++ 객체의 수명은 Python 변수의 수명을 따른다.**
**해결**: `jets, cs = cluster(...)`로 시퀀스를 변수에 붙잡아 둔다.

**증상 2**: `GhostedAreaSpec.set_random_status([s1, s2])` → `TypeError: argument 2 of type
'std::vector<int>'`. PyPI `fastjet` SWIG 바인딩에 정수 벡터 래퍼(`vectori`)가 없음.
**결과**: 고스트 위치가 실행마다 다르다 → 제트 면적이 ~1e-3 수준에서 흔들리고, Module 4
그림은 바이트 단위로 재현되지 않는다 (숫자 결론에는 영향 없음). `hicpy.jets.area_definition`은
바인딩이 받아 주면 seed를 쓰고, 아니면 그냥 넘어간다.
