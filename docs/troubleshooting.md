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
