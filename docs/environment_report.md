# Phase 0 — 환경 조사 보고서

작성: 2026-08-27

## 1. 이 세션에서 실제로 접근 가능한 계산 자원

Claude가 명령을 실행할 수 있는 곳은 두 곳이며, **둘 다 Windows 본체가 아닙니다.**

| | (A) 클라우드 컨테이너 | (B) Cowork 로컬 리눅스 VM |
|---|---|---|
| OS | Ubuntu 24.04.4 LTS | Ubuntu 22.04.5 LTS |
| 커널 | 6.18.44 | 6.8.0 |
| CPU | Intel Xeon @2.80 GHz, 2 논리코어 | AMD Ryzen 7 7800X3D 중 2 논리코어 |
| RAM | 7.8 GiB | 3.8 GiB |
| 디스크 여유 | ~30 GB | ~8.6 GB (`/sessions`) |
| Python | 3.11.15 | 3.10.12 |
| gcc/g++ | 13.3.0 | 11.x |
| cmake | 3.28.3 | **없음** |
| 네트워크 | PyPI ✅ / GitHub ✅ / conda-forge ❌ / root.cern ❌ / opendata.cern.ch ❌ / hepdata.net ❌ | **전면 차단** (프록시 403) |
| 지속성 | 세션 종료 시 소멸 | 세션 종료 시 소멸 (`$HOME`), `mnt/`만 Windows 디스크에 남음 |

핵심 결론 두 가지:

1. **클라우드 컨테이너에는 ROOT를 설치할 수 없습니다.** conda-forge와 root.cern이
   egress allowlist에서 막혀 있고, ROOT는 PyPI wheel로 배포되지 않습니다.
   같은 이유로 CERN Open Data / HEPData 직접 다운로드도 불가능합니다.
2. **Windows 본체의 8코어 Ryzen과 실제 디스크는 WSL2에서만 온전히 쓸 수 있습니다.**
   위 두 환경은 각각 2코어만 노출합니다.

## 2. 채택한 아키텍처: 이중 트랙

```
   [Track A] Windows 11 + WSL2 Ubuntu + micromamba   ← 장기 연구 환경 (본체)
        │     ROOT · PYTHIA8 · FastJet · HepMC3 · (나중에 LHAPDF/MadGraph/Delphes)
        │     8코어 전부 사용, 디스크 제한 없음, 실제 데이터 다운로드 가능
        │
        ├── 같은 git 저장소 ──────────────────────────┐
        │                                             │
   [Track B] 클라우드 컨테이너 (Claude 작업 공간)      │
              PyPI 기반 Scikit-HEP 스택으로 코드를 실제로 실행·검증
              uproot · awkward · vector · hist · fastjet(3.5.1) · pyhepmc(HepMC3)
```

- Claude가 짜는 모든 코드는 Track B에서 **먼저 실행해 보고** 저장소에 올립니다.
  실행하지 못한 코드(ROOT C++ 매크로 등)는 파일 상단에 `STATUS: NOT executed`로
  명시합니다.
- Track A는 Unoh님이 `bash environment/setup_wsl2.sh` 한 번으로 만듭니다.
  출력 로그를 주시면 그 자리에서 디버깅합니다.

이렇게 나눈 이유: "설치는 사용자가, 검증은 Claude가"가 아니라, **코드의 정확성은
지금 즉시 검증하고, 무거운 바이너리 스택은 진짜 쓸 곳에 한 번만 설치**하기
위해서입니다.

## 3. Track B에서 실제로 검증된 스택 (2026-08-27)

```
python 3.11.15   numpy 2.4.6   scipy 1.17.1   matplotlib 3.11.1
uproot 5.7.6     awkward 2.13.0   vector 1.8.1   hist 2.11.0
mplhep 1.3.3     particle 1.0.0   fastjet 3.5.1.4   pyhepmc 2.16.1 (HepMC3)
gcc/g++ 13.3.0   cmake 3.28.3     git 2.43.0
```

전체 기록: `environment/versions_installed.txt` (`verify_environment.py`가 생성).

**Track B에서 안 되는 것**: ROOT, PYTHIA 8, LHAPDF, MadGraph, Delphes, Rivet.
PYTHIA는 PyPI에 유지되는 바인딩이 없어(`numpythia` 빌드 실패, 오답노트 #3 참조)
Module 2부터는 Track A가 필수입니다.

## 4. Track A 설치 시 확인할 것

1. WSL2가 이미 있는지: PowerShell에서 `wsl -l -v` → `VERSION 2`인 Ubuntu가 있으면 재사용
2. 디스크: ROOT+PYTHIA+FastJet conda 환경이 약 3–5 GB, 이후 데이터까지 감안해 30 GB 이상 여유 권장
3. `.wslconfig`로 메모리/코어 제한이 걸려 있지 않은지 (기본은 호스트의 50 %)
4. 저장소는 **WSL 파일시스템 안**(`~/projects/...`)에 두는 것이 `/mnt/c/...`보다 I/O가 훨씬 빠릅니다.
   다만 이 저장소는 Windows 쪽(`C:\Users\Administrator\Claude\Projects\HIC analysis`)에
   있어 Claude가 접근할 수 있으므로, 당분간은 `/mnt/c/...` 경로로 작업하고
   대용량 이벤트 생성만 WSL 내부 디스크에서 하는 편이 현실적입니다.
