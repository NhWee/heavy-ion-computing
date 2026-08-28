# GitHub 연동

원격: `https://github.com/NhWee/heavy-ion-computing.git` (이미 `origin`으로 등록됨)

## 1. 저장소 생성 (한 번만, 직접)

Claude의 API 토큰은 **세션에 연결된 저장소로만 스코프가 제한**되어 있어 새 저장소를
만들 수 없습니다. 아래에서 직접 만들어 주세요.

<https://github.com/new> 에서:

| 항목 | 값 |
|---|---|
| Repository name | `heavy-ion-computing` |
| Description | Building a heavy-ion / collider computational workflow from scratch: ROOT, PYTHIA 8, FastJet, HepMC, and comparison with public LHC data |
| Visibility | **Private** 권장 (나중에 Settings에서 Public 전환 가능) |
| Initialize this repository with | **아무것도 체크하지 않기** (README·.gitignore·license 모두 해제) |

마지막 항목이 중요합니다. 초기 커밋이 생기면 로컬 히스토리와 갈라져서
`git pull --allow-unrelated-histories`가 필요해집니다.

## 2. 첫 push (WSL2에서)

```bash
cd "/mnt/c/Users/Administrator/Claude/Projects/HIC analysis/heavy-ion-computing"
git push -u origin main
```

인증이 필요하면 셋 중 하나:

**(a) Windows Git Credential Manager 재사용** — 가장 간단합니다.
```bash
git config --global credential.helper "/mnt/c/Program\\ Files/Git/mingw64/bin/git-credential-manager.exe"
```
Windows에 Git for Windows가 설치되어 있어야 합니다.

**(b) gh CLI**
```bash
micromamba install -y -n hic gh   # 또는 시스템 apt
gh auth login                     # → GitHub.com → HTTPS → 브라우저 인증
git push -u origin main
```

**(c) Personal Access Token** — <https://github.com/settings/tokens> 에서
`repo` 권한 fine-grained token을 만들고, push 시 비밀번호 자리에 붙여넣습니다.
```bash
git config --global credential.helper store   # 평문 저장이므로 개인 머신에서만
```

## 3. 커밋 규약

- 커밋 메시지는 **무엇을 왜** 바꿨는지 씁니다. "fix", "update" 단독 금지.
- 물리 결과가 바뀌는 커밋에는 바뀐 숫자를 본문에 적습니다
  (예: `pion closure fit chi2/ndf 0.91 -> 0.50`).
- 생성된 ROOT 파일(`data/generated/*.root`)은 커밋하지 않습니다 —
  `config/analysis.yaml` + seed + 코드로 재생성됩니다.
  대신 provenance JSON(`*.json`)은 커밋합니다.
- 그림(`results/figures/*.pdf|png`)은 커밋합니다. 결과의 기록이기 때문입니다.

## 4. 무엇을 어디에 남기나

| | GitHub | Notion |
|---|---|---|
| 내용 | 코드, 설정, 재현 방법, 결과 수치 | 연구 노트, 판단 근거, 실패 기록 |
| 이유 | 실행 가능한 상태를 정확히 보존 | 왜 그렇게 했는지는 코드에 안 남는다 |

Notion 허브: <https://app.notion.com/p/3c96974402a681d5bbfdd050e0af6ebe>
- 작업 로그 (Work Log) — 세션별 기록, 커밋 해시와 seed 포함
- 실패 보고 (Failure Report) — 오답노트, 심각도/영역/가드 위치 포함
