#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# One-shot environment bootstrap for heavy-ion-computing on WSL2 (Ubuntu).
#
#   bash environment/setup_wsl2.sh
#
# What it does, and why:
#   1. installs a handful of system libraries ROOT's graphics/X11 layer needs
#      (ROOT itself comes from conda-forge, but it links against these);
#   2. installs micromamba -- a single static binary, no base environment, no
#      `conda init` pollution of your shell;
#   3. creates the `hic` environment from environment/environment.yml;
#   4. runs the validation script, which is the only thing that proves the
#      installation actually works.
#
# It is idempotent: re-running it is safe.
# ---------------------------------------------------------------------------
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV_NAME="hic"

echo "=== [1/4] system libraries (sudo required) ==="
sudo apt-get update -qq
sudo apt-get install -y --no-install-recommends \
    build-essential curl git ca-certificates bzip2 \
    libx11-6 libxext6 libxft2 libxpm4 libglu1-mesa libgl1

echo "=== [2/4] micromamba ==="
if ! command -v micromamba >/dev/null 2>&1 && [ ! -x "$HOME/.local/bin/micromamba" ]; then
    mkdir -p "$HOME/.local/bin"
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xvj -C "$HOME/.local" bin/micromamba
fi
export PATH="$HOME/.local/bin:$PATH"
export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
micromamba --version

echo "=== [3/4] create/update the '$ENV_NAME' environment ==="
if micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    micromamba install -y -n "$ENV_NAME" -f "$REPO/environment/environment.yml"
else
    micromamba create -y -f "$REPO/environment/environment.yml"
fi

echo "=== [4/4] validation ==="
micromamba run -n "$ENV_NAME" python "$REPO/environment/verify_environment.py" \
    --write "$REPO/environment/versions_installed.txt"

cat <<'MSG'

-----------------------------------------------------------------------
Done.  Add these two lines to ~/.bashrc so the environment is available
in every new shell:

    export PATH="$HOME/.local/bin:$PATH"
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
    eval "$(micromamba shell hook --shell bash)"

Then:   micromamba activate hic
-----------------------------------------------------------------------
MSG
