#!/usr/bin/env bash
# Module 10: fetch and build the viscous hydrodynamics code MUSIC and the
# Cooper-Frye particle sampler iSS at pinned commits, into external/ (gitignored).
#
#   bash scripts/m10_build_hydro.sh          # inside the hic conda env (cmake, g++, gsl, zlib)
#
# MUSIC  (GPL-2.0)  https://github.com/MUSIC-fluid/MUSIC
#        Schenke, Jeon, Gale, PRC 82 014903 (2010) [arXiv:1004.1408];
#        PRL 106 042301 (2011) [arXiv:1009.3244]; Paquet et al. PRC 93 044906 (2016) [arXiv:1509.06738]
# iSS    (MIT)      https://github.com/chunshen1987/iSS   (C. Shen et al.)
#
# Writes external/hydro_versions.json (commit hashes) which the analysis copies
# into every result file.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
EXT=$ROOT/external
mkdir -p "$EXT"

MUSIC_URL=https://github.com/MUSIC-fluid/MUSIC
MUSIC_REV=85f435139d15bd67375872c490120a9b2cde6155
ISS_URL=https://github.com/chunshen1987/iSS
ISS_REV=8b3e43b26dade8252e2671907518a0d111625fd1

for t in git cmake make; do
  command -v $t >/dev/null || { echo "ERROR: $t not found (activate the hic env: conda activate hic)"; exit 1; }
done
command -v c++ >/dev/null || command -v g++ >/dev/null || { echo "ERROR: no C++ compiler (conda activate hic)"; exit 1; }

# MUSIC and iSS declare cmake_minimum_required below 3.5; CMake >= 4 refuses that
# unless told the minimum policy version (older CMake ignores this variable).
# conda provides GSL/zlib under $CONDA_PREFIX; cmake's FindGSL looks at GSL_ROOT_DIR
CMAKE_EXTRA=("-DCMAKE_POLICY_VERSION_MINIMUM=3.5")
if [[ -n "${CONDA_PREFIX:-}" && -f "$CONDA_PREFIX/include/gsl/gsl_math.h" ]]; then
  CMAKE_EXTRA+=("-DGSL_ROOT_DIR=$CONDA_PREFIX" "-DCMAKE_PREFIX_PATH=$CONDA_PREFIX")
elif [[ ! -f /usr/include/gsl/gsl_math.h ]]; then
  echo "ERROR: GSL headers not found (conda install -c conda-forge gsl, or apt install libgsl-dev)"; exit 1
fi

fetch() {  # url rev dir
  local url=$1 rev=$2 dir=$3
  if [[ ! -d "$dir/.git" ]]; then
    git clone -q "$url" "$dir"
  fi
  (cd "$dir" && { git cat-file -e "$rev^{commit}" 2>/dev/null || git fetch -q origin; } && git -c advice.detachedHead=false checkout -q "$rev")
  echo "$(basename "$dir") at $(cd "$dir" && git rev-parse HEAD)"
}

build() {  # dir target
  local dir=$1 target=$2
  mkdir -p "$dir/build"
  (cd "$dir/build" && cmake .. -DCMAKE_BUILD_TYPE=Release "${CMAKE_EXTRA[@]}" > cmake.log 2>&1 \
     || { tail -20 cmake.log; exit 1; })
  (cd "$dir/build" && make -j"$(nproc)" > make.log 2>&1 || { tail -30 make.log; exit 1; })
  test -x "$dir/build/src/$target" || { echo "ERROR: $target not built"; exit 1; }
  echo "built $dir/build/src/$target"
}

fetch "$MUSIC_URL" "$MUSIC_REV" "$EXT/MUSIC"
fetch "$ISS_URL"   "$ISS_REV"   "$EXT/iSS"
build "$EXT/MUSIC" MUSIChydro
build "$EXT/iSS"   iSS.e

cat > "$EXT/hydro_versions.json" <<EOF
{
  "MUSIC": {"url": "$MUSIC_URL", "commit": "$(cd "$EXT/MUSIC" && git rev-parse HEAD)", "license": "GPL-2.0"},
  "iSS":   {"url": "$ISS_URL",   "commit": "$(cd "$EXT/iSS" && git rev-parse HEAD)", "license": "MIT"},
  "built": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "compiler": "$( (c++ --version || g++ --version) 2>/dev/null | head -1)"
}
EOF
echo "OK: external/hydro_versions.json"
