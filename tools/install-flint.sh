#!/usr/bin/env bash
#
# install-flint.sh — build and install FLINT >= 3.0 from source.
#
# WHY THIS EXISTS
#   Mathilda's accelerated algebraic-extension arithmetic (multivariate GCD/
#   factoring over Q, number-field GCD/factoring over Q(alpha), RootReduce via
#   qqbar, rigorous acb numerics) needs FLINT >= 3.0 — the release that merged
#   ANTIC, Arb, and Calcium into FLINT and introduced the generic-ring `gr`
#   layer. Ubuntu 22.04 / Debian Bullseye only package FLINT 2.x through apt,
#   which lacks all of that, so `sudo apt install libflint-dev` gives a version
#   Mathilda's `pkg-config --exists 'flint >= 3.0'` gate rejects (it then builds
#   fine with USE_FLINT=0 and the classical fallback — just slower).
#
#   The good news: Ubuntu 22.04 already ships GMP 6.2.1 and MPFR 4.1.0, which
#   meet FLINT 3's minimums, so a source build is quick and clean. This script
#   automates it. (Alternative with no build at all: conda-forge —
#   `conda install -c conda-forge libflint` — see docs/building.md.)
#
# USAGE
#   ./tools/install-flint.sh                 # build the default version into /usr/local (needs sudo)
#   FLINT_VERSION=3.1.3 ./tools/install-flint.sh
#   PREFIX="$HOME/.local" ./tools/install-flint.sh   # user-local install, no sudo
#
# After it finishes it prints the exact PKG_CONFIG_PATH / LD_LIBRARY_PATH to
# export if you installed outside the default system paths, then you just
# `make` Mathilda as usual and the FLINT paths light up automatically.

set -euo pipefail

# A conservative, widely-published default; override for a newer one. Mathilda
# supports the whole 3.x line (>= 3.0), so any recent release works.
FLINT_VERSION="${FLINT_VERSION:-3.1.3}"
PREFIX="${PREFIX:-/usr/local}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
WORKDIR="${WORKDIR:-$(mktemp -d)}"
TARBALL="flint-${FLINT_VERSION}.tar.gz"
# flintlib.org is the canonical host; the GitHub release mirror is the fallback.
URLS=(
  "https://flintlib.org/download/${TARBALL}"
  "https://github.com/flintlib/flint/releases/download/v${FLINT_VERSION}/${TARBALL}"
)

# sudo only when the prefix is not writable by the current user.
SUDO=""
if [ ! -w "$PREFIX" ] && [ "$(id -u)" -ne 0 ]; then SUDO="sudo"; fi

say() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
die() { printf '\033[1;31merror:\033[0m %s\n' "$*" >&2; exit 1; }

command -v cc >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1 \
  || die "no C compiler found — install build-essential first"
command -v make >/dev/null 2>&1 || die "make not found — install build-essential first"

# FLINT 3 hard-requires GMP >= 6.2.1 and MPFR >= 4.1.0. Ubuntu 22.04 satisfies
# both; warn early rather than after a long configure if the -dev headers are
# missing (apt: libgmp-dev libmpfr-dev).
for pc in gmp mpfr; do
  if command -v pkg-config >/dev/null 2>&1 && ! pkg-config --exists "$pc"; then
    printf '\033[1;33mwarning:\033[0m %s development files not found; ' "$pc"
    printf 'FLINT needs them. Try: sudo apt-get install -y libgmp-dev libmpfr-dev\n'
  fi
done

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists 'flint >= 3.0'; then
  say "FLINT $(pkg-config --modversion flint) (>= 3.0) is already visible to pkg-config — nothing to do."
  say "If Mathilda still reports it missing, check PKG_CONFIG_PATH."
  exit 0
fi

say "Building FLINT ${FLINT_VERSION} -> ${PREFIX} (${JOBS} jobs) in ${WORKDIR}"
cd "$WORKDIR"

fetched=""
for url in "${URLS[@]}"; do
  say "Fetching ${url}"
  if command -v curl >/dev/null 2>&1; then
    curl -fSL -o "$TARBALL" "$url" && { fetched=1; break; } || true
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$TARBALL" "$url" && { fetched=1; break; } || true
  else
    die "need curl or wget to download FLINT"
  fi
done
[ -n "$fetched" ] || die "could not download FLINT ${FLINT_VERSION} from any mirror (check the version number / network)"

tar xzf "$TARBALL"
cd "flint-${FLINT_VERSION}"

# --disable-static keeps the build fast and the install small; the shared lib is
# all Mathilda links against.
say "Configuring"
./configure --prefix="$PREFIX" --disable-static
say "Compiling (this is the slow part — a few minutes)"
make -j"$JOBS"
say "Installing to ${PREFIX} ${SUDO:+(with sudo)}"
$SUDO make install

# Refresh the runtime linker cache for a system prefix so the shared lib is found
# at run time without setting LD_LIBRARY_PATH.
if [ "$PREFIX" = "/usr/local" ] || [ "$PREFIX" = "/usr" ]; then
  $SUDO ldconfig 2>/dev/null || true
fi

say "Done. Installed FLINT ${FLINT_VERSION} to ${PREFIX}."

# Tell the user exactly what to export if the prefix is non-standard, so
# Mathilda's makefile (which uses pkg-config) can see it.
if [ "$PREFIX" != "/usr/local" ] && [ "$PREFIX" != "/usr" ]; then
  cat <<EOF

Add these to your shell profile so pkg-config and the loader find this install:

  export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:\${PKG_CONFIG_PATH:-}"
  export LD_LIBRARY_PATH="${PREFIX}/lib:\${LD_LIBRARY_PATH:-}"

EOF
fi

if command -v pkg-config >/dev/null 2>&1; then
  export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  if pkg-config --exists 'flint >= 3.0'; then
    say "Verified: pkg-config now sees flint $(pkg-config --modversion flint). Rebuild Mathilda with 'make' and USE_FLINT lights up automatically."
  else
    say "Installed, but pkg-config still can't see 'flint >= 3.0'. Export the PKG_CONFIG_PATH above, then re-check with: pkg-config --modversion flint"
  fi
fi
