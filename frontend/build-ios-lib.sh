#!/usr/bin/env bash
# build-ios-lib.sh — cross-compile the Mathilda kernel (+ GMP) as static
# libraries for iOS, so the notebook can run the CAS in-process (no sidecar,
# which iOS forbids). Produces, per target:
#
#   frontend/src-tauri/gen/ios-libs/<target>/libmathilda.a
#   frontend/src-tauri/gen/ios-libs/<target>/libgmp.a
#   frontend/src-tauri/gen/ios-libs/<target>/libmpfr.a   (when WITH_MPFR=1, default)
#
# `build.rs` links these when CARGO_CFG_TARGET_OS=ios; point it at the right
# directory with MATHILDA_IOS_LIBDIR (see the tauri ios wrapper in package.json).
#
# WHY MPFR IS ON BY DEFAULT: the kernel's linalg/numerical subsystems reference
# MPFR unconditionally today (a USE_MPFR=0 build does not yet compile — tracked
# as follow-up work). So iOS uses the same GMP+MPFR config that the desktop
# build and the FFI smoke test already validate. Set WITH_MPFR=0 only once the
# minimal-config build is fixed upstream.
#
# REQUIREMENTS (cannot be satisfied by Command Line Tools alone):
#   * Full Xcode + an iOS Simulator runtime.
#       sudo xcode-select -s /Applications/Xcode.app
#       sudo xcodebuild -license accept
#       xcodebuild -downloadPlatform iOS
#   * GMP and MPFR source tarballs (neither cross-compiles via brew).
#       curl -LO https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
#       curl -LO https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
#     Point GMP_SRC / MPFR_SRC at the extracted dirs, or set GMP_TARBALL /
#     MPFR_TARBALL and this script extracts them.
#
# USAGE:
#   ./build-ios-lib.sh sim      # arm64 iOS Simulator (Apple Silicon macs)
#   ./build-ios-lib.sh device   # arm64 iOS device (for TestFlight/App Store)
#   TARGET=sim WITH_MPFR=0 ./build-ios-lib.sh   # (only after minimal build fixed)
#
set -euo pipefail

# --- resolve paths ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"          # the Mathilda repo root
SRC_TAURI="$SCRIPT_DIR/src-tauri"
OUT_ROOT="$SRC_TAURI/gen/ios-libs"

# --- pick target ------------------------------------------------------------
TARGET="${1:-${TARGET:-sim}}"
case "$TARGET" in
  sim)
    SDK="iphonesimulator"
    ARCH="arm64"
    MINFLAG="-mios-simulator-version-min=13.0"
    GMP_HOST="aarch64-apple-darwin"      # sim is arm64-darwin-like
    OUT_DIR="$OUT_ROOT/aarch64-apple-ios-sim" ;;
  device)
    SDK="iphoneos"
    ARCH="arm64"
    MINFLAG="-miphoneos-version-min=13.0"
    GMP_HOST="aarch64-apple-darwin"
    OUT_DIR="$OUT_ROOT/aarch64-apple-ios" ;;
  *)
    echo "usage: $0 [sim|device]" >&2; exit 2 ;;
esac

WITH_MPFR="${WITH_MPFR:-1}"

# --- sanity: full Xcode present --------------------------------------------
if ! xcrun --sdk "$SDK" --show-sdk-path >/dev/null 2>&1; then
  cat >&2 <<EOF
ERROR: the "$SDK" SDK cannot be located.

This needs full Xcode (Command Line Tools are not enough). Install Xcode from
the App Store, then:
  sudo xcode-select -s /Applications/Xcode.app
  sudo xcodebuild -license accept
  xcodebuild -downloadPlatform iOS
EOF
  exit 1
fi

SYSROOT="$(xcrun --sdk "$SDK" --show-sdk-path)"
CC_BIN="$(xcrun --sdk "$SDK" --find clang)"
AR_BIN="$(xcrun --sdk "$SDK" --find ar)"
RANLIB_BIN="$(xcrun --sdk "$SDK" --find ranlib)"
IOS_CFLAGS="-arch $ARCH -isysroot $SYSROOT $MINFLAG -O2 -fembed-bitcode-marker"

echo ">>> target=$TARGET sdk=$SDK arch=$ARCH"
echo ">>> sysroot=$SYSROOT"
echo ">>> out=$OUT_DIR"
mkdir -p "$OUT_DIR"

# --- 1. GMP -----------------------------------------------------------------
# GMP must be cross-compiled; its assembly does not build for the iOS toolchain,
# so we disable it (--disable-assembly) at a modest perf cost.
build_gmp() {
  local src="${GMP_SRC:-}"
  if [ -z "$src" ] && [ -n "${GMP_TARBALL:-}" ]; then
    src="$(mktemp -d)/gmp"
    mkdir -p "$src"
    tar -xf "$GMP_TARBALL" -C "$src" --strip-components=1
  fi
  if [ -z "$src" ]; then
    echo "ERROR: set GMP_SRC (extracted gmp dir) or GMP_TARBALL (gmp-*.tar.xz)." >&2
    echo "       Download: https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz" >&2
    exit 1
  fi
  echo ">>> building GMP from $src"
  local build_dir="$OUT_DIR/gmp-build"
  rm -rf "$build_dir"; mkdir -p "$build_dir"
  ( cd "$build_dir"
    CC="$CC_BIN" CFLAGS="$IOS_CFLAGS" AR="$AR_BIN" RANLIB="$RANLIB_BIN" \
      "$src/configure" \
        --host="$GMP_HOST" \
        --disable-shared --enable-static --disable-assembly \
        --prefix="$build_dir/install"
    make -j"$(sysctl -n hw.ncpu)"
    make install
  )
  cp "$build_dir/install/lib/libgmp.a" "$OUT_DIR/libgmp.a"
  cp "$build_dir/install/include/gmp.h" "$OUT_DIR/gmp.h"
  echo ">>> libgmp.a ready"
}

# --- 1b. MPFR (against the cross-built GMP) ---------------------------------
build_mpfr() {
  [ "$WITH_MPFR" = "1" ] || return 0
  local src="${MPFR_SRC:-}"
  if [ -z "$src" ] && [ -n "${MPFR_TARBALL:-}" ]; then
    src="$(mktemp -d)/mpfr"
    mkdir -p "$src"
    tar -xf "$MPFR_TARBALL" -C "$src" --strip-components=1
  fi
  if [ -z "$src" ]; then
    echo "ERROR: WITH_MPFR=1 but MPFR_SRC / MPFR_TARBALL not set." >&2
    echo "       Download: https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz" >&2
    exit 1
  fi
  echo ">>> building MPFR from $src (against iOS GMP)"
  local gmp_install="$OUT_DIR/gmp-build/install"
  local build_dir="$OUT_DIR/mpfr-build"
  rm -rf "$build_dir"; mkdir -p "$build_dir"
  ( cd "$build_dir"
    CC="$CC_BIN" CFLAGS="$IOS_CFLAGS" AR="$AR_BIN" RANLIB="$RANLIB_BIN" \
      "$src/configure" \
        --host="$GMP_HOST" \
        --disable-shared --enable-static \
        --with-gmp="$gmp_install" \
        --prefix="$build_dir/install"
    make -j"$(sysctl -n hw.ncpu)"
    make install
  )
  cp "$build_dir/install/lib/libmpfr.a" "$OUT_DIR/libmpfr.a"
  cp "$build_dir/install/include/mpfr.h" "$OUT_DIR/mpfr.h"
  echo ">>> libmpfr.a ready"
}

# --- 2. libmathilda ---------------------------------------------------------
# Minimal, dependency-light kernel: GMP required, everything else off. The C
# FFI entry points live in src/ffi/mathilda_ffi.c (archived via the makefile's
# `libmathilda.a` target, which excludes repl.o).
build_kernel() {
  echo ">>> building libmathilda.a for iOS"
  local mpfr_flag="USE_MPFR=0"
  local extra_inc="-I$OUT_DIR"          # our cross-built gmp.h / mpfr.h live here
  if [ "$WITH_MPFR" = "1" ]; then
    mpfr_flag="USE_MPFR=1"
    if [ ! -f "$OUT_DIR/libmpfr.a" ]; then
      echo "ERROR: WITH_MPFR=1 but $OUT_DIR/libmpfr.a is missing (build_mpfr failed?)." >&2
      exit 1
    fi
  fi
  ( cd "$REPO_ROOT"
    make clean >/dev/null 2>&1 || true
    make libmathilda.a \
      CC="$CC_BIN" \
      AR="$AR_BIN" \
      CFLAGS="$IOS_CFLAGS -std=c99 -I./src -I./src/ffi $extra_inc" \
      USE_READLINE=0 USE_THREADS=0 USE_ECM=0 $mpfr_flag \
      USE_LAPACK=0 USE_GRAPHICS=0 USE_FLINT=0 USE_REGEX=0 USE_FFTW=0 \
      -j"$(sysctl -n hw.ncpu)"
  )
  cp "$REPO_ROOT/libmathilda.a" "$OUT_DIR/libmathilda.a"
  echo ">>> libmathilda.a ready"
}

build_gmp
build_mpfr
build_kernel

echo
echo "DONE. Artifacts in $OUT_DIR:"
ls -la "$OUT_DIR"/*.a
echo
echo "Next: run the iOS app with the lib dir exported, e.g."
echo "  MATHILDA_IOS_LIBDIR=\"$OUT_DIR\" WITH_MPFR=$WITH_MPFR \\"
echo "    npm --prefix \"$SCRIPT_DIR\" run tauri ios dev"
