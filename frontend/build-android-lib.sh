#!/usr/bin/env bash
# build-android-lib.sh — cross-compile the Mathilda kernel (+ GMP, MPFR) as
# static libraries for Android via the NDK, so the notebook runs the CAS
# in-process (Android, like iOS, forbids spawning the sidecar). Produces, per
# ABI:
#
#   frontend/src-tauri/gen/android-libs/<rust-target>/libmathilda.a
#   frontend/src-tauri/gen/android-libs/<rust-target>/libgmp.a
#   frontend/src-tauri/gen/android-libs/<rust-target>/libmpfr.a   (WITH_MPFR=1, default)
#
# `build.rs` links these when CARGO_CFG_TARGET_OS=android.
#
# REQUIREMENTS:
#   * Android NDK (set NDK_HOME, or it is auto-detected under
#     $ANDROID_HOME/ndk/*). Tested with r27.
#   * GMP and MPFR source tarballs (same as the iOS build):
#       GMP_TARBALL=gmp-6.3.0.tar.xz  MPFR_TARBALL=mpfr-4.2.1.tar.xz
#
# USAGE:
#   ./build-android-lib.sh arm64      # arm64-v8a  (aarch64-linux-android) — default
#   ./build-android-lib.sh x86_64     # x86_64     (emulator on Intel hosts)
#   ABI=arm64 API=24 ./build-android-lib.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC_TAURI="$SCRIPT_DIR/src-tauri"
OUT_ROOT="$SRC_TAURI/gen/android-libs"

# API 26 (Android 8.0): the kernel uses <complex.h> clog/cpow/csqrt, which
# Bionic gates behind __INTRODUCED_IN(26). The app's minSdk must be >= 26 too.
API="${API:-26}"
WITH_MPFR="${WITH_MPFR:-1}"            # see build-ios-lib.sh for why MPFR is on

# --- locate the NDK ---------------------------------------------------------
NDK="${NDK_HOME:-${ANDROID_NDK_HOME:-}}"
if [ -z "$NDK" ]; then
  sdk="${ANDROID_HOME:-${ANDROID_SDK_ROOT:-$HOME/Library/Android/sdk}}"
  NDK="$(ls -d "$sdk"/ndk/* 2>/dev/null | sort -V | tail -1 || true)"
fi
if [ -z "$NDK" ] || [ ! -d "$NDK" ]; then
  echo "ERROR: Android NDK not found. Set NDK_HOME (e.g. \$ANDROID_HOME/ndk/<ver>)." >&2
  exit 1
fi
HOSTTAG="$(ls -d "$NDK"/toolchains/llvm/prebuilt/* | head -1)"
TOOLBIN="$HOSTTAG/bin"

# --- pick ABI ---------------------------------------------------------------
ABI="${1:-${ABI:-arm64}}"
case "$ABI" in
  arm64|arm64-v8a|aarch64)
    TRIPLE="aarch64-linux-android";      RUST_TARGET="aarch64-linux-android" ;;
  x86_64)
    TRIPLE="x86_64-linux-android";       RUST_TARGET="x86_64-linux-android" ;;
  arm|armeabi-v7a|armv7)
    TRIPLE="armv7a-linux-androideabi";   RUST_TARGET="armv7-linux-androideabi" ;;
  x86)
    TRIPLE="i686-linux-android";         RUST_TARGET="i686-linux-android" ;;
  *) echo "usage: $0 [arm64|x86_64|arm|x86]" >&2; exit 2 ;;
esac

CC_BIN="$TOOLBIN/${TRIPLE}${API}-clang"
AR_BIN="$TOOLBIN/llvm-ar"
RANLIB_BIN="$TOOLBIN/llvm-ranlib"
OUT_DIR="$OUT_ROOT/$RUST_TARGET"

[ -x "$CC_BIN" ] || { echo "ERROR: NDK compiler not found: $CC_BIN" >&2; exit 1; }

echo ">>> NDK=$NDK"
echo ">>> abi=$ABI triple=$TRIPLE api=$API"
echo ">>> cc=$CC_BIN"
echo ">>> out=$OUT_DIR"
mkdir -p "$OUT_DIR"

# GMP/MPFR configure needs --host to differ from --build so autoconf enters
# cross mode (it must NOT try to run the Android test binaries on the host).
CONFIG_HOST="$TRIPLE"
ANDROID_CFLAGS="-O2 -fPIC"

# --- 1. GMP -----------------------------------------------------------------
build_gmp() {
  if [ -f "$OUT_DIR/libgmp.a" ] && [ "${FORCE_DEPS:-0}" != "1" ]; then
    echo ">>> libgmp.a already present (FORCE_DEPS=1 to rebuild)"; return 0
  fi
  local src="${GMP_SRC:-}"
  if [ -z "$src" ] && [ -n "${GMP_TARBALL:-}" ]; then
    src="$(mktemp -d)/gmp"; mkdir -p "$src"
    tar -xf "$GMP_TARBALL" -C "$src" --strip-components=1
  fi
  [ -n "$src" ] || { echo "ERROR: set GMP_SRC or GMP_TARBALL (gmp-6.3.0.tar.xz)." >&2; exit 1; }
  echo ">>> building GMP from $src"
  local bd="$OUT_DIR/gmp-build"; rm -rf "$bd"; mkdir -p "$bd"
  ( cd "$bd"
    CC="$CC_BIN" CFLAGS="$ANDROID_CFLAGS" AR="$AR_BIN" RANLIB="$RANLIB_BIN" \
      "$src/configure" --host="$CONFIG_HOST" \
        --disable-shared --enable-static --disable-assembly \
        --prefix="$bd/install"
    make -j"$(command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu)"; make install )
  cp "$bd/install/lib/libgmp.a" "$OUT_DIR/libgmp.a"
  cp "$bd/install/include/gmp.h" "$OUT_DIR/gmp.h"
  echo ">>> libgmp.a ready"
}

# --- 1b. MPFR ---------------------------------------------------------------
build_mpfr() {
  [ "$WITH_MPFR" = "1" ] || return 0
  if [ -f "$OUT_DIR/libmpfr.a" ] && [ "${FORCE_DEPS:-0}" != "1" ]; then
    echo ">>> libmpfr.a already present (FORCE_DEPS=1 to rebuild)"; return 0
  fi
  local src="${MPFR_SRC:-}"
  if [ -z "$src" ] && [ -n "${MPFR_TARBALL:-}" ]; then
    src="$(mktemp -d)/mpfr"; mkdir -p "$src"
    tar -xf "$MPFR_TARBALL" -C "$src" --strip-components=1
  fi
  [ -n "$src" ] || { echo "ERROR: WITH_MPFR=1 but MPFR_SRC/MPFR_TARBALL unset." >&2; exit 1; }
  echo ">>> building MPFR from $src"
  local gmp_install="$OUT_DIR/gmp-build/install"
  local bd="$OUT_DIR/mpfr-build"; rm -rf "$bd"; mkdir -p "$bd"
  ( cd "$bd"
    CC="$CC_BIN" CFLAGS="$ANDROID_CFLAGS" AR="$AR_BIN" RANLIB="$RANLIB_BIN" \
      "$src/configure" --host="$CONFIG_HOST" \
        --disable-shared --enable-static --with-gmp="$gmp_install" \
        --prefix="$bd/install"
    make -j"$(command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu)"; make install )
  cp "$bd/install/lib/libmpfr.a" "$OUT_DIR/libmpfr.a"
  cp "$bd/install/include/mpfr.h" "$OUT_DIR/mpfr.h"
  echo ">>> libmpfr.a ready"
}

# --- 2. libmathilda ---------------------------------------------------------
build_kernel() {
  echo ">>> building libmathilda.a for Android/$ABI"
  local mpfr_flag="USE_MPFR=0"
  if [ "$WITH_MPFR" = "1" ]; then
    mpfr_flag="USE_MPFR=1"
    [ -f "$OUT_DIR/libmpfr.a" ] || { echo "ERROR: libmpfr.a missing." >&2; exit 1; }
  fi
  # CC wrapper injects the NDK target + our cross-built headers, without
  # clobbering the makefile's -I include list (same pattern as the iOS build).
  local wrapper="$OUT_DIR/cc-android.sh"
  cat > "$wrapper" <<EOF
#!/bin/sh
exec "$CC_BIN" $ANDROID_CFLAGS -I"$OUT_DIR" "\$@"
EOF
  chmod +x "$wrapper"
  ( cd "$REPO_ROOT"
    make clean >/dev/null 2>&1 || true
    make libmathilda.a \
      CC="$wrapper" AR="$AR_BIN" \
      USE_READLINE=0 USE_THREADS=0 USE_ECM=0 $mpfr_flag \
      USE_LAPACK=0 USE_GRAPHICS=0 USE_FLINT=0 USE_REGEX=0 USE_FFTW=0 \
      -j"$(command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu)" )
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
echo "Next: npm --prefix \"$SCRIPT_DIR\" run android:dev"
