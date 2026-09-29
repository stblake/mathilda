#!/bin/sh
# build-sidecar.sh
# Builds the Mathilda C binary and copies it into the Tauri sidecar
# binaries directory with the required target-triple suffix.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARIES_DIR="$SCRIPT_DIR/src-tauri/binaries"

# Get the current Rust target triple.
TARGET=$(rustc -vV | grep host | awk '{print $2}')
if [ -z "$TARGET" ]; then
    echo "ERROR: could not determine Rust target triple" >&2
    exit 1
fi

echo "Building Mathilda for target: $TARGET"

# Build the Mathilda binary (without ECM: the bundled app must not need a libecm
# the user has not installed).
#
# The two sources that read -DNO_ECM are forced to recompile first. Object files
# do not record the flags they were built with, so after a plain `make` (which
# defaults to USE_ECM=1) these two are up to date as far as make is concerned,
# and the USE_ECM=0 link then fails on a missing _ecm_init -- the objects still
# carry the calls the flag was supposed to compile out. Touching them is the
# whole fix; everything else in the tree is flag-independent.
touch "$REPO_ROOT/src/facint.c" "$REPO_ROOT/src/version.c"
make -C "$REPO_ROOT" USE_ECM=0 -j4

mkdir -p "$BINARIES_DIR"
cp "$REPO_ROOT/Mathilda" "$BINARIES_DIR/mathilda-$TARGET"
echo "Sidecar installed: $BINARIES_DIR/mathilda-$TARGET"

# The makefile writes one output path, so the USE_ECM=0 build above has REPLACED
# the repo's ./Mathilda. Said out loud rather than left to be discovered: the next
# person to use ./Mathilda for a factorisation would get the degraded path with no
# indication of why, and `$Version` is the only visible difference.
echo "NOTE: $REPO_ROOT/Mathilda is now the USE_ECM=0 build (the app must not need libecm)."
echo "      Restore the default with:  make -C '$REPO_ROOT' -j"

# ALSO refresh the copies cargo already made, because they are what a running dev app executes.
#
# Tauri copies the sidecar into target/{debug,release}/mathilda during the cargo build. If the Rust
# source has not changed, cargo skips the build, the copy step never re-runs, and `cargo tauri dev`
# keeps launching a kernel from whenever the Rust last compiled. That is not a hypothetical: a
# volume-rendering change was verified against a freshly built binary and reported as working while
# the app was running a kernel over an hour old, so every Image3D result on screen came from the
# previous payload format and looked like a bug in the new code. Rebuilding the kernel has to update
# what the app will actually run, not only what a future packaging step would pick up.
# Replace ATOMICALLY, via a temporary file and mv. Copying onto the path in place is not merely
# untidy: on macOS, overwriting a Mach-O that a running process has mapped invalidates its code
# signature, and the result is a binary the kernel SIGKILLs on exec (exit 137) -- so the "refreshed"
# kernel would not start at all, which is a worse failure than the stale one it was fixing. A mv
# gives the new file its own inode and leaves the running image alone.
for D in "$SCRIPT_DIR/src-tauri/target/debug" "$SCRIPT_DIR/src-tauri/target/release"; do
    if [ -f "$D/mathilda" ]; then
        cp "$REPO_ROOT/Mathilda" "$D/mathilda.new"
        # Ad-hoc sign the replacement: an unsigned or stale-signature binary is killed on exec on
        # Apple Silicon, and codesign is a no-op on platforms without it.
        command -v codesign >/dev/null 2>&1 && codesign --force --sign - "$D/mathilda.new" 2>/dev/null || true
        mv -f "$D/mathilda.new" "$D/mathilda"
        echo "Refreshed dev copy:  $D/mathilda"
    fi
done
