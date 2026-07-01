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

# Build the Mathilda binary (without ECM by default).
make -C "$REPO_ROOT" USE_ECM=0 -j4

mkdir -p "$BINARIES_DIR"
cp "$REPO_ROOT/Mathilda" "$BINARIES_DIR/mathilda-$TARGET"

echo "Sidecar installed: $BINARIES_DIR/mathilda-$TARGET"
