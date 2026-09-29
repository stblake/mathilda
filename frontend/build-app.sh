#!/bin/sh
# build-app.sh
# Builds the whole clickable desktop app, in the order the pieces depend on each
# other, and leaves the double-clickable launcher in place.
#
#   1. the Mathilda kernel, installed as the Tauri sidecar (build-sidecar.sh)
#   2. the JS dependencies, if they are not already there
#   3. the release bundle (`tauri build`), which runs vite via beforeBuildCommand
#   4. MathildaNotebook.command at the repo root, made executable, as the thing
#      you click
#
# ONE path, called from two places: `npm run build:app`, and the launcher itself
# when the bundle is missing. Keeping the sequence here rather than duplicating it
# in the launcher means the two cannot drift — the launcher's job is to open the
# app, not to know how it is built.
#
# `--bundles app` builds the .app only, not the .dmg. The .dmg is for shipping to
# someone else and costs minutes more; `npm run tauri build` (no flags) still
# builds every target when that is what you want.
set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
LAUNCHER="$REPO_ROOT/MathildaNotebook.command"

command -v cargo >/dev/null 2>&1 || { echo "ERROR: cargo not found; install Rust (https://rustup.rs)" >&2; exit 1; }
command -v npm   >/dev/null 2>&1 || { echo "ERROR: npm not found; install Node" >&2; exit 1; }

"$SCRIPT_DIR/build-sidecar.sh"

cd "$SCRIPT_DIR"
[ -d node_modules ] || npm install

npm run tauri build -- --bundles app

# The launcher is checked in and git tracks its executable bit, so this is
# belt-and-braces: an editor or a copy that drops the bit turns a double-click
# into "open in TextEdit", which looks like the launcher is broken.
if [ -f "$LAUNCHER" ]; then
    chmod +x "$LAUNCHER"
else
    echo "WARNING: $LAUNCHER is missing — the clickable launcher will not be there." >&2
fi

APP="$SCRIPT_DIR/src-tauri/target/release/bundle/macos/Mathilda.app"
echo
if [ -d "$APP" ]; then
    echo "Built: $APP"
else
    echo "Built (see $SCRIPT_DIR/src-tauri/target/release/bundle/ for the artefact)."
fi
[ -f "$LAUNCHER" ] && echo "Open it by double-clicking: ${LAUNCHER#"$REPO_ROOT"/}"
exit 0
