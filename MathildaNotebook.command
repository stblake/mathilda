#!/bin/sh
# MathildaNotebook.command — double-click to open the Mathilda notebook.
#
# It lives at the REPOSITORY ROOT, beside the makefile — the first thing you see
# when you open the folder, and it resolves everything from its own location, so
# it works however it is invoked and wherever the repo is checked out.
#
# WHY THE .command SUFFIX. It is what makes the file clickable: Finder runs a
# .command file in Terminal, where a plain executable with no extension opens in
# a text editor instead. Finder hides the extension when "show all filename
# extensions" is off, so it reads as "MathildaNotebook" in the folder.
#
# WHY A SCRIPT AND NOT A SYMLINK to the built app. The bundle lives under
# frontend/src-tauri/target/, which is ephemeral: `cargo clean` or a fresh clone
# removes it, and a symlink would then be a dead file that reports nothing when
# clicked. This (re)builds the bundle whenever it is missing OR older than the
# source, so a click ALWAYS opens the latest version — and says so while it does.
#
# LATEST VERSION ALWAYS. A click means "open the current app", so this launcher
# does not open a stale bundle: if any source is newer than the built binary it
# rebuilds first, then opens. (It used to only print a warning and open the stale
# bundle — the trap that let a months-old binary answer as if it were current.)
# An unchanged bundle opens straight away; only a real source edit pays a rebuild.
#
# It launches the RELEASE app — the real, self-contained bundle with its own
# kernel sidecar — not `cargo tauri dev`. For development (hot reload, the vite
# dev server) run `npm run tauri dev` in frontend/ instead; see frontend/README.md.

set -e

ROOT=$(cd "$(dirname "$0")" && pwd)
APP="$ROOT/frontend/src-tauri/target/release/bundle/macos/Mathilda.app"

if [ "$(uname -s)" != "Darwin" ]; then
    echo "This launcher is macOS-only (it uses \`open\`)."
    echo "On Linux, build with:  cd '$ROOT/frontend' && npm run tauri build"
    echo "then run the binary in frontend/src-tauri/target/release/."
    exit 1
fi

# Finder starts Terminal with the login shell's PATH, which usually has node and
# cargo on it — but not always, and a missing cargo halfway through a build is a
# confusing failure. Put the usual homes on PATH explicitly.
PATH="$HOME/.cargo/bin:/opt/homebrew/bin:/usr/local/bin:$PATH"
export PATH

# Decide whether the bundle must be (re)built before opening: it is missing, or a
# source file is newer than the built binary. Either way the rebuild happens here,
# so what opens is always the latest version rather than a stale bundle.
#
# The comparison is against the executable INSIDE the bundle, not the bundle
# directory: a directory's mtime tracks only its immediate entries, so touching a
# nested file during a rebuild need not move Mathilda.app itself. The executable is
# NOT named after the bundle -- productName is "Mathilda" but the binary is the
# crate, mathilda-notebook -- so it is read from Info.plist rather than guessed.
#
# The file set spans everything a build turns into behaviour on screen: the C
# kernel and its `.m` bootstrap (the sidecar), the Svelte/TS/CSS front end (CSS is
# where the syntax-scheme colour ramps live, so a colour-only edit must count), and
# the Rust/Tauri shell (where the native menus are compiled in).
NEED_BUILD=""
REASON=""
if [ ! -d "$APP" ]; then
    NEED_BUILD=1
    REASON="not built yet"
else
    STAMP=""
    EXE=$(/usr/libexec/PlistBuddy -c 'Print :CFBundleExecutable' \
            "$APP/Contents/Info.plist" 2>/dev/null || true)
    [ -n "$EXE" ] && [ -f "$APP/Contents/MacOS/$EXE" ] && STAMP="$APP/Contents/MacOS/$EXE"
    [ -n "$STAMP" ] || STAMP="$APP"
    NEWER=$(find "$ROOT/src" "$ROOT/frontend/src" "$ROOT/frontend/src-tauri/src" \
                 -type f \( -name '*.c' -o -name '*.h' -o -name '*.m' \
                            -o -name '*.rs' -o -name '*.ts' -o -name '*.svelte' \
                            -o -name '*.css' \) \
                 -newer "$STAMP" -print 2>/dev/null | head -1)
    if [ -n "$NEWER" ]; then
        NEED_BUILD=1
        REASON="sources changed since it was built (e.g. ${NEWER#"$ROOT"/})"
    fi
fi

if [ -n "$NEED_BUILD" ]; then
    echo "Building Mathilda.app — $REASON."
    echo "The first build takes a few minutes (the kernel, then the Rust release"
    echo "build); an incremental rebuild after a small edit is faster. Only a"
    echo "source change triggers this, so an unchanged app opens straight away."
    echo
    # The ONE build path, shared with `npm run build:app`: kernel sidecar, JS
    # dependencies, release bundle. Deliberately not reimplemented here — a
    # launcher that knew how to build the app would drift from the build itself.
    "$ROOT/frontend/build-app.sh"
    [ -d "$APP" ] || { echo "Build finished but $APP is missing."; exit 1; }
    echo
fi

echo "Opening Mathilda Notebook…"
open "$APP"
