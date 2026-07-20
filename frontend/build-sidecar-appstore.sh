#!/bin/sh
# build-sidecar-appstore.sh
#
# Builds a SELF-CONTAINED Mac App Store (MAS) sidecar for the Mathilda
# notebook. The standard build-sidecar.sh drops a binary that links against
# Homebrew dylibs under /opt/homebrew — those paths do not exist on end-user
# machines and are disallowed under the App Sandbox that MAS requires.
#
# This script:
#   1. Builds the C engine (make -j at repo root).
#   2. Copies the sidecar to src-tauri/binaries/mathilda-<triple>.
#   3. Copies every required Homebrew dylib (gmp, mpfr, pcre2, raylib) AND any
#      transitive Homebrew deps into src-tauri/binaries/libs/.
#   4. Rewrites, with install_name_tool:
#        (a) each dylib's own LC_ID_DYLIB   -> @rpath/<lib>
#        (b) the sidecar's LC_LOAD_DYLIB    -> @rpath/<lib>
#        (c) inter-dylib refs (mpfr->gmp)   -> @loader_path/<lib>
#      and adds two LC_RPATHs to the sidecar so @rpath resolves both when the
#      dylibs sit in Contents/Frameworks (MAS bundle) and in ./libs next to the
#      sidecar (dev / manual runs).
#   5. Ad-hoc re-signs (codesign -f -s -) every dylib and the sidecar, because
#      install_name_tool invalidates any existing signature. (The real MAS
#      build re-signs everything with the Apple Distribution identity; ad-hoc
#      signing here just keeps the artifacts loadable for local verification.)
#
# NOTE: arm64-only for now. Universal (arm64 + x86_64) requires x86_64 builds of
# GMP/MPFR/pcre2, which is future work (see docs/appstore.md).
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARIES_DIR="$SCRIPT_DIR/src-tauri/binaries"
LIBS_DIR="$BINARIES_DIR/libs"

# --- target triple -----------------------------------------------------------
TARGET=$(rustc -vV | grep host | awk '{print $2}')
if [ -z "$TARGET" ]; then
    echo "ERROR: could not determine Rust target triple" >&2
    exit 1
fi
case "$TARGET" in
    aarch64-apple-darwin) ;;
    *)
        echo "ERROR: this script is arm64-only for now (got $TARGET)." >&2
        echo "       Universal builds are documented as future work in docs/appstore.md." >&2
        exit 1
        ;;
esac

SIDECAR="$BINARIES_DIR/mathilda-$TARGET"

echo "==> Building Mathilda for target: $TARGET"
make -C "$REPO_ROOT" USE_ECM=0 -j4

mkdir -p "$BINARIES_DIR" "$LIBS_DIR"
cp "$REPO_ROOT/Mathilda" "$SIDECAR"
echo "==> Sidecar copied: $SIDECAR"

# --- collect Homebrew dylibs, recursively (transitive deps) ------------------
# Walk the dependency graph starting from the sidecar. Any load command whose
# path begins with /opt/homebrew is a dylib we must bundle. libSystem, libedit,
# /usr/lib/* and /System/Library/* frameworks stay as OS-provided references.
collect() {
    # $1 = mach-o file to inspect
    otool -L "$1" | tail -n +2 | awk '{print $1}' | while read -r dep; do
        case "$dep" in
            /opt/homebrew/*)
                base=$(basename "$dep")
                if [ ! -f "$LIBS_DIR/$base" ]; then
                    echo "    bundling $base  (from $dep)"
                    cp "$dep" "$LIBS_DIR/$base"
                    chmod u+w "$LIBS_DIR/$base"
                    # Recurse into the freshly copied dylib to catch its own
                    # Homebrew dependencies.
                    collect "$LIBS_DIR/$base"
                fi
                ;;
        esac
    done
}

echo "==> Collecting Homebrew dylibs (recursively)"
collect "$SIDECAR"

echo "==> Bundled dylibs:"
ls -1 "$LIBS_DIR"

# --- rewrite install names ---------------------------------------------------
# (a) each dylib's own LC_ID_DYLIB -> @rpath/<lib>
echo "==> Rewriting dylib IDs (LC_ID_DYLIB -> @rpath/...)"
for lib in "$LIBS_DIR"/*.dylib; do
    base=$(basename "$lib")
    install_name_tool -id "@rpath/$base" "$lib"
done

# (c) inter-dylib references -> @loader_path/<lib>
# For every bundled dylib, any load command still pointing at /opt/homebrew must
# point at a sibling in the same directory. @loader_path resolves relative to
# the dylib doing the loading, so it works wherever the libs/ dir is placed.
echo "==> Rewriting inter-dylib references (-> @loader_path/...)"
for lib in "$LIBS_DIR"/*.dylib; do
    otool -L "$lib" | tail -n +2 | awk '{print $1}' | while read -r dep; do
        case "$dep" in
            /opt/homebrew/*)
                base=$(basename "$dep")
                install_name_tool -change "$dep" "@loader_path/$base" "$lib"
                ;;
        esac
    done
done

# (b) sidecar's LC_LOAD_DYLIB entries -> @rpath/<lib>, plus LC_RPATHs
echo "==> Rewriting sidecar load commands (-> @rpath/...)"
otool -L "$SIDECAR" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    case "$dep" in
        /opt/homebrew/*)
            base=$(basename "$dep")
            install_name_tool -change "$dep" "@rpath/$base" "$SIDECAR"
            ;;
    esac
done

# Add rpaths so @rpath resolves in both deployment layouts:
#   - MAS bundle:   sidecar in Contents/MacOS, dylibs in Contents/Frameworks
#   - dev / manual: dylibs in ./libs next to the sidecar
add_rpath() {
    # add $2 to $1 only if not already present (install_name_tool errors on dupes)
    if ! otool -l "$1" | grep -A2 LC_RPATH | grep -q " $2$"; then
        install_name_tool -add_rpath "$2" "$1"
    fi
}
add_rpath "$SIDECAR" "@executable_path/../Frameworks"
add_rpath "$SIDECAR" "@executable_path/libs"

# --- re-sign (ad-hoc) --------------------------------------------------------
# install_name_tool invalidates any existing code signature. Re-sign ad-hoc so
# the artifacts load for local verification. The production MAS build re-signs
# with "Apple Distribution: ..." (see docs/appstore.md).
echo "==> Ad-hoc re-signing dylibs and sidecar"
for lib in "$LIBS_DIR"/*.dylib; do
    codesign -f -s - "$lib"
done
codesign -f -s - "$SIDECAR"

echo ""
echo "==> DONE. Verifying no /opt/homebrew paths remain:"
echo "--- sidecar ---"
otool -L "$SIDECAR"
for lib in "$LIBS_DIR"/*.dylib; do
    echo "--- $(basename "$lib") ---"
    otool -L "$lib"
done

if otool -L "$SIDECAR" "$LIBS_DIR"/*.dylib | grep -q "/opt/homebrew"; then
    echo "ERROR: /opt/homebrew references still present!" >&2
    exit 1
fi
echo ""
echo "OK: sidecar and all bundled dylibs are free of /opt/homebrew references."
