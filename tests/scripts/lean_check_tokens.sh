#!/bin/sh
# lean_check_tokens.sh -- keep LeanVerify's command-keyword list in sync with
# Lean's own parser tables.
#
# Usage: sh tests/scripts/lean_check_tokens.sh [--mathlib] [--words] <Mathilda binary>
#
# Runs tools/lean_command_tokens.lean (a trusted script of this repository, not
# a LeanVerify source) to print every command-leading token of Lean's parser
# and the token table, then checks with tests/scripts/lean_check_tokens.m that
# Lean`Internal`CheckSource refuses each of them (see there).
#
#   default    `import Lean` (a superset of the core Init that a LeanVerify
#              source sees), with the newest Lean found: MATHILDA_LEAN, else
#              the newest ~/.elan/toolchains/*/bin/lean, else lean on PATH.
#   --mathlib  `import Mathlib`, elaborated with `lake env` in the project of
#              MATHILDA_LEAN_PROJECT (read only: the script is copied to a
#              private temp directory, nothing is written to the project).
#   --words    print the words for src/lean_check_tokens.h instead of checking.
#
# Exit 77 (ctest SKIP) when no Lean (or, with --mathlib, no project) is found.

MATHLIB=0
WORDS=0
while [ $# -gt 1 ]; do
    case $1 in
        --mathlib) MATHLIB=1 ;;
        --words) WORDS=1 ;;
        *) echo "usage: $0 [--mathlib] [--words] <Mathilda binary>"; exit 2 ;;
    esac
    shift
done
BIN=$1
HERE=$(cd "$(dirname "$0")" && pwd)
TOOL=$HERE/../../tools/lean_command_tokens.lean

if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "SKIP: Mathilda binary '$BIN' not found (build it with make)"
    exit 77
fi
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")

LEAN=${MATHILDA_LEAN:-}
if [ -z "$LEAN" ]; then
    for d in "$HOME"/.elan/toolchains/*/bin/lean; do
        [ -x "$d" ] && LEAN=$d
    done
fi
[ -z "$LEAN" ] && LEAN=$(command -v lean 2>/dev/null)
if [ -z "$LEAN" ] || [ ! -x "$LEAN" ]; then
    echo "SKIP: no Lean 4 found (set MATHILDA_LEAN)"
    exit 77
fi

TMP=$(mktemp -d "${TMPDIR:-/tmp}/mathilda-tokens-XXXXXX") || exit 1
trap 'rm -rf "$TMP"' EXIT INT TERM

if [ "$MATHLIB" = 1 ]; then
    PROJ=${MATHILDA_LEAN_PROJECT:-}
    if [ -z "$PROJ" ] || [ ! -d "$PROJ" ]; then
        echo "SKIP: MATHILDA_LEAN_PROJECT is not set to a Lake project with Mathlib"
        exit 77
    fi
    LAKE=$(dirname "$LEAN")/lake
    [ -x "$LAKE" ] || { echo "SKIP: no lake beside $LEAN"; exit 77; }
    sed 's/^import Lean$/import Mathlib/' "$TOOL" > "$TMP/tokens.lean"
    (cd "$PROJ" && "$LAKE" env "$LEAN" "$TMP/tokens.lean") > "$TMP/tokens.txt" 2> "$TMP/err.txt"
else
    cp "$TOOL" "$TMP/tokens.lean"
    "$LEAN" "$TMP/tokens.lean" > "$TMP/tokens.txt" 2> "$TMP/err.txt"
fi
st=$?
if [ $st -ne 0 ] || ! grep -q '^lead ' "$TMP/tokens.txt"; then
    echo "FAIL: tools/lean_command_tokens.lean did not run (exit $st)"
    head -20 "$TMP/tokens.txt" "$TMP/err.txt"
    exit 1
fi
echo "lean: $LEAN"
MATHILDA_TOKENS_FILE=$TMP/tokens.txt MATHILDA_TOKENS_MATHLIB=$MATHLIB MATHILDA_TOKENS_WORDS=$WORDS \
    "$BIN" -file "$HERE/lean_check_tokens.m"
