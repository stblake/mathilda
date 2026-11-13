#!/bin/sh
# lean_tactics.sh -- test the Lean -> Mathilda package in lean/.
#
# Usage: sh tests/scripts/lean_tactics.sh <path to the Mathilda binary>
#
# Builds the core-only Lake package lean/ and elaborates
#   lean/test/Tactics.lean   (MATHILDA_BIN = the real binary: every tactic closes
#                             its goal, only standard axioms, false goals fail),
#   lean/test/Untrusted.lean (MATHILDA_BIN = lean/test/bad_oracle.sh, a fake
#                             oracle with wrong, malformed and hostile answers:
#                             Lean rejects all, and runs none of them),
#   lean/test/Missing.lean   (MATHILDA_BIN = a missing path: a clear error).
# Each file pins its expected messages with #guard_msgs, so `lean` exiting 0
# means every check passed.
#
# lake is MATHILDA_LAKE, else the elan toolchain named in lean/lean-toolchain,
# else the newest ~/.elan/toolchains/*/bin/lake.  Exit 77 (ctest SKIP) when no
# lake or no built Mathilda binary is found.  Nothing is downloaded.

BIN=$1
HERE=$(cd "$(dirname "$0")" && pwd)
PKG=$(cd "$HERE/../../lean" && pwd)

if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "SKIP: Mathilda binary '$BIN' not found (build it with make)"
    exit 77
fi
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")

LAKE=${MATHILDA_LAKE:-}
if [ -z "$LAKE" ]; then
    TC=$(sed -e 's|/|--|' -e 's|:|---|' "$PKG/lean-toolchain" | tr -d '\n')
    if [ -x "$HOME/.elan/toolchains/$TC/bin/lake" ]; then
        LAKE=$HOME/.elan/toolchains/$TC/bin/lake
    else
        for d in "$HOME"/.elan/toolchains/*/bin/lake; do
            [ -x "$d" ] && LAKE=$d
        done
    fi
fi
if [ -z "$LAKE" ] || [ ! -x "$LAKE" ]; then
    echo "SKIP: no lake found (set MATHILDA_LAKE or install the toolchain in lean/lean-toolchain)"
    exit 77
fi
echo "lake: $LAKE"
echo "Mathilda: $BIN"

cd "$PKG" || exit 1
"$LAKE" build || { echo "FAIL: lake build"; exit 1; }

fail=0
run() {
    # run <label> <MATHILDA_BIN> <file>
    t0=$(date +%s)
    if MATHILDA_BIN=$2 "$LAKE" env lean "$3"; then
        echo "PASS: $1 ($3, $(( $(date +%s) - t0 )) s)"
    else
        echo "FAIL: $1 ($3)"
        fail=1
    fi
}
run "tactics close their goals" "$BIN" test/Tactics.lean
# Untrusted.lean checks that a hostile answer's file write (relative to the
# package directory, the working directory here) never happened.
rm -f mathilda_untrusted_pwned
run "wrong certificates are rejected" "$PKG/test/bad_oracle.sh" test/Untrusted.lean
run "missing binary is a clear error" /nonexistent/Mathilda test/Missing.lean
exit $fail
