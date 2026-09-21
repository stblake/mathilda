#!/bin/sh
# algebraicnumberpolynomial_leakcheck.sh -- leak gate for AlgebraicNumberPolynomial.
#
# The C test binary cannot detect a leak on its own, so this is the real gate.
# Drives AlgebraicNumberPolynomial over its code paths -- integer/bigint/rational
# passthrough, the polynomial build (basic, rational coeff, zero-drop, higher
# degree), Listable threading, and the error (naobj) / wrong-argc declines --
# under a leak checker and FAILS on any leaked byte.
#
#   macOS : leaks --atExit
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 -- an unavailable checker is not
#   a pass, so it says so rather than printing nothing (docs/adr/0004).
#
# Usage: bash tests/scripts/algebraicnumberpolynomial_leakcheck.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "algebraicnumberpolynomial_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t anpleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
Do[AlgebraicNumberPolynomial[2, x], {200}];
Do[AlgebraicNumberPolynomial[-3, x], {200}];
Do[AlgebraicNumberPolynomial[1/2, x], {200}];
Do[AlgebraicNumberPolynomial[100!, x], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 1/2}], x], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[2^(1/3), {1, 0, 3}], x], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,2,3,4}], x], {200}];
Do[AlgebraicNumberPolynomial[{2, AlgebraicNumber[Sqrt[2], {1, 2}]}, x], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], {x, y}], {200}];
Do[AlgebraicNumberPolynomial[Sqrt[2], x], {200}];
Do[AlgebraicNumberPolynomial[y^2, z], {200}];
Do[AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}]], {200}];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    OUT="$(leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null | grep 'total leaked bytes' | tail -1)"
    rm -f "$WORK"
    echo "algebraicnumberpolynomial_leakcheck: $OUT"
    case "$OUT" in
        *"0 leaks for 0 total leaked bytes"*) echo "algebraicnumberpolynomial_leakcheck: PASS"; exit 0 ;;
        "") echo "algebraicnumberpolynomial_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "algebraicnumberpolynomial_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>/tmp/anpleak.valgrind
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { echo "algebraicnumberpolynomial_leakcheck: PASS (valgrind)"; exit 0; }
    echo "algebraicnumberpolynomial_leakcheck: FAIL (valgrind) — see /tmp/anpleak.valgrind"; exit 1
else
    rm -f "$WORK"
    echo "algebraicnumberpolynomial_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "algebraicnumberpolynomial_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi
