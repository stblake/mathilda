#!/bin/sh
# field_kernels_leakcheck.sh -- leak gate for the native field-coefficient
# kernels (flint_field_monomials for CoefficientRules/MonomialList, and — once
# landed — flint_field_together/_cancel).
#
# The C test binary cannot detect a leak on its own, so this is the real gate.
# Drives the field kernels in Do[...,{200}] loops over every code path: the
# successful read-off (Q(sqrt2), Q(i), degree-3 Root field, compositum, with and
# without mod-M reduction and rational collapse), MonomialList, and the DECLINE
# paths (stray parameter symbol, non-field/pure-rational input) which must free
# their scratch and fall through cleanly.  FAILS on any leaked byte.
#
#   macOS : leaks --atExit
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 (docs/adr/0004).
#
# Usage: sh tests/scripts/field_kernels_leakcheck.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "field_kernels_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t fieldleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
(* successful field read-off *)
Do[CoefficientRules[AlgebraicNumber[Sqrt[2],{0,1}] x^2 + 3, {x}], {200}];
Do[CoefficientRules[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + AlgebraicNumber[Sqrt[2],{1,0}])^3], {x}], {200}];
Do[CoefficientRules[AlgebraicNumber[I,{0,1}] x + AlgebraicNumber[I,{2,0}], {x}], {200}];
Do[CoefficientRules[AlgebraicNumber[2^(1/3),{0,1}] x^3 + AlgebraicNumber[2^(1/3),{0,0,1}], {x}], {200}];
Do[CoefficientRules[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + 3 y - 1)^3], {x, y}], {200}];
Do[MonomialList[AlgebraicNumber[Sqrt[2],{0,1}] x^2 + 3, {x}], {200}];
(* rational-collapse group *)
Do[CoefficientRules[AlgebraicNumber[Sqrt[2],{0,1}] x^2 + AlgebraicNumber[Sqrt[2],{5,0}], {x}], {200}];
(* decline paths: stray symbol, pure-rational (no field), non-symbol var *)
Do[CoefficientRules[AlgebraicNumber[Sqrt[2],{0,1}] x^2 + a, {x}], {200}];
Do[CoefficientRules[3 x^2 + 2 x + 1, {x}], {200}];
Do[CoefficientRules[AlgebraicNumber[Sqrt[2],{0,1}] x[1]^2 + 1, {x[1]}], {200}];
(* field Together / Cancel: sqrt2, Q(i), degree-4 Root field, cancellation *)
Do[Together[AlgebraicNumber[Sqrt[2],{0,1}]/(x - AlgebraicNumber[Sqrt[2],{0,1}]) + 1/x], {200}];
Do[Together[AlgebraicNumber[I,{0,1}]/(x - AlgebraicNumber[I,{0,1}]) + 1/(x + AlgebraicNumber[I,{0,1}])], {200}];
Do[Together[AlgebraicNumber[Root[#^4-#^2-1&,1],{0,1,0,0}]/(x - AlgebraicNumber[Root[#^4-#^2-1&,1],{0,1,0,0}]) + 1/x], {200}];
Do[Cancel[Expand[(x - AlgebraicNumber[Root[#^4-#^2-1&,1],{0,1,0,0}])(x - 1)]/(x - AlgebraicNumber[Root[#^4-#^2-1&,1],{0,1,0,0}])], {200}];
(* field Together decline path: multivariate (two symbols) falls to generic *)
Do[Together[AlgebraicNumber[Sqrt[2],{0,1}]/(x y) + 1/x], {200}];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    OUT="$(leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null | grep 'total leaked bytes' | tail -1)"
    rm -f "$WORK"
    echo "field_kernels_leakcheck: $OUT"
    case "$OUT" in
        *"0 leaks for 0 total leaked bytes"*) echo "field_kernels_leakcheck: PASS"; exit 0 ;;
        "") echo "field_kernels_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "field_kernels_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>/tmp/fieldleak.valgrind
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { echo "field_kernels_leakcheck: PASS (valgrind)"; exit 0; }
    echo "field_kernels_leakcheck: FAIL (valgrind) — see /tmp/fieldleak.valgrind"; exit 1
else
    rm -f "$WORK"
    echo "field_kernels_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "field_kernels_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi
