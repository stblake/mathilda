#!/bin/sh
# field_kernels_ab.sh -- A/B byte-identity gate for the native field-coefficient
# kernels (CoefficientRules / MonomialList, and — once landed — Together/Cancel).
#
# The field fast paths (flint_field_monomials, flint_field_together/_cancel) are
# an ACCELERATION of the generic evaluator path: they must produce byte-identical
# output.  This runs a battery of AlgebraicNumber[theta,{..}]-coefficient inputs
# (the ParallelMixedTower assembly's representation) twice — once with the field
# path ON (default) and once with MATHILDA_NO_FIELD_KERNELS=1 forcing the generic
# path — and FAILS on any difference.  It also confirms (via the field diag) that
# the fast path actually fired, so an all-decline regression cannot masquerade as
# a pass.
#
# Usage: sh tests/scripts/field_kernels_ab.sh [path-to-Mathilda]
set -e
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "field_kernels_ab: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t fieldab.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
(* Q(sqrt 2) *)
Print[FullForm[CoefficientRules[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + 3 y - 1)^3], {x, y}]]];
Print[FullForm[MonomialList[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + 3 y - 1)^3], {x, y}]]];
(* reduction modulo the minimal polynomial *)
Print[FullForm[CoefficientRules[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + AlgebraicNumber[Sqrt[2],{1,0}])^4], {x}]]];
(* Gaussian field Q(i), including a bare Complex mixed with AlgebraicNumber[I,..] *)
Print[FullForm[CoefficientRules[Expand[(AlgebraicNumber[I,{0,1}] x + (2 - I) y)^2 + I], {x, y}]]];
(* degree-3 field Q(2^(1/3)) via a Root generator *)
Print[FullForm[CoefficientRules[Expand[(AlgebraicNumber[2^(1/3),{0,1}] x + 1)^3], {x}]]];
(* compositum Q(sqrt2, sqrt3) built by ToNumberField, then a polynomial over it *)
nf = ToNumberField[{Sqrt[2], Sqrt[3]}];
p = Expand[(nf[[1]] x + nf[[2]] y)^3 + nf[[1]] nf[[2]] x y];
Print[FullForm[CoefficientRules[p, {x, y}]]];
Print[FullForm[MonomialList[p, {x, y}]]];
(* explicit monomial order *)
Print[FullForm[CoefficientRules[Expand[(AlgebraicNumber[Sqrt[2],{0,1}] x + 3 y - 1)^3], {x, y}, "DegreeReverseLexicographic"]]];
MSCRIPT

ON="$(MATHILDA_FIELD_DIAG=1 "$BIN" -file "$WORK" 2>/tmp/fieldab_diag.txt)"
OFF="$(MATHILDA_NO_FIELD_KERNELS=1 "$BIN" -file "$WORK" 2>/dev/null)"
rm -f "$WORK"

FIRED="$(grep -c 'read-off' /tmp/fieldab_diag.txt 2>/dev/null || echo 0)"

if [ "$ON" != "$OFF" ]; then
    echo "field_kernels_ab: FAIL — CoefficientRules field-path output differs from generic path"
    printf '%s\n' "$ON" > /tmp/fieldab_on.txt
    printf '%s\n' "$OFF" > /tmp/fieldab_off.txt
    diff /tmp/fieldab_off.txt /tmp/fieldab_on.txt | head -40
    exit 1
fi
if [ "$FIRED" -lt 1 ]; then
    echo "field_kernels_ab: FAIL — the field CoefficientRules fast path never fired (all declines); outputs match only vacuously"
    exit 1
fi

# --- Together/Cancel: byte-identity need not hold (canonical form may differ),
#     so gate on VALUE-equality at a generic rational point + that the field path
#     actually fired. ---
TC="$(mktemp -t fieldabtc.XXXXXX)"
cat > "$TC" <<'MSCRIPT'
th = Root[#^4 - #^2 - 1 &, 1];
a = AlgebraicNumber[th, {0,1,0,0}]; b = AlgebraicNumber[th, {0,0,1,0}];
Print[Abs[N[(Together[a/(x - a) + b/(x^2 - a)] - (a/(x - a) + b/(x^2 - a))) /. x -> 7/13, 30]] < 10^-25];
Print[Abs[N[(Together[a/(x - a) + 1/(x + a) - a/x] - (a/(x - a) + 1/(x + a) - a/x)) /. x -> 5/17, 30]] < 10^-25];
i = AlgebraicNumber[I, {0,1}];
Print[Abs[N[(Together[i/(x - i) + 3/x] - (i/(x - i) + 3/x)) /. x -> 9/11, 30]] < 10^-25];
Print[Simplify[Cancel[Expand[(x - a)(x - b)]/(x - a)] - (x - b)] === 0];
MSCRIPT
TCOUT="$(MATHILDA_FIELD_DIAG=1 "$BIN" -file "$TC" 2>/tmp/fieldab_tcdiag.txt)"
rm -f "$TC"
TCFIRED="$(grep -c 'Together/Cancel candidate' /tmp/fieldab_tcdiag.txt 2>/dev/null || echo 0)"
NTRUE="$(printf '%s\n' "$TCOUT" | grep -c '^True$')"
if [ "$NTRUE" -ne 4 ]; then
    echo "field_kernels_ab: FAIL — a field Together/Cancel result is not value-equal to its input:"
    printf '%s\n' "$TCOUT"
    exit 1
fi
if [ "$TCFIRED" -lt 1 ]; then
    echo "field_kernels_ab: FAIL — the field Together/Cancel fast path never fired"
    exit 1
fi

echo "field_kernels_ab: PASS — CoefficientRules byte-identical (fired $FIRED); Together/Cancel value-equal (fired $TCFIRED)"
exit 0
