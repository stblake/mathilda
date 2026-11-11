#!/bin/sh
# lean_leakcheck.sh -- leak gate for the Lean 4 integration (LeanForm, the
# strict translator, every Certify generator, Reduce`SOSCertificate,
# Lean`Oracle, and the Lean runner).
#
# The C test binaries cannot detect a leak on their own, so this is the real
# gate.  It drives, in loops:
#   - LeanForm / ToString[_, LeanForm] and Lean`Internal`ToLean (with refusals),
#   - every certificate generator with Verify -> False (no Lean needed),
#   - Reduce`SOSCertificate (success and $Failed),
#   - Lean`Oracle, all four query kinds plus their failure answers,
# and, when a Lean executable is found ($LeanVersion is a string), also
#   - real Certify runs (True and a refuted False), a batch, cache hits,
#     LeanVerify (accepted, sorry, native_decide) and Lean`Internal`RunProcess.
# It FAILS on any leaked byte.
#
#   macOS : leaks --atExit (with MallocStackLogging=1 for stack traces)
#   Linux : valgrind --error-exitcode
#   neither present: SKIPS, loudly, with exit 0 -- an unavailable checker is not
#   a pass, so it says so rather than printing nothing (docs/adr/0004).
#
# Usage: sh tests/scripts/lean_leakcheck.sh [path-to-Mathilda]
BIN="${1:-./Mathilda}"
[ -x "$BIN" ] || { echo "lean_leakcheck: no Mathilda binary at $BIN — build first"; exit 2; }

WORK="$(mktemp -t leanleak.XXXXXX)"
cat > "$WORK" <<'MSCRIPT'
(* --- C printer and strict translator --- *)
Do[ToString[LeanForm[ForAll[{x, y}, (x + y)^2 == x^2 + 2 x y + y^2]]];
   ToString[{Function[{x}, x + 1], 3/4 a, Log[b], f[-3], "s", 1.5, Pi, Sqrt[x]}, LeanForm];
   ToString[LeanForm[Implies[a < b <= c && Not[d], Exists[z, z > 0, z^2 == 2]]]];
   Lean`Internal`ToLean[x^2/3 + 1/2, "Rat", True, {}];
   Quiet[Lean`Internal`ToLean[x/2, "Int", True, {}]];
   Quiet[Lean`Internal`ToLean[Log[x], "Rat", True, {}]];
   Quiet[Lean`Internal`ToLean[1.5 x, "Rat", True, {}]];
   Quiet[Lean`Internal`ToLean[g[x], "Rat", True, {}]];
   Quiet[Lean`Internal`ToLean[x^(1/3), "Rat", True, {}]];
   Quiet[Lean`Internal`ToLean[#^2 &, "Rat", True, {}]], {50}];

(* --- every generator, no Lean run --- *)
gens = Hold[
  Certify[Factor[x^10 - 1], Verify -> False],
  Certify[Expand[(x + y)^5], Verify -> False],
  Certify[Together[1/x + 1/(x + 1)], Verify -> False],
  Certify[(x + 1)^2 == x^2 + 1, Verify -> False],
  Certify[Integrate[(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4), x], Verify -> False],
  Certify[Integrate[1/(x^2 + 2), x], Verify -> False],
  Certify[Integrate[1/(x^2 - 1), x], Verify -> False],
  Certify[Integrate[(1 + x)/(1 + x^2), x], Verify -> False, Tier -> 1, LeanProject -> "/nonexistent"],
  Certify[PolynomialExtendedGCD[x^4 - 1, x^3 + 2 x^2 + 2 x + 1, x], Verify -> False],
  Certify[PolynomialReduce[x^3 y - x y^2, {x^2 - y, x y - 1}, {x, y}], Verify -> False],
  Certify[Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1], Verify -> False],
  Certify[Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x], Verify -> False],
  Certify[Solve[{x + y == 3, x - y == 1}, {x, y}], Verify -> False],
  Certify[Inverse[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}], Verify -> False],
  Certify[Det[{{2, -1, 0}, {1, 4, 2}, {0, 5, -2}}], Verify -> False],
  Certify[LinearSolve[{{2, 1}, {1, 3}}, {1, 2}], Verify -> False],
  Certify[NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}], Verify -> False],
  Certify[ForAll[{x, y}, x^2 + y^2 >= 2 x y], Verify -> False],
  Certify[ForAll[{x, y}, x^4 - 2 x^2 + 1 + y^2 >= 0], Verify -> False],
  Certify[Implies[x >= 0, x^3 >= 0], Verify -> False],
  Certify[Minimize[{x^2 + y^2, x + y >= 1}, {x, y}], Verify -> False],
  Certify[Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}], Verify -> False],
  Certify[FactorInteger[2^64 + 1], Verify -> False],
  Certify[PrimeQ[274177], Verify -> False],
  Certify[PrimeQ[274179], Verify -> False],
  Certify[PrimeQ[2], Verify -> False],
  Certify[Sum[k^3, {k, 1, n}], Verify -> False],
  Certify[Sum[k 2^k, {k, 0, n}], Verify -> False],
  Certify[{Factor[x^6 - 1], Det[{{1, 2}, {3, 4}}]}, Verify -> False],
  Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, Engine -> "Mathilda"],
  Certify[Factor[x^2 - 1], Engine -> "Bogus"],
  Certify[Sin[x] == Cos[x], Verify -> False],
  ProofObject[Factor[x^4 - 1], Verify -> False]];
Do[Quiet[ReleaseHold[Map[Function[c, With[{o = c}, If[Head[o] === ProofObject,
     {o["LeanSource"], o["Properties"], ToString[o]}, o]]], gens]]], {5}];

(* --- the SOS engine --- *)
Do[Reduce`SOSCertificate[x^4 + y^4 - x^2 y^2 + 1, {}, {x, y}];
   Reduce`SOSCertificate[x^4 - 2 x^2 + 1 + y^2, {}, {x, y}];
   Reduce`SOSCertificate[x^2 - 1, {}, {x}], {3}];

(* --- Lean`Oracle: all four kinds and their failures --- *)
Do[Lean`Oracle["LinearCombination", {2 v2^2, 1}, {{v1^2 + v2^2, 1}, {v1, v2}}, {v1, v2}];
   Lean`Oracle["LinearCombination", {v1, 2}, {{v1, 1}}, {v1}];
   Lean`Oracle["SOS", {2 v1 v2, v1^2 + v2^2}, {}, {v1, v2}];
   Lean`Oracle["SOS", {0, v1^3}, {v1}, {v1}];
   Lean`Oracle["SOS", {v1, 1}, {}, {v1}];
   Lean`Oracle["Pratt", 274177];
   Lean`Oracle["Pratt", 274179];
   Lean`Oracle["Eval", "Factor[x^10 - 1]"];
   Lean`Oracle["Eval", "Factor["];
   Lean`Oracle["Nope", 1], {10}];

(* --- real Lean runs, when Lean is installed --- *)
If[StringQ[$LeanVersion],
  Quiet[
    Certify[Integrate[1/(x^2 + 2), x]];
    Certify[ForAll[{x, y}, x^2 + y^2 >= 2 x y]];
    Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x];
    Certify[{Factor[x^10 - 1], Inverse[{{2, 1}, {1, 1}}], FactorInteger[2^64 + 1],
             (x + 1)^2 == x^2 + 1}];
    Do[Certify[Factor[x^10 - 1]], {5}];
    LeanVerify["theorem ok : 2 + 2 = 4 := by decide"];
    LeanVerify["theorem s : 2 + 2 = 5 := by sorry"];
    LeanVerify["theorem nd : 2 + 2 = 4 := by native_decide"];
    Lean`Internal`RunProcess[{"sh", "-c", "echo out; echo err 1>&2; exit 3"}, None]];
  Print["lean_leakcheck: real Lean runs included"],
  Print["lean_leakcheck: no Lean found; real Lean runs skipped"]];
MSCRIPT

if command -v leaks >/dev/null 2>&1; then
    RAW="$(MallocStackLogging=1 leaks --atExit -- "$BIN" -file "$WORK" 2>/dev/null)"
    rm -f "$WORK"
    echo "$RAW" | grep '^lean_leakcheck:'
    OUT="$(echo "$RAW" | grep 'total leaked bytes' | tail -1)"
    echo "lean_leakcheck: $OUT"
    case "$OUT" in
        *" 0 leaks for 0 total leaked bytes"*) echo "lean_leakcheck: PASS"; exit 0 ;;
        "") echo "lean_leakcheck: FAIL — no verdict line from leaks"; exit 1 ;;
        *) echo "$RAW" | sed -n '/^leaks Report/,$p' | head -80
           echo "lean_leakcheck: FAIL — leaked memory above"; exit 1 ;;
    esac
elif command -v valgrind >/dev/null 2>&1; then
    LOG="$(mktemp -t leanleak-vg.XXXXXX)"
    valgrind --leak-check=full --errors-for-leak-kinds=definite \
             --error-exitcode=1 "$BIN" -file "$WORK" >/dev/null 2>"$LOG"
    rc=$?
    rm -f "$WORK"
    [ $rc -eq 0 ] && { rm -f "$LOG"; echo "lean_leakcheck: PASS (valgrind)"; exit 0; }
    echo "lean_leakcheck: FAIL (valgrind) — see $LOG"; exit 1
else
    rm -f "$WORK"
    echo "lean_leakcheck: SKIPPED — neither 'leaks' nor 'valgrind' found."
    echo "lean_leakcheck: this is NOT a pass; the leak gate did not run."
    exit 0
fi
