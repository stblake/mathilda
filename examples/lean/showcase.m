(* ====================================================================== *)
(* examples/lean/showcase.m -- Mathilda <-> Lean 4, end to end.            *)
(*                                                                          *)
(*   ./Mathilda -file examples/lean/showcase.m                              *)
(*                                                                          *)
(* Mathilda computes an answer and a certificate; the Lean 4 kernel judges  *)
(* it.  Every row prints the Mathilda input, the result, Lean's verdict,    *)
(* the tier (0 = core Lean, 1 = Mathlib), the Lean theorem (truncated),     *)
(* how many steps Lean did NOT check ("TrustedSteps"), and the timings      *)
(* (Lean's own wall time and the whole Certify call).                       *)
(*                                                                          *)
(* Needs a Lean 4 toolchain (v4.35 or newer; found through $LeanExecutable, *)
(* MATHILDA_LEAN or ~/.elan).  Row 2 (Mathlib HasDerivAt) runs only when     *)
(* MATHILDA_LEAN_PROJECT names a Lake project with Mathlib built; row 20     *)
(* (Lean -> Mathilda tactics) only when `lake` is found.  Both say SKIPPED   *)
(* otherwise.  Rows 17 and 18 are False by design.  The summary's Exact     *)
(* column re-runs each row with Engine -> "Exact" (Mathilda's exact re-check *)
(* of the certificate, Missing["Exact", b]: b is shown; never a kernel      *)
(* verdict).                                                                *)
(* ====================================================================== *)

$demoRoot = FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
$rows = {};
$width = 150;

clip[s_String] := If[StringLength[s] > $width, StringTake[s, $width - 3] <> "...", s];
clip[e_] := clip[ToString[e, InputForm]];
secs[t_?NumericQ] := ToString[NumberForm[N[t], {6, 2}]] <> " s";
secs[_] := "-";

header[n_, title_] := (Print[""];
  Print["[", n, "] ", title];
  Print["    ", StringJoin[Table["-", {StringLength[title] + 4}]]]);

(* Record a row for the summary table at the end (exact: the row's
   verdict under Engine -> "Exact", "-" when not run). *)
record[n_, title_, verdict_, tier_, t_, trusted_, exact_:"-"] :=
  AppendTo[$rows, {n, title, verdict, tier, t, trusted, exact}];

(* The held call of a row re-run with Engine -> "Exact" (replacing any
   Engine option), for the summary's Exact column. *)
SetAttributes[exactVerdict, HoldAll];
exactVerdict[h_[a_, opts___]] := With[{os = DeleteCases[{opts}, (Engine | "Engine") -> _]},
  exactShow[Quiet[ReleaseHold[Hold[h[a, "Engine" -> "Exact", $opts]] /.
    $opts -> Sequence @@ os]]["Verified"]]];
exactShow[Missing["Exact", b_]] := ToString[b];
exactShow[Missing[w_]] := w;
exactShow[v_] := ToString[v];

(* Print one certified row.  `input` is the literal Mathilda input (held). *)
report[n_, title_, input_String, obj_ProofObject, wall_, expect_:True] := Module[{v},
  v = obj["Verified"];
  header[n, title];
  Print["    In:        ", clip[input]];
  Print["    Result:    ", clip[obj["Result"]]];
  Print["    Verified:  ", v, If[v === expect, "", "   <-- UNEXPECTED (expected " <> ToString[expect] <> ")"],
    "    Tier: ", obj["Tier"], "    Kind: ", obj["Kind"],
    "    TrustedSteps: ", Length[obj["TrustedSteps"]]];
  Print["    Lean:      ", clip[obj["LeanTheorem"]]];
  Print["    Time:      Lean ", secs[obj["Timing"]], If[TrueQ[obj["Cached"]], " (cached)", ""],
    ", whole Certify call ", secs[wall]];
  record[n, title, v, obj["Tier"], wall, obj["TrustedSteps"]];
  obj];

SetAttributes[demo, HoldAll];
demo[n_, title_, call_, expect_:True] := Module[{t, obj},
  t = AbsoluteTiming[call];
  obj = report[n, title, ToString[Unevaluated[call], InputForm], t[[2]], t[[1]], expect];
  $rows = Append[Most[$rows], ReplacePart[Last[$rows], 7 -> exactVerdict[call]]];
  obj];

skipped[n_, title_, why_] := (header[n, title];
  Print["    SKIPPED: ", why];
  record[n, title, "SKIPPED", "-", Missing[], {}]);

Print["Mathilda <-> Lean 4 showcase"];
Print["Checker: ", $LeanVersion];
Print["Mathlib project (Tier 1): ", If[$LeanProject === None, "none", "configured"]];
$t0 = AbsoluteTime[];

(* ---------------------------------------------------------------------- *)
(* 1-3. Antiderivatives: Risch/Rioboo output, certified.                   *)
(* ---------------------------------------------------------------------- *)

r1 = demo[1, "Risch/Rioboo antiderivative, certified in core Lean (Tier 0)",
  Certify[Integrate[(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4), x]]];
Print["    Proved:    the integrand equals the formal derivative of the result, cleared of denominators"];

If[$LeanProject === None,
  skipped[2, "Mathlib Tier 1: HasDerivAt (fun x => F) f x over the reals",
    "set MATHILDA_LEAN_PROJECT to a Lake project with Mathlib built"],
  demo[2, "Mathlib Tier 1: HasDerivAt (fun x => F) f x over the reals",
    Certify[Integrate[(1 + x)/(1 + x^2), x], Tier -> 1]]];

demo[3, "Antiderivative with a surd: s1 stands for Sqrt[2], modulo s1^2 - 2",
  Certify[Integrate[1/(x^2 + 2), x]]];

(* ---------------------------------------------------------------------- *)
(* 4-10. Algebra and linear algebra.                                       *)
(* ---------------------------------------------------------------------- *)

demo[4, "Factorization, through the ProofObject[..., Engine -> \"Lean\"] form",
  ProofObject[Factor[x^10 - 1], Engine -> "Lean"]];

demo[5, "Bezout: s a + t b = g, and g divides a and b (so g is a gcd)",
  Certify[PolynomialExtendedGCD[x^4 - 1, x^3 + 2 x^2 + 2 x + 1, x]]];

demo[6, "Ideal membership: hypotheses imply the goal (Mathilda finds cofactors)",
  Certify[Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1]]];

demo[7, "Solve is complete: p = 0 <-> x = 1 \\/ x = 2 \\/ x = 3",
  Certify[Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]]];

demo[8, "Matrix inverse: A B = I and B A = I, evaluated by the kernel",
  Certify[Inverse[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}]]];

demo[9, "Determinant (Laplace expansion in the kernel)",
  Certify[Det[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, {3, 0, 1, 1}}]]];

demo[10, "NullSpace: A v = 0 for every basis vector",
  Certify[NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]];

(* ---------------------------------------------------------------------- *)
(* 11-13. Inequalities and optimality.                                     *)
(* ---------------------------------------------------------------------- *)

demo[11, "Sum-of-squares inequality (two real zeros, so no interior Gram point)",
  Certify[ForAll[{x, y}, x^4 - 2 x^2 + 1 + y^2 >= 0]]];

r12 = demo[12, "Minimize a QP: KKT multiplier plus an SOS remainder",
  Certify[Minimize[{x^2 + y^2, x + y >= 1}, {x, y}]]];
Print["    Multiplier: ", r12["Certificate"]["Multipliers"]];

r13 = demo[13, "LP optimality via duality: the multipliers are the dual solution",
  Certify[Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}]]];
Print["    Dual multipliers: ", r13["Certificate"]["Multipliers"],
  "    LinearProgram: ", r13["Certificate"]["LinearProgram"]];

(* ---------------------------------------------------------------------- *)
(* 14-16. Number theory and sums.                                          *)
(* ---------------------------------------------------------------------- *)

demo[14, "Pratt primality certificates for both factors of 2^64 + 1",
  Certify[FactorInteger[2^64 + 1]]];

demo[15, "Telescoping sum: Sum k^3 = (n (n+1) / 2)^2 by induction",
  Certify[Sum[k^3, {k, 1, n}]]];

demo[16, "Closed form of Sum k 2^k by induction",
  Certify[Sum[k 2^k, {k, 0, n}]]];

(* ---------------------------------------------------------------------- *)
(* 17. Lean catches the CAS.                                               *)
(* ---------------------------------------------------------------------- *)

header[17, "Lean catches the CAS: a PossibleZeroQ false positive"];
pz = PossibleZeroQ[(1 + x)^20 - Expand[(1 + x)^20] + 10^-25 x];
Print["    In:        PossibleZeroQ[(1 + x)^20 - Expand[(1 + x)^20] + 10^-25 x]"];
Print["    Result:    ", pz, "   (wrong: the expression is the polynomial x/10^25)"];
(* The claim is written out in each call: a symbol holding it has a value,
   which Certify refuses (it certifies claims as written). *)
heur = Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, Engine -> "Mathilda"];
Print["    Engine -> \"Mathilda\":  Verified: ", heur["Verified"], "    Tier: ", heur["Tier"]];
t17 = AbsoluteTiming[Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x,
  Engine -> "Lean"]]];
lean17 = t17[[2]];
Print["    Engine -> \"Lean\":      Verified: ", lean17["Verified"], "    Tier: ", lean17["Tier"],
  If[lean17["Verified"] === False, "", "   <-- UNEXPECTED (expected False)"]];
Print["    Lean:      ", clip[lean17["LeanTheorem"]]];
Print["    Time:      Lean ", secs[lean17["Timing"]], ", whole Certify call ", secs[t17[[1]]]];
record[17, "Lean catches the CAS (the heuristic said " <> ToString[Last[heur["Verified"]]] <> ")",
  lean17["Verified"], lean17["Tier"], t17[[1]], lean17["TrustedSteps"],
  exactVerdict[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x]]];

(* ---------------------------------------------------------------------- *)
(* 18. A tampered certificate is rejected.                                 *)
(* ---------------------------------------------------------------------- *)

header[18, "A deliberately wrong certificate is rejected"];
good = Certify[Factor[x^10 - 1], Verify -> False];
(* LeanVerify is given the source without the trusted footer's #print axioms
   lines (they are Certify's own; a LeanVerify source has no # commands). *)
src = StringReplace[good["LeanSource"], RegularExpression["#print axioms [^\n]*\n"] -> ""];
bad = StringReplace[src, "(1 + x + x ^ 2 + x ^ 3 + x ^ 4)" -> "(1 + x + 2 * x ^ 2 + x ^ 3 + x ^ 4)"];
Print["    In:        LeanVerify[<the Factor[x^10 - 1] certificate with one coefficient changed>]"];
Print["    Changed:   (1 + x + x ^ 2 + x ^ 3 + x ^ 4)  ->  (1 + x + 2 * x ^ 2 + x ^ 3 + x ^ 4)"];
t18 = AbsoluteTiming[LeanVerify[bad]];
v18 = t18[[2]];
Print["    Verified:  ", v18["Success"], If[v18["Success"] === False, "", "   <-- UNEXPECTED (expected False)"],
  "    Tier: 0    Lean errors: ",
  Length[Select[v18["Messages"], #["Severity"] === "error" &]]];
firstErr = Select[v18["Messages"], #["Severity"] === "error" &];
If[firstErr =!= {},
  Print["    Lean says: ", clip[StringReplace[First[StringSplit[firstErr[[1]]["Text"], "\n"]], "\t" -> " "]]]];
Print["    Time:      Lean ", secs[v18["Timing"]], ", whole call ", secs[t18[[1]]]];
record[18, "Tampered factorization certificate",
  If[v18["ExitCode"] === Missing["NotRun"], "not run", v18["Success"]], 0, t18[[1]],
  "n/a: a hand-edited source, checked as is by LeanVerify"];

(* ---------------------------------------------------------------------- *)
(* 19. Batch: rows 4-10 in one Lean run.                                   *)
(* ---------------------------------------------------------------------- *)

header[19, "Batch: rows 4-10 in ONE Lean run (session cache cleared first)"];
Lean`Private`$certCache = <||>;
singles = Total[Select[#[[5]] & /@ Select[$rows, 4 <= #[[1]] <= 10 &], NumericQ]];
t19 = AbsoluteTiming[batch = Certify[{
    Factor[x^10 - 1],
    PolynomialExtendedGCD[x^4 - 1, x^3 + 2 x^2 + 2 x + 1, x],
    Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1],
    Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x],
    Inverse[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}],
    Det[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, {3, 0, 1, 1}}],
    NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]}]];
Print["    In:        Certify[{Factor[...], PolynomialExtendedGCD[...], Implies[...], Solve[...], Inverse[...], Det[...], NullSpace[...]}]"];
Print["    Verified:  ", #["Verified"] & /@ batch];
Print["    Kinds:     ", #["Kind"] & /@ batch];
Print["    Time:      one Lean run ", secs[batch[[1]]["Timing"]], ", whole batch ", secs[t19[[1]]],
  "  vs  ", secs[singles], " for the 7 separate calls above (speedup ",
  ToString[NumberForm[singles/t19[[1]], {4, 1}]], "x)"];
again = AbsoluteTiming[Certify[Factor[x^10 - 1]]];
Print["    Repeat:    Certify[Factor[x^10 - 1]] again: Cached -> ", again[[2]]["Cached"],
  " in ", secs[again[[1]]]];
record[19, "Batch of rows 4-10, one Lean run",
  If[AllTrue[#["Verified"] & /@ batch, TrueQ], True, False], 0, t19[[1]], "as rows 4-10"];

(* ---------------------------------------------------------------------- *)
(* 20. Lean -> Mathilda: tactics that ask Mathilda and check in Lean.      *)
(* ---------------------------------------------------------------------- *)

$tacTrust = "nothing beyond each tactic's goal: Mathilda supplies only data, the kernel checks every step";
header[20, "Lean -> Mathilda: tactics in lean/ ask Mathilda, Lean checks the answer"];
Print["    import Mathilda"];
Print[""];
Print["    example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 1 := by"];
Print["      mathilda_linear_combination"];
Print["    example (x y : Rat) : 2*x*y ≤ x^2 + y^2 := by"];
Print["      mathilda_sos"];
Print["    example : ∃ c, Mathilda.PrattCertifies 274177 c := by"];
Print["      mathilda_prime"];
Print["    #mathilda \"Factor[x^10-1]\""];
Print[""];
Print["    The oracle answers (Lean`Internal`Oracle) for the three tactics above:"];
Print["      ", Lean`Internal`Oracle["LinearCombination", {2 v2^2, 1}, {{v1^2 + v2^2, 1}, {v1, v2}}, {v1, v2}]];
Print["      ", Lean`Internal`Oracle["SOS", {2 v1 v2, v1^2 + v2^2}, {}, {v1, v2}]];
Print["      ", Lean`Internal`Oracle["Pratt", 274177]];
Print["    Running: sh tests/scripts/lean_tactics.sh ./Mathilda"];
t20 = AbsoluteTiming[tac = Lean`Internal`RunProcess[{"sh",
     FileNameJoin[{$demoRoot, "tests", "scripts", "lean_tactics.sh"}],
     FileNameJoin[{$demoRoot, "Mathilda"}]}, None]];
tacLines = Select[StringSplit[tac["StandardOutput"], "\n"],
  StringMatchQ[#, ("PASS" | "FAIL" | "SKIP") ~~ ___] &];
(* lean_tactics.sh prints paths relative to the package; keep the report path-free. *)
Scan[Print["      ", #] &, tacLines];
Which[
  tac["ExitCode"] === 0,
    record[20, "Lean -> Mathilda tactics (lean_tactics.sh)", True, 0, t20[[1]], $tacTrust],
  tac["ExitCode"] === 77,
    record[20, "Lean -> Mathilda tactics (lean_tactics.sh)", "SKIPPED", "-", Missing[], {}],
  True,
    (Print["      exit code ", tac["ExitCode"], "   <-- UNEXPECTED"];
     record[20, "Lean -> Mathilda tactics (lean_tactics.sh)", False, 0, t20[[1]], $tacTrust])];
Print["    Time:      ", secs[t20[[1]]], " (lake build and the three test files)"];

(* ---------------------------------------------------------------------- *)
(* Summary.                                                                *)
(* ---------------------------------------------------------------------- *)

Print[""];
Print["Summary"];
Print["======="];
Print[StringPadLeft["#", 3], "  ", StringPadRight["Showcase", 74], StringPadRight["Verdict", 9],
  StringPadRight["Tier", 6], StringPadRight["Exact", 7], "seconds"];
Scan[Print[StringPadLeft[ToString[#[[1]]], 3], "  ",
    StringPadRight[StringTake[#[[2]], UpTo[72]], 74],
    StringPadRight[ToString[#[[3]]], 9], StringPadRight[ToString[#[[4]]], 6],
    StringPadRight[ToString[#[[7]]], 7], secs[#[[5]]]] &,
  $rows];
Print[""];
Print["What Lean did NOT check (TrustedSteps), per row:"];
$seen = <||>;
Scan[Function[r, Print[StringPadLeft[ToString[r[[1]]], 3], "  ",
    Which[
      StringQ[r[[6]]], r[[6]],
      r[[6]] === {}, "nothing: the theorem is the whole claim",
      KeyExistsQ[$seen, r[[6]]], "same as row " <> ToString[$seen[r[[6]]]],
      True, ($seen[r[[6]]] = r[[1]]; StringRiffle[r[[6]], "; "])]]],
  Select[$rows, #[[3]] =!= "SKIPPED" &]];

expected = Table[If[MemberQ[{17, 18}, r[[1]]], False, True], {r, $rows}];
ok = And @@ MapThread[#1[[3]] === "SKIPPED" || #1[[3]] === #2 &, {$rows, expected}];
Print[""];
Print["All rows as expected (17 and 18 False by design): ", ok];
Print["Total wall time: ", secs[AbsoluteTime[] - $t0]];
If[!ok, Exit[1]];
