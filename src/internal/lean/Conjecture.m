(* ====================================================================== *)
(* Conjecture.m -- experiment, guess, test, certify (loaded by Lean.m      *)
(* after Refute.m).                                                         *)
(*                                                                          *)
(* Conjecture[expr, n] evaluates the held expr exactly at n = s, s + 1, ... *)
(* ("FitPoints" values), guesses a closed form g(n) from them, tests g at   *)
(* "TestPoints" further values and then certifies the claim expr == g AS    *)
(* WRITTEN (the claim walker sees the user's expr, iterator names           *)
(* included) with Certify's pipeline.  The guess is only ever a candidate:  *)
(* "Status" is "Proved" only when a kernel accepted the certificate,        *)
(* "Checked" when Mathilda's exact re-check passed without a kernel          *)
(* verdict, "Refuted" when the certification found a counterexample,        *)
(* "Tested" otherwise (no certificate kind, or not verified).               *)
(*                                                                          *)
(* Guessers, in order (each reproduces every fit value exactly or gives     *)
(* $Failed): "Polynomial" (InterpolatingPolynomial, degree <= 10),          *)
(* "Geometric" (an exact Hankel NullSpace recurrence whose characteristic   *)
(* polynomial has rational roots; polynomial-times-r^n ansatz),             *)
(* "Fibonacci" (only when expr has Fibonacci: the smallest subset of        *)
(* {1, n, (-1)^n, F(n), F(n+1), F(n)^2, F(n) F(n+1), F(n+1)^2} that fits),   *)
(* "Rational" (p(n) - v(n) q(n) = 0, degrees <= 4, exact NullSpace) and     *)
(* "Hypergeometric" (named forms Binomial[2n, n], Binomial[2n, n]/(n+1),    *)
(* n!, (2n)! times c^n; else a rational term ratio with rational roots, as  *)
(* a product of Pochhammer symbols).                                        *)
(*                                                                          *)
(* The messages of this file are defined here (Lean.m's line budget); the   *)
(* file opens the package itself (Conjecture is exported by Lean.m).        *)
(* ====================================================================== *)

BeginPackage["Lean`"];
Begin["`Private`"];

Conjecture::var = "`1` is not a plain symbol without a value; Conjecture[expr, n] needs the variable n as a symbol (Clear it first)";
Conjecture::nonexact = "`1` at `2` == `3` is `4`, not an exact rational number; Conjecture guesses only exact sequences";
Conjecture::noguess = "no closed form of the forms `1` fits `2` at `3` == `4`, ..., `5` and its test points";
Conjecture::badopt = "`1` -> `2` is not valid: use `3`";
Conjecture::optx = "unknown option `1`; the options of Conjecture are `2`";
Conjecture::optvalue = "the option name `1` is a symbol with the value `2`, which would replace the name; Clear[`1`] or write the name as the string \"`3`\"";

Options[Conjecture] = {"Engine" -> {"Lean", "Exact"}, "FitPoints" -> 12, "TestPoints" -> 40,
  "Start" -> Automatic, "Forms" -> All, "Verify" -> True, TimeConstraint -> 60,
  "LeanExecutable" -> Automatic, "LeanProject" -> Automatic};

$conjectureForms = {"Polynomial", "Geometric", "Fibonacci", "Rational", "Hypergeometric"};

SetAttributes[Conjecture, HoldAll];

Conjecture[expr_, n_, opts:OptionsPattern[]] := Which[
  heldOptionProblemFor[Conjecture, HoldComplete[opts]] === $Failed, $Failed,
  !MatchQ[Unevaluated[n], _Symbol] || symbolProblem[n] =!= None,
    Message[Conjecture::var, HoldForm[n]]; $Failed,
  True, conjecture[HoldComplete[expr], n, {opts}]];

(* ---------------------------------------------------------------------- *)
(* Options and values                                                       *)
(* ---------------------------------------------------------------------- *)

conjOptions[opts_List] := Module[{keys = {"Engine", "FitPoints", "TestPoints", "Start",
    "Forms", "Verify", "TimeConstraint", "LeanExecutable", "LeanProject"}, o, bad},
  If[optionNameProblem[Conjecture, opts] =!= None, Return[$Failed, Module]];
  o = AssociationThread[keys, OptionValue[Conjecture, opts, #] & /@ keys];
  bad = Which[
    !(IntegerQ[o["FitPoints"]] && o["FitPoints"] >= 3),
      {"FitPoints", o["FitPoints"], "an integer >= 3"},
    !(IntegerQ[o["TestPoints"]] && o["TestPoints"] >= 0),
      {"TestPoints", o["TestPoints"], "an integer >= 0"},
    !(o["Start"] === Automatic || IntegerQ[o["Start"]]),
      {"Start", o["Start"], "Automatic or an integer"},
    !(o["Forms"] === All || ListQ[o["Forms"]] && o["Forms"] =!= {} &&
        SubsetQ[$conjectureForms, o["Forms"]]),
      {"Forms", o["Forms"], "All or a list of " <> StringRiffle[$conjectureForms, ", "]},
    engineList[o["Engine"]] === $Failed,
      {"Engine", o["Engine"], "\"Lean\", \"Exact\", \"Mathilda\" or a list of them"},
    True, None];
  If[bad =!= None, Message[Conjecture::badopt, Sequence @@ bad]; Return[$Failed, Module]];
  o];

(* The lower limit of an outer Sum[t, {k, a, n}] (1 for {k, n}), else 0. *)
conjStart[hc_HoldComplete, n_] := Replace[hc, {HoldComplete[Sum[_, {_, a_Integer, n}]] :> a,
  HoldComplete[Sum[_, {_, n}]] :> 1, _ :> 0}];

(* expr at n = m, evaluated (quietly) from the held input. *)
conjValue[hc_HoldComplete, n_, m_Integer] := Quiet[Replace[hc /. n -> m, HoldComplete[e_] :> e]];

(* ---------------------------------------------------------------------- *)
(* Guessers: vals are the values at n = s, s + 1, ...; each gives a closed  *)
(* form in n reproducing every value exactly, or $Failed.                   *)
(* ---------------------------------------------------------------------- *)

fitsQ[g_, n_, s_, vals_] := And @@ MapIndexed[Quiet[g /. n -> s + First[#2] - 1] === #1 &, vals];
mpow[m_, 0] := 1;  (* m^0 without 0^0 *)
mpow[m_, j_] := m^j;

guessPolynomial[vals_, s_, n_] := Module[{pts = Transpose[{s - 1 + Range[Length[vals]], vals}], p},
  Do[p = Expand[InterpolatingPolynomial[Take[pts, d + 1], n]];
    If[fitsQ[p, n, s, vals], Return[Factor[p], Module]],
    {d, 0, Min[10, Length[vals] - 2]}];
  $Failed];

(* The coefficients {c0, ..., cr} of the smallest-order linear recurrence
   c0 v(i) + ... + cr v(i + r) == 0 holding on every value (an exact
   NullSpace of the Hankel matrix, with at least r + 2 equations), or
   $Failed. *)
findRecurrence[vals_, maxOrder_] := Module[{l = Length[vals], ns},
  Do[ns = NullSpace[Table[vals[[i + j]], {i, 1, l - r}, {j, 0, r}]];
    If[ns =!= {}, Return[First[ns], Module]],
    {r, 1, Min[maxOrder, Floor[(l - 2)/2]]}];
  $Failed];

(* The roots {{r, mult}, ...} of a polynomial in x when all are nonzero
   rationals, else $Failed. *)
rationalRoots[poly_, x_] := Module[{rs = {}},
  Do[Which[FreeQ[f[[1]], x], Null,
      Exponent[f[[1]], x] === 1 && Coefficient[f[[1]], x, 0] =!= 0,
        AppendTo[rs, {-Coefficient[f[[1]], x, 0]/Coefficient[f[[1]], x, 1], f[[2]]}],
      True, Return[$Failed, Module]],
    {f, FactorList[poly]}];
  rs];

(* An exact linear fit vals == Sum c_i basis_i(m) (basis: functions of m)
   by the normal equations; the coefficients, or $Failed when the basis is
   degenerate on the points or does not reproduce every value. *)
linearFit[basis_List, vals_, s_] := Module[{a, ata, c},
  a = Table[#[s + i - 1] & /@ basis, {i, Length[vals]}];
  ata = Transpose[a] . a;
  If[Det[ata] === 0, Return[$Failed, Module]];
  c = LinearSolve[ata, Transpose[a] . vals];
  If[a . c === vals, c, $Failed]];

guessGeometric[vals_, s_, n_] := Module[{c, x, rs, terms, cs},
  c = findRecurrence[vals, 6];
  If[c === $Failed, Return[$Failed, Module]];
  rs = rationalRoots[c . (mpow[x, #] & /@ Range[0, Length[c] - 1]), x];
  If[rs === $Failed || rs === {}, Return[$Failed, Module]];
  terms = Flatten[Table[{r[[1]], j}, {r, rs}, {j, 0, r[[2]] - 1}], 1];
  cs = linearFit[With[{r = #[[1]], j = #[[2]]}, Function[m, mpow[m, j] r^m]] & /@ terms, vals, s];
  If[cs === $Failed, Return[$Failed, Module]];
  Total[Table[With[{r = rt[[1]]},
      Factor[Total[MapThread[If[#1[[1]] === r, #2 mpow[n, #1[[2]]], 0] &, {terms, cs}]]] *
        If[r === 1, 1, r^n]],
    {rt, rs}]]];

(* The Fibonacci basis, as functions of m and as expressions in n. *)
fibBasis[n_] := {1, n, (-1)^n, Fibonacci[n], Fibonacci[n + 1], Fibonacci[n]^2,
  Fibonacci[n] Fibonacci[n + 1], Fibonacci[n + 1]^2};
$fibFunctions = {Function[m, 1], Function[m, m], Function[m, (-1)^m],
  Function[m, Fibonacci[m]], Function[m, Fibonacci[m + 1]], Function[m, Fibonacci[m]^2],
  Function[m, Fibonacci[m] Fibonacci[m + 1]], Function[m, Fibonacci[m + 1]^2]};
guessFibonacci[vals_, s_, n_] := Module[{c},
  Do[If[Length[sub] <= Length[vals] - 2,
      c = linearFit[$fibFunctions[[sub]], vals, s];
      If[c =!= $Failed, Return[c . fibBasis[n][[sub]], Module]]],
    {sub, Rest[Subsets[Range[8]]]}];
  $Failed];

(* p(n) / q(n), deg p <= 4, 1 <= deg q <= 4, smallest total degree first. *)
guessRational[vals_, s_, n_] := Module[{l = Length[vals], ms, ns, p, q, g},
  ms = s - 1 + Range[l];
  Do[Do[With[{dp = d - dq},
      If[0 <= dp <= 4 && dp + dq + 2 <= l - 1,
        ns = NullSpace[MapThread[Function[{m, v}, Join[mpow[m, #] & /@ Range[0, dp],
            -v (mpow[m, #] & /@ Range[0, dq])]], {ms, vals}]];
        Do[p = Take[v, dp + 1] . (mpow[n, #] & /@ Range[0, dp]);
          q = Drop[v, dp + 1] . (mpow[n, #] & /@ Range[0, dq]);
          If[q =!= 0 && AllTrue[ms, (q /. n -> #) =!= 0 &],
            g = Factor[Cancel[p/q]];
            If[fitsQ[g, n, s, vals], Return[g, Module]]],
          {v, ns}]]],
      {dq, 1, Min[4, d]}],
    {d, 1, 8}];
  $Failed];

$namedHyper = {Function[m, Binomial[2 m, m]], Function[m, Binomial[2 m, m]/(m + 1)],
  Function[m, m!], Function[m, (2 m)!]};

guessHypergeometric[vals_, s_, n_] := Module[{l = Length[vals], w, c, rho, nu, de, fac},
  If[MemberQ[vals, 0] || s < 0, Return[$Failed, Module]];
  (* Named forms times c^n: vals / f geometric with a constant ratio. *)
  Do[w = vals / Table[f[s + i - 1], {i, l}];
    c = w[[2]]/w[[1]];
    If[Rest[w] === c Most[w], Return[w[[1]] c^(n - s) f[n], Module]],
    {f, $namedHyper}];
  (* The term ratio as a rational function with rational roots. *)
  rho = guessPolynomial[Rest[vals]/Most[vals], s, n];
  If[rho === $Failed, rho = guessRational[Rest[vals]/Most[vals], s, n]];
  If[rho === $Failed, Return[$Failed, Module]];
  rho = Together[rho];
  nu = rationalRoots[Numerator[rho], n]; de = rationalRoots[Denominator[rho], n];
  If[nu === $Failed || de === $Failed, Return[$Failed, Module]];
  (* rho = C Prod (n - alpha)^mu / Prod (n - beta)^nu, so v(n) = v(s) C^(n - s)
     Prod Pochhammer[s - alpha, n - s]^mu / Prod Pochhammer[s - beta, n - s]^nu. *)
  fac = Together[rho / ((Times @@ ((n - #[[1]])^#[[2]] & /@ nu)) /
    (Times @@ ((n - #[[1]])^#[[2]] & /@ de)))];
  If[!FreeQ[fac, n], Return[$Failed, Module]];
  vals[[1]] fac^(n - s) (Times @@ (Pochhammer[s - #[[1]], n - s]^#[[2]] & /@ nu)) /
    (Times @@ (Pochhammer[s - #[[1]], n - s]^#[[2]] & /@ de))];

guessForm["Polynomial", vals_, s_, n_, _] := guessPolynomial[vals, s, n];
guessForm["Geometric", vals_, s_, n_, _] := guessGeometric[vals, s, n];
guessForm["Fibonacci", vals_, s_, n_, fibQ_] := If[fibQ, guessFibonacci[vals, s, n], $Failed];
guessForm["Rational", vals_, s_, n_, _] := guessRational[vals, s, n];
guessForm["Hypergeometric", vals_, s_, n_, _] := guessHypergeometric[vals, s, n];

(* ---------------------------------------------------------------------- *)
(* Certification of expr == guess                                           *)
(* ---------------------------------------------------------------------- *)

(* The messages whose name is recorded as the reason when Certify gives no
   ProofObject, highest priority first.  $MessageList does not see quieted
   messages, so each is caught by its own Check inside the Quiet (the
   outermost Check runs last, so the first name of the list that fired
   wins). *)
$conjReasons = {HoldForm[Certify::eval], HoldForm[Certify::value], HoldForm[Certify::context],
  HoldForm[Certify::fracpow], HoldForm[Certify::tier], HoldForm[Certify::nomathlib],
  HoldForm[Certify::unsafe], HoldForm[Certify::nocert]};
$conjWhy = None;
SetAttributes[checkChain, HoldFirst];
checkChain[e_, {}] := e;
checkChain[e_, {h_HoldForm, rest___}] := Replace[h, HoldForm[mn_] :>
  Check[checkChain[e, {rest}], $conjWhy = ToString[h], mn]];

(* The claim expr == g as written: Join puts the held input and the guess
   side by side (no literal Sum and no With around the held input, which
   would rename its iterator). *)
conjClaim[hc_HoldComplete, g_] := Replace[Join[hc, HoldComplete[g]],
  HoldComplete[a_, b_] :> HoldComplete[Equal[a, b]]];

(* {ProofObject or $Failed, reason}: Certify's pipeline on the claim, every
   message quieted. *)
conjCertify[claim_HoldComplete, o_] := Module[{r = $Failed, copts},
  copts = {"Engine" -> o["Engine"], "Verify" -> o["Verify"],
    TimeConstraint -> o["TimeConstraint"], "LeanExecutable" -> o["LeanExecutable"],
    "LeanProject" -> o["LeanProject"]};
  Block[{$conjWhy = None},
    Quiet[checkChain[r = certifyItems[{claim}, copts], $conjReasons]];
    {If[ListQ[r], First[r], r], If[$conjWhy === None, "Certify::nocert", $conjWhy]}]];

(* The status of a certification (plan section 3.7, item 6). *)
conjStatus[p_ProofObject] := With[{a = First[p], v = p["Verified"]}, Which[
  KeyExistsQ[a, "Refuted"], "Refuted",
  v === True, "Proved",
  v =!= False && (v === Missing["Exact", True] ||
    Lookup[Lookup[a, "Engines", <||>], "Exact", None] === Missing["Exact", True]), "Checked",
  True, "Tested"]];
conjStatus[_] := "Tested";

(* ---------------------------------------------------------------------- *)
(* Conjecture                                                               *)
(* ---------------------------------------------------------------------- *)

conjecture[hc_HoldComplete, n_, opts_List] := Module[{t0 = AbsoluteTime[], o, s, fit, tot,
    vals, bad, forms, fibQ, g = $Failed, form, proof, why, claim, a, input, pq},
  o = conjOptions[opts];
  If[o === $Failed, Return[$Failed, Module]];
  input = Replace[hc, HoldComplete[e_] :> HoldForm[e]];
  s = If[o["Start"] === Automatic, conjStart[hc, n], o["Start"]];
  fit = o["FitPoints"]; tot = fit + o["TestPoints"];
  (* The experiment: every value exact. *)
  vals = Table[conjValue[hc, n, s + i - 1], {i, tot}];
  bad = SelectFirst[Range[tot], !MatchQ[vals[[#]], _Integer | _Rational] &, None];
  If[bad =!= None,
    Message[Conjecture::nonexact, input, n, s + bad - 1, vals[[bad]]];
    Return[$Failed, Module]];
  forms = If[o["Forms"] === All, $conjectureForms,
    Select[$conjectureForms, MemberQ[o["Forms"], #] &]];
  fibQ = !FreeQ[hc, Fibonacci];
  (* Guess on the fit points, test on all of them; the first form that passes. *)
  Do[With[{c = guessForm[f, Take[vals, fit], s, n, fibQ]},
      If[c =!= $Failed && fitsQ[c, n, s, vals], g = c; form = f; Break[]]],
    {f, forms}];
  If[g === $Failed,
    Message[Conjecture::noguess, forms, input, n, s, s + fit - 1]; Return[$Failed, Module]];
  claim = conjClaim[hc, g];
  {proof, why} = conjCertify[claim, o];
  pq = Head[proof] === ProofObject;
  a = <|"Kind" -> "Conjecture", "Input" -> input, "Variable" -> n, "Guess" -> g, "Form" -> form,
    "FittedOn" -> {s, s + fit - 1}, "TestedOn" -> {s, s + tot - 1},
    "Status" -> conjStatus[proof],
    "Verified" -> If[pq, proof["Verified"], Missing["NoCertificate", why]],
    "Proof" -> If[pq, proof, Missing["NoCertificate"]],
    "Theorem" -> If[pq, proof["Theorem"],
      Replace[Join[HoldComplete @@ {n, s}, claim], HoldComplete[v_, s0_, Equal[l_, r_]] :>
        HoldForm[ForAll[v, v >= s0, l == r]]]],
    "Engine" -> o["Engine"], "Timing" -> AbsoluteTime[] - t0|>;
  issueObject[a]];

(* The summary box (Lean.m's Format):
   ProofObject[«Conjecture» Proved: expr == guess | Lean Tier 0 + Exact | 2.4 s]. *)
conjectureText[a_] := Module[{p = a["Proof"]}, StringJoin[
  "ProofObject[\[LeftGuillemet]Conjecture\[RightGuillemet] ",
  If[authenticQ[a], a["Status"] <> ": " <> ToString[a["Input"]] <> " == " <> ToString[a["Guess"]],
    "Verified: " <> ToString[Missing["Unauthenticated"]]], " | ",
  If[Head[p] === ProofObject,
    StringRiffle[engineLabel[#, First[p]] & /@ Flatten[{First[p]["Engine"]}], " + "],
    "no certificate"],
  " | ", timingText[a["Timing"]], "]"]];

End[];

EndPackage[];
