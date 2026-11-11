(* ====================================================================== *)
(* ExactBackend.m -- the "Exact" backend of Certify (loaded by Lean.m).    *)
(*                                                                          *)
(* An always-available, in-tree re-check of a certificate by Mathilda's     *)
(* own exact rational arithmetic and polynomial Expand.  It is NOT a kernel *)
(* proof: the checker is the system that produced the certificate, so its  *)
(* verdict is Missing["Exact", True | False] (never True; the pipeline      *)
(* enforces this for every non-kernel backend), or Missing["NotChecked"]    *)
(* when no decisive exact check exists for the shape (an identity with an   *)
(* inequality hypothesis, a check that timed out).  False is given only on  *)
(* a concrete failing exact computation.                                    *)
(*                                                                          *)
(* Every check is recomputed from the obligation's held "Statement" plus    *)
(* the certificate data, never from the generator's own "PreCheck".  The    *)
(* encoder (exactEncode) only builds the checks: each is                    *)
(* <|"Check" -> readable text, "Test" -> Hold[predicate]|>, evaluated by    *)
(* exactRun under TimeConstrained, so Verify -> False runs nothing.  All    *)
(* messages are in Lean.m.  The file opens the package itself (every symbol *)
(* is in Lean`Private`).                                                    *)
(* ====================================================================== *)

BeginPackage["Lean`"];
Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Helpers                                                                 *)
(* ---------------------------------------------------------------------- *)

(* Readable text of a held expression, clipped. *)
exText[e_] := With[{s = ToString[e]},
  If[StringLength[s] > 160, StringTake[s, 157] <> "...", s]];

exEqText[Hold[l_, r_]] := exText[HoldForm[l == r]];
exCheck[text_String, test_Hold] := <|"Check" -> text, "Test" -> test|>;
exUndecided[why_String] := exCheck[why, Hold[Missing["NotChecked"]]];

(* {Hold[hypothesis], Hold[conclusion]} of a held ForAll[vs, (h,) s],
   Implies[h, s] or s (hypothesis True when there is none).  The statement
   is never evaluated as a whole: Equal of identical sides would be True. *)
exParts[HoldForm[ForAll[_, h_, s_]]] := {Hold[h], Hold[s]};
exParts[HoldForm[ForAll[_, s_]]] := {Hold[True], Hold[s]};
exParts[HoldForm[Implies[h_, s_]]] := {Hold[h], Hold[s]};
exParts[HoldForm[s_]] := {Hold[True], Hold[s]};

(* The held equations l == r of a held conclusion (a conjunction or one). *)
exEquations[c_Hold] := Cases[c, HoldPattern[Equal[l_, r_]] :> Hold[l, r], {1, Infinity}];

(* l - r of a held pair Hold[l, r], evaluated. *)
exDiff[Hold[l_, r_]] := l - r;

(* The curried hypotheses and goal of ForAll[vs, h1 -> h2 -> ... -> l == r]:
   {{Hold[l1, r1], ...}, Hold[l, r]}, or $Failed. *)
SetAttributes[exImp, HoldAllComplete];
exImp[ForAll[_, b_]] := exImp[b];
exImp[Implies[Equal[a_, b_], rest_]] := With[{t = exImp[rest]},
  If[t === $Failed, $Failed, {Prepend[t[[1]], Hold[a, b]], t[[2]]}]];
exImp[Equal[a_, b_]] := {{}, Hold[a, b]};
exImp[_] := $Failed;

(* A nonzero rational function of x with rational coefficients: a concrete
   nonzero value, as opposed to a normal form Together could not reach. *)
exRationalNonzeroQ[d_, x_] := d =!= 0 && FreeQ[d, _Real | _Complex] &&
  FreeQ[d, Power[_?NumberQ, _Rational]] && SubsetQ[{x}, Variables[d]] &&
  polyQ[Numerator[d], {x}] && polyQ[Denominator[d], {x}];

(* ---------------------------------------------------------------------- *)
(* Check predicates (evaluated by exactRun)                                *)
(* ---------------------------------------------------------------------- *)

(* q lies in the ideal of the ps: cofactors found by Mathilda's search, then
   the identity checked; no cofactors is "not decided", never False. *)
exIdealQ[q_, ps_List, vars_List] := With[{c = cofactors[Expand[q], Expand[ps], vars]},
  If[c === $Failed, Missing["NotChecked"], Expand[q - c . ps] === 0]];

exIdealCofactorsQ[hs_List, goal_Hold, cs_List] := Length[hs] === Length[cs] &&
  Expand[exDiff[goal] - Sum[cs[[i]] exDiff[hs[[i]]], {i, Length[hs]}]] === 0;

(* lc * Prod (x - r_i)^m_i (inert) equals l - r, lc a nonzero number, and its
   roots are exactly the alternatives of the statement. *)
exRootsQ[eq_Hold, roots_List, fact_, x_] := Module[{fs, lin},
  fs = If[Head[fact] === iTimes, List @@ fact, {fact}];
  lin = Replace[fs, {c_?ratNumberQ /; c =!= 0 :> "const",
      iPlus[x, c_?ratNumberQ] :> -c,
      iPower[iPlus[x, c_?ratNumberQ], m_Integer?Positive] :> -c,
      _ -> $Failed}, {1}];
  lin = DeleteCases[lin, "const"];
  !MemberQ[lin, $Failed] && Expand[exDiff[eq] - inertValue[fact]] === 0 &&
    Sort[DeleteDuplicates[lin]] === Sort[DeleteDuplicates[roots]]];

(* The closed statement with Exact-side meanings of the prelude placeholders
   (det by a cofactor expansion independent of the Det builtin). *)
exHasShape[m_, r_, c_] := MatrixQ[m] && Dimensions[m] === {r, c};
exHasLength[v_, c_] := VectorQ[v] && Length[v] === c;
exactLaplaceDet[m_?MatrixQ] /; Length[m] === Length[First[m]] && Length[m] <= 6 :=
  If[Length[m] === 1, m[[1, 1]],
    Sum[If[m[[1, j]] === 0, 0,
      (-1)^(j + 1) m[[1, j]] exactLaplaceDet[Drop[#, {j}] & /@ Rest[m]]], {j, Length[m]}]];
(* (Solve's soundness statement is a HoldForm inside the HoldForm.) *)
exClosed[held_HoldForm] := ReleaseHold[ReleaseHold[held /. {hasShape -> exHasShape,
  matMul -> Dot, idMat -> IdentityMatrix, matVec -> Dot, hasLength -> exHasLength,
  det -> exactLaplaceDet}]];

(* The target of a held relation a >= b (a - b) or a <= b (b - a). *)
exTarget[Hold[GreaterEqual[a_, b_]]] := a - b;
exTarget[Hold[LessEqual[a_, b_]]] := b - a;
exTarget[_] := $Failed;

(* target == Sum d l^2 mult, every d > 0 and every mult 1 or a hypothesis
   polynomial of the statement. *)
exSquaresQ[target_, terms_List, gs_List] := target =!= $Failed &&
  AllTrue[terms, ratNumberQ[#[[1]]] && #[[1]] > 0 &] &&
  AllTrue[terms, #[[3]] === 1 || MemberQ[gs, Expand[#[[3]]]] &] &&
  Expand[target - Total[(#[[1]] #[[2]]^2 #[[3]]) & /@ terms]] === 0;

exSOSQ[conc_Hold, hyp_Hold, terms_List] := Module[{gs},
  gs = Expand[ineqPoly /@ condList[ReleaseHold[hyp]]];
  !MemberQ[gs, $Failed] && exSquaresQ[exTarget[conc], terms, gs]];

(* KKT: constraints of the statement as {kind, polynomial}. *)
exCons[c_] := Map[Function[r, Which[
    Length[r] =!= 2, $Failed,
    Head[r] === GreaterEqual, {"ineq", r[[1]] - r[[2]]},
    Head[r] === LessEqual, {"ineq", r[[2]] - r[[1]]},
    Head[r] === Equal, {"eq", r[[1]] - r[[2]]},
    True, $Failed]], condList[c]];

(* Each eliminated variable v -> sl is forced by an equality of the
   statement (h = c (v - sl) after the earlier eliminations, c a nonzero
   number), so every feasible point is a point of the reduced problem; then
   F - M = Sum lambda_j g_j + Sum d l^2 identically, lambda >= 0. *)
exKKTQ[fa_Hold, ob_] := Module[{vars, cons, rel, sgn, f, m, eqs, ineqs, subs = {}, ok, gr,
    lam = ob["Multipliers"], terms = ob["Squares"]},
  {cons, rel} = Replace[fa, {Hold[ForAll[_, c_, r_]] :> {c, Hold[r]},
      Hold[ForAll[_, r_]] :> {True, Hold[r]}, _ -> {$Failed, $Failed}}];
  If[rel === $Failed, Return[False, Module]];
  {sgn, f, m} = Replace[rel, {Hold[GreaterEqual[a_, b_]] :> {1, a, b},
      Hold[LessEqual[a_, b_]] :> {-1, a, b}, _ -> {$Failed, 0, 0}}];
  cons = exCons[cons];
  If[sgn === $Failed || MemberQ[cons, $Failed], Return[False, Module]];
  eqs = Select[cons, #[[1]] === "eq" &][[All, 2]];
  ineqs = Select[cons, #[[1]] === "ineq" &][[All, 2]];
  Do[
    ok = AnyTrue[eqs, Function[h, With[{hh = Expand[h /. subs]},
      With[{c = Coefficient[hh, e[[1]], 1]},
        ratNumberQ[c] && c =!= 0 && Expand[hh - c (e[[1]] - e[[2]])] === 0]]]];
    If[!ok, Return[False, Module]];
    subs = Append[subs /. e[[1]] -> e[[2]], e[[1]] -> e[[2]]],
    {e, ob["Eliminated"]}];
  gr = Expand[ineqs /. subs];
  Length[lam] === Length[gr] && AllTrue[lam, ratNumberQ[#] && # >= 0 &] &&
    AllTrue[terms, #[[3]] === 1 &] &&
    exSquaresQ[Expand[sgn f /. subs] - sgn m - lam . gr, terms, {}]];

(* The closed conjuncts of a held And[..., facts...] all True. *)
exFacts[facts_Hold] := With[{vs = List @@ facts},
  Which[AllTrue[vs, TrueQ], True, MemberQ[vs, False], False, True, Missing["NotChecked"]]];

(* A flat Pratt certificate for p, replayed: each step {q, a, q - 1
   factored} has the factorization right, every odd prime factor certified
   by an earlier step (2 is the base), and a Lucas witness a. *)
exLucasQ[q_, a_, fs_] := PowerMod[a, q - 1, q] === 1 &&
  AllTrue[fs[[All, 1]], PowerMod[a, (q - 1)/#, q] =!= 1 &];
exPrattQ[p_Integer, steps_List] := steps =!= {} && Last[steps][[1]] === p &&
  AllTrue[Range[Length[steps]], Function[i, With[{q = steps[[i, 1]], a = steps[[i, 2]],
      fs = steps[[i, 3]]},
    IntegerQ[q] && q > 2 && IntegerQ[a] && Times @@ (#[[1]]^#[[2]] & /@ fs) === q - 1 &&
      AllTrue[fs[[All, 1]], # === 2 || MemberQ[steps[[1 ;; i - 1, 1]], #] &] &&
      exLucasQ[q, a, fs]]]];

(* Sums: an exponential polynomial in n (powers c^(alpha n + beta), c
   rational, over prime bases z_p = p^n, (-1)^n reduced modulo z^2 - 1) is
   zero as a formal expression. *)
exPow[c_, ex_, n_] := Module[{al = Coefficient[ex, n, 1], be, fac},
  be = Expand[ex - al n];
  If[!IntegerQ[al] || !IntegerQ[be] || !polyQ[ex, {n}] || Exponent[ex, n] > 1 || c === 0,
    Return[$Failed, Module]];
  fac[k_Integer, s_] := If[k === 1, 1, Times @@ (ez[#[[1]]]^(s #[[2]] al) & /@ FactorInteger[k])];
  c^be fac[Numerator[Abs[c]], 1] fac[Denominator[Abs[c]], -1] If[c < 0, ez[-1]^al, 1]];
expPolyZeroQ[e_, n_] := Module[{r},
  r = Expand[e] /. Power[c_?ratNumberQ, ex_] /; !FreeQ[ex, n] :> exPow[c, ex, n];
  If[!FreeQ[r, $Failed] || !FreeQ[r, Power[_, ex_] /; !FreeQ[ex, n]], Return[False, Module]];
  r = Expand[r] /. Power[ez[-1], k_Integer] :> ez[-1]^Mod[k, 2];
  Together[Expand[r]] === 0];

(* Fibonacci[n + j] = F(j - 1) f0 + F(j) f1 for every integer j (F(-1) = 1),
   f0 = F(n) and f1 = F(n + 1) fresh: a Fibonacci step becomes a polynomial
   identity in f0, f1.  Other Fibonacci arguments are left alone. *)
exFibLin[e_, n_] := e /. Fibonacci[u_] /; IntegerQ[Expand[u - n]] :>
  With[{j = Expand[u - n]}, Fibonacci[j - 1] exFib0 + Fibonacci[j] exFib1];

(* Base and step of the induction; when the step is not an exponential-
   polynomial identity, the claim at n = a .. a + 10 decides: a concrete
   mismatch is False, otherwise not decided. *)
exSumBaseQ[t_, k_, a_, n_, c_] := Together[(c /. n -> a) - (t /. k -> a)] === 0;
exSumStepQ[t_, k_, a_, n_, c_] := Module[{bad},
  If[expPolyZeroQ[exFibLin[(c /. n -> n + 1) - c - (t /. k -> n + 1), n], n],
    Return[True, Module]];
  bad = SelectFirst[Range[a, a + 10],
    Together[Total[Table[t /. k -> j, {j, a, #}]] - (c /. n -> #)] =!= 0 &, None];
  If[bad === None, Missing["NotChecked"], False]];

(* D[F, x] == f: zero by Together; a nonzero rational function is False. *)
exDerivQ[F_, f_, x_] := With[{d = Together[D[F, x] - f]},
  Which[d === 0, True, exRationalNonzeroQ[d, x], False, True, Missing["NotChecked"]]];

(* ---------------------------------------------------------------------- *)
(* Encoder: the neutral IR -> checks, one rule per obligation "Type"       *)
(* ---------------------------------------------------------------------- *)

exactEncode[ir_, o_] := <|"Checks" -> exactChecks[ir["Obligation"]["Type"], ir]|>;

(* An identity: no assumptions, the polynomial identity; equality
   assumptions, ideal membership; any other assumption has no decision
   procedure here.  Cleared denominators must not vanish identically. *)
exactChecks["Identity", ir_] := Module[{ob = ir["Obligation"], conc, eq, assum, main, dens},
  conc = exParts[ob["Statement"]][[2]];
  eq = With[{es = exEquations[conc]}, If[es === {}, $Failed, First[es]]];
  assum = ob["Assumptions"];
  If[eq === $Failed, Return[{exUndecided["the statement is not an equation"]}, Module]];
  main = Which[
    assum === {}, With[{e = eq}, exCheck["Expand[lhs - rhs] === 0 for " <>
      exEqText[e], Hold[Expand[exDiff[e]] === 0]]],
    AllTrue[assum, Head[#] === Equal &],
      With[{e = eq, ps = (#[[1]] - #[[2]]) & /@ assum, vs = ob["Variables"]},
        exCheck["lhs - rhs lies in the ideal of the hypotheses " <> exText[assum] <>
          " (cofactors c with Expand[lhs - rhs - c . p] === 0)", Hold[exIdealQ[exDiff[e], ps, vs]]]],
    True, exUndecided["an identity under the assumptions " <> exText[assum] <>
      " (not only equations): no exact decision procedure"]];
  dens = Cases[ir["Hypotheses"], HoldPattern[HoldForm[Unequal[d_, 0]]] :> Hold[d]];
  Prepend[With[{d = #}, exCheck["the cleared denominator " <> exText[HoldForm @@ d] <>
      " is not identically zero", Hold[Expand[ReleaseHold[d]] =!= 0]]] & /@ dens, main]];

(* Polynomial identities (Bezout, PolynomialReduce, the Tier 0 cleared
   antiderivative identity): each equation of the statement expands to 0. *)
exactChecks["PolyIdentities", ir_] := With[{eqs = exEquations[exParts[ir["Obligation"]["Statement"]][[2]]]},
  If[eqs === {}, {exUndecided["no equation in the statement"]},
    With[{e = #}, exCheck["Expand[lhs - rhs] === 0 for " <> exEqText[e],
      Hold[Expand[exDiff[e]] === 0]]] & /@ eqs]];

exactChecks["IdealImplication", ir_] := With[{p = exImp @@ ir["Obligation"]["Statement"],
    cs = ir["Obligation"]["Cofactors"]},
  If[p === $Failed, {exUndecided["the statement is not a chain of equations"]},
    {With[{hs = p[[1]], g = p[[2]]}, exCheck["Expand[(l - r) - Sum c_i (l_i - r_i)] === 0 " <>
      "for the cofactors c = " <> exText[cs], Hold[exIdealCofactorsQ[hs, g, cs]]]]}]];

exactChecks["RootsComplete", ir_] := Module[{ob = ir["Obligation"], parts},
  parts = Replace[exParts[ob["Statement"]][[2]], {Hold[Equivalent[Equal[l_, r_], alts_]] :>
    {Hold[l, r], Cases[Hold[alts], HoldPattern[Equal[_, v_]] :> v, {1, Infinity}]},
    _ -> $Failed}];
  If[parts === $Failed, Return[{exUndecided["the statement is not an equivalence"]}, Module]];
  {With[{e = parts[[1]], rs = parts[[2]], fact = ob["Factorization"], x = ob["Variable"]},
    exCheck["l - r == " <> exText[inertHeld[fact]] <> " (a nonzero constant times linear " <>
      "factors) and its roots are exactly " <> exText[rs], Hold[exRootsQ[e, rs, fact, x]]]]}];

exactChecks["ClosedFact", ir_] := With[{s = ir["Obligation"]["Statement"]},
  {exCheck["the closed statement evaluates to True with exact rationals: " <> exText[s],
    Hold[exClosed[s]]]}];

exactChecks["WeightedSOS", ir_] := Module[{ob = ir["Obligation"], p},
  p = exParts[ob["Statement"]];
  {With[{c = p[[2]], h = p[[1]], terms = ob["Squares"]},
    exCheck["target == Sum d l^2 mult with every d > 0 and every mult 1 or a hypothesis " <>
      "polynomial (" <> ToString[Length[terms]] <> " squares)", Hold[exSOSQ[c, h, terms]]]]}];

exactChecks["KKT", ir_] := Module[{ob = ir["Obligation"], parts},
  parts = Replace[ob["Statement"], {HoldForm[And[fa_, rest___]] :> {Hold[fa], Hold[rest]},
    _ -> $Failed}];
  If[parts === $Failed, Return[{exUndecided["the statement is not a conjunction"]}, Module]];
  {With[{f = parts[[2]]}, exCheck["the optimum is feasible and attains the value: " <>
      exText[Replace[f, Hold[x___] :> HoldForm[And[x]]]], Hold[exFacts[f]]]],
   With[{fa = parts[[1]], o2 = ob}, exCheck["F - M - lambda . g == Sum d l^2 identically " <>
      "(lambda = " <> exText[ob["Multipliers"]] <> " >= 0, eliminated equalities forced by " <>
      "the constraints)", Hold[exKKTQ[fa, o2]]]]}];

exactChecks["Primality", ir_] := With[{ob = ir["Obligation"]}, With[{nv = ob["N"]},
  Switch[ob["Form"],
    "Factorization", With[{fi = ob["Factorization"], certs = ob["Certificates"],
        ps = Cases[ob["Statement"], HoldPattern[PrimeQ[p_]] :> p, {1, Infinity}]},
      Join[{exCheck[ToString[nv] <> " == " <> exText[fi] <> " and the statement lists " <>
          "these primes", Hold[Times @@ (#[[1]]^#[[2]] & /@ fi) === nv &&
          Sort[ps] === Sort[fi[[All, 1]]]]]},
        MapThread[With[{p = #1, c = #2}, exCheck["Pratt certificate of " <> ToString[p] <>
          " replayed (" <> ToString[Length[c]] <> " Lucas steps)", Hold[exPrattQ[p, c]]]] &,
          {DeleteCases[fi[[All, 1]], 2], certs}]]],
    "Prime", With[{c = First[ob["Certificates"]]}, {exCheck["Pratt certificate of " <>
      ToString[nv] <> " replayed (" <> ToString[Length[c]] <> " Lucas steps)",
      Hold[exPrattQ[nv, c]]]}],
    "Two", {exCheck["the divisors of 2 are 1 and 2", Hold[Select[Range[2],
      Mod[2, #] === 0 &] === {1, 2}]]},
    "Small", {exCheck[ToString[nv] <> " < 2", Hold[nv < 2]]},
    "Composite", With[{a = ob["Divisor"]}, {exCheck["1 < " <> ToString[a] <> " < " <>
      ToString[nv] <> " and " <> ToString[a] <> " divides " <> ToString[nv],
      Hold[1 < a < nv && Mod[nv, a] === 0]]}],
    _, {exUndecided["unknown primality form"]}]]];

exactChecks["Induction", ir_] := With[{ob = ir["Obligation"]},
  With[{t = ob["Term"], k = ob["Variable"], a = ob["LowerLimit"], n = ob["SumVariable"],
      c = ob["ClosedForm"]},
    {exCheck["base: closed form at n = " <> ToString[a] <> " equals the term at k = " <>
       ToString[a], Hold[exSumBaseQ[t, k, a, n, c]]],
     exCheck["step: C(n + 1) - C(n) == t(n + 1) as an exponential-polynomial identity " <>
       "(C = " <> exText[c] <> ")", Hold[exSumStepQ[t, k, a, n, c]]]}]];

exactChecks["HasDerivAt", ir_] := With[{ob = ir["Obligation"]},
  With[{F = ob["F"], f = ob["f"], x = ob["Variable"]},
    {exCheck["Together[D[F, x] - f] === 0 for F = " <> exText[F], Hold[exDerivQ[F, f, x]]]}]];

exactChecks[_, ir_] := {exUndecided["no exact check for this certificate kind"]};

(* ---------------------------------------------------------------------- *)
(* The "Exact" backend: runs, trusted steps and presentation               *)
(* ---------------------------------------------------------------------- *)

(* One verdict per item: every check passed, Missing["Exact", True]; one
   failed, Missing["Exact", False] (Certify::exactfail); otherwise
   Missing["NotChecked"].  Each check runs quietly (a 0/0 in a refuted
   claim is a False, not a message) under
   TimeConstrained[..., Min[TimeConstraint, 10]]; a timeout is not decided. *)
exactRun[items_List, o_] := Module[{tl = Min[o["TimeConstraint"], 10]},
  Map[Function[item, Module[{cs = item["G"]["Checks"], rs, logs = {}, t, v},
    If[!TrueQ[o["Verify"]],
      Return[<|"Verified" -> Missing["NotChecked"], "Timing" -> 0., "Logs" -> {},
        "Cached" -> False, "Checks" -> ((#["Check"] -> Missing["NotChecked"]) & /@ cs)|>,
        Module]];
    t = First[AbsoluteTiming[rs = Map[Function[c, With[{r = TimeConstrained[
          Quiet[ReleaseHold[c["Test"]]], tl, $exactTimeout]},
        Which[BooleanQ[r], r,
          r === $exactTimeout, AppendTo[logs, "timed out after " <> ToString[tl] <>
            " s: " <> c["Check"]]; Missing["NotChecked"],
          True, Missing["NotChecked"]]]], cs]]];
    v = Which[MemberQ[rs, False], Missing["Exact", False],
      rs =!= {} && AllTrue[rs, TrueQ], Missing["Exact", True],
      True, Missing["NotChecked"]];
    If[v === Missing["Exact", False], Message[Certify::exactfail, item["G"]["Input"]]];
    <|"Verified" -> v, "Timing" -> t, "Logs" -> logs, "Cached" -> False,
      "Checks" -> Thread[(#["Check"] & /@ cs) -> rs]|>]], items]];

$exactTrust = "Mathilda's exact rational arithmetic and polynomial Expand: an exact \
re-check by the system that produced the certificate, not a kernel proof";
$exactTypeTrust = <|
  "Primality" -> "Lucas's primality theorem (Mathilda replays the Pratt certificate with \
PowerMod; 2 is the base of the scheme)",
  "Induction" -> "the induction principle over the integers n >= the lower limit (Mathilda \
checks the base and the step)",
  "HasDerivAt" -> "Mathilda's D and Together (the derivative identity holds where both \
sides are defined)",
  "RootsComplete" -> "no zero divisors: a nonzero constant times linear factors vanishes \
exactly at their roots",
  "ClosedFact" -> "the meanings of the statement's matrix operations as Mathilda's \
Dimensions, Dot, IdentityMatrix, Length and a cofactor expansion of det"|>;

(* Exact's own steps; the Lean-specific Lucas step of the Pratt generators
   ("Lean evaluates Mathilda.prattCheck") is replaced by Exact's. *)
exactSteps[g_] := DeleteDuplicates[Join[DeleteCases[g["TrustedSteps"], $lucasStep],
  exactTrusted[g]]];
exactTrusted[g_] := Join[{$exactTrust},
  With[{s = Lookup[$exactTypeTrust, g["Obligation"]["Type"], None]}, If[StringQ[s], {s}, {}]]];

exactPresent[item_, v_, o_] := With[{g = item["G"]}, With[{a = <|
  "Kind" -> g["Kind"], "Input" -> g["Input"], "Result" -> g["Result"],
  "Theorem" -> g["Theorem"], "Certificate" -> g["Certificate"], "Checks" -> v["Checks"],
  "Verified" -> v["Verified"], "Engine" -> "Exact", "Tier" -> "Exact",
  "Domain" -> g["Domain"], "Hypotheses" -> g["Hypotheses"],
  "TrustedSteps" -> exactSteps[g],
  "Timing" -> v["Timing"], "Cached" -> False|>},
  $issued[Hash[a]] = True;
  ProofObject[a]]];

(* A refutation (Refute.m): the witness re-checked exactly (domain,
   denominators, hypotheses True and conclusion False), under the same time
   limit as the other checks; Missing["Exact", True] when it holds, else
   Missing["NotChecked"] (never a refutation on a doubtful check). *)
exactRefuteCheck[rir_, o_] := Module[{text, ok, t},
  text = "the statement is False at " <> ToString[rir["Witness"], InputForm] <>
    " (every hypothesis True, the conclusion False, every denominator nonzero, the \
witness in the statement's domain), by exact arithmetic";
  If[!TrueQ[o["Verify"]], Return[<|"Verified" -> Missing["NotChecked"], "Timing" -> 0.,
    "Logs" -> {}, "Cached" -> False, "Checks" -> {text -> Missing["NotChecked"]}|>, Module]];
  t = First[AbsoluteTiming[ok = TimeConstrained[witnessOK[rir["Problem"], rir["Witness"]],
    Min[o["TimeConstraint"], 10], $exactTimeout]]];
  <|"Verified" -> If[ok === True, Missing["Exact", True], Missing["NotChecked"]],
    "Timing" -> t, "Logs" -> {}, "Cached" -> False,
    "Checks" -> {text -> If[BooleanQ[ok], ok, Missing["NotChecked"]]}|>];

registerBackend["Exact", <|"Name" -> "Exact", "Trust" -> "ExactArithmetic",
  "Encode" -> exactEncode, "Run" -> exactRun, "Present" -> exactPresent, "Refute" -> exactRefuteCheck,
  "Label" -> ("Exact" &), "TrustedSteps" -> exactTrusted,
  "Options" -> {"TimeConstraint"}|>];

End[];

EndPackage[];
