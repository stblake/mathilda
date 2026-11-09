(* tools/gr/verify.m -- the numeric oracle for the Gradshteyn & Ryzhik corpus,
   in Mathilda.

   Loaded by tools/gr/verify.py, which writes one call per transcribed record:

       grVerify[id, lhs, {rhs1, rhs2, ...}, {cond1, cond2, ...}]

   (all four as InputForm strings, exactly as the transcription batches hold
   them; see tools/gr/TRANSCRIBE.md) and runs the file with ./Mathilda -file.
   The oracle decides whether the book's value -- as read -- is the value of
   the integral.  It is the only source of the reference values stored in the
   corpus, so a misreading or a misprint in the book never becomes an expected
   answer.

   Per record:
    1. Parse lhs/rhs/conds held (a syntax error is BADEXPR, not a crash).
    2. Expand elliptic notation placeholders (MatharK[k] etc., G&R's modulus
       convention) into the four candidate conventions; keep the first variant
       that verifies.
    3. Isolate the integral: lhs may be c*Integrate[f, spec] + d; it is solved
       for the integral, which itself is never evaluated.
    4. Draw parameter instantiations: parameters the conditions declare integer
       from the integer pool; conventionally-integer names (n, m, k, l, j) as
       integers first; the rest from a pool of positive rationals.  An
       instantiation must satisfy every condition.
    5. Definite: compare N[value, 20] with NIntegrate at two instantiations.
         strict  NIntegrate finished without a single message and
                 |num - val| / max(1, |val|) < 1e-7;
         loose   NIntegrate (default or LevinRule) only reached its answer with
                 a convergence message, and agrees to 1e-5.  An oscillatory
                 kernel over [0, Infinity) -- most of chapter 6 -- often ends
                 here; a misreading is still rejected (it is off by O(1)), and
                 the record is flagged so the corpus can tell the two apart.
       A value below 1e-6 proves nothing under a relative test; skipped.
       Indefinite: D[F, x] - f must vanish (relative 1e-8) at x = 0.31, 0.77,
       1.43 for each instantiation -- differentiation, independent of
       NIntegrate entirely.
    6. Reference value N[value, 30] per instantiation.

   Output, one tab-separated line per fact (parsed by verify.py):
     @@R  id status kind var integrand limits value params variant why
     @@G  id inst ref num tier inst_integrand inst_limits
     @@F  id inst f var lo hi value N[value] verdict   (FullForm probes)
   status: VERIFIED (2 instantiations), VERIFIED1 (only one admissible),
           MISMATCH, NOCONV, NOHEAD (uses a function Mathilda lacks), BADEXPR,
           TIMEOUT.

   All helpers are gr-prefixed and every local is a Module variable: the
   records' own symbols (f, x, a, ...) live in Global` and must never meet a
   definition made here. *)

grIntNames = {"n", "m", "k", "l", "j"};
grRealPool = {3/8, 5/8, 5/6, 5/4, 7/4, 7/3, 1/5, 3};
grIntPool  = {1, 2, 3, 4};

(* G&R writes elliptic integrals with the modulus k; Mathilda's heads take the
   parameter m = k^2.  The transcription uses placeholder heads and the oracle
   tries each convention (modulus amplitude-first, modulus modulus-first, then
   the parameter convention twice), as for the Mathar corpus. *)
grVariants = {
  {MatharK[grK_] :> EllipticK[grK^2], MatharE[grK_] :> EllipticE[grK^2],
   MatharK1[grK_] :> EllipticK[1 - grK^2], MatharE1[grK_] :> EllipticE[1 - grK^2],
   MatharF2[grA_, grB_] :> EllipticF[grA, grB^2], MatharE2[grA_, grB_] :> EllipticE[grA, grB^2],
   MatharPi[grN_, grK_] :> EllipticPi[grN, grK^2]},
  {MatharK[grK_] :> EllipticK[grK^2], MatharE[grK_] :> EllipticE[grK^2],
   MatharK1[grK_] :> EllipticK[1 - grK^2], MatharE1[grK_] :> EllipticE[1 - grK^2],
   MatharF2[grA_, grB_] :> EllipticF[grB, grA^2], MatharE2[grA_, grB_] :> EllipticE[grB, grA^2],
   MatharPi[grN_, grK_] :> EllipticPi[grN^2, grK^2]},
  {MatharK[grK_] :> EllipticK[grK], MatharE[grK_] :> EllipticE[grK],
   MatharK1[grK_] :> EllipticK[1 - grK], MatharE1[grK_] :> EllipticE[1 - grK],
   MatharF2[grA_, grB_] :> EllipticF[grA, grB], MatharE2[grA_, grB_] :> EllipticE[grA, grB],
   MatharPi[grN_, grK_] :> EllipticPi[grN, grK]},
  {MatharK[grK_] :> EllipticK[grK], MatharE[grK_] :> EllipticE[grK],
   MatharK1[grK_] :> EllipticK[1 - grK], MatharE1[grK_] :> EllipticE[1 - grK],
   MatharF2[grA_, grB_] :> EllipticF[grB, grA], MatharE2[grA_, grB_] :> EllipticE[grB, grA],
   MatharPi[grN_, grK_] :> EllipticPi[grN, grK]}};

grHasPlaceholder[e_] := ! FreeQ[e, MatharK | MatharE | MatharK1 | MatharE1 | MatharF2 | MatharE2 | MatharPi];

(* Free Global` symbol names of an expression, without the dummy indices of
   Sum / Product / Table. *)
grSymNames[e_] := Complement[
  Union[Cases[e, grS_Symbol /; Context[grS] === "Global`" :> SymbolName[grS], {0, Infinity}, Heads -> True]],
  Union[Flatten[Cases[e, (Sum | Product | Table)[_, grIts__List] :>
      Cases[{grIts}, {grV_Symbol, ___} :> SymbolName[grV]], {0, Infinity}, Heads -> True]]]];

grRelErr[a_, b_] := Abs[a - b]/Max[1, Abs[b]];
grFiniteQ[v_] := NumericQ[v] && Abs[v] < 10^30;

(* Mathilda's NIntegrate can return a wrong number WITHOUT a message (seen on
   oscillatory Bessel products over [0, Infinity): 0.63 for the 1/(2 beta) =
   0.286 of G&R 6.512.3).  So no single NIntegrate run may reject a reading:
     grNIntA  the default strategy; {value, "strict"} when it finished without a
              message, {value, "loose"} when only with one, else {$Failed, ..};
     grNIntB  an independent second strategy -- LevinRule on an infinite range
              with an oscillatory kernel, else a 20-digit MPFR run.
   grJudge (below) verifies on whichever agrees with the book, and calls a
   MISMATCH only when the two strategies agree with EACH OTHER and not with the
   book; when they disagree with each other the instantiation is inconclusive. *)
grOscHeads = BesselJ | BesselY | Sin | Cos | StruveH | AngerJ | WeberE | Exp | Power;
grNIntA[f_, x_, lo_, hi_] := Module[{r},
  r = Check[TimeConstrained[NIntegrate[f, {x, lo, hi}], 25, $Failed], $Failed];
  If[NumericQ[r], Return[{r, "strict"}]];
  r = Quiet[TimeConstrained[NIntegrate[f, {x, lo, hi}], 25, $Failed]];
  If[NumericQ[r], {r, "loose"}, {$Failed, "none"}]];
grNIntB[f_, x_, lo_, hi_] := Module[{r},
  r = If[(hi === Infinity || lo === -Infinity) && ! FreeQ[f, BesselJ | BesselY | Sin | Cos | StruveH | AngerJ | WeberE | Complex],
    Quiet[TimeConstrained[NIntegrate[f, {x, lo, hi}, Method -> "LevinRule"], 25, $Failed]],
    Quiet[TimeConstrained[NIntegrate[f, {x, lo, hi}, WorkingPrecision -> 20], 40, $Failed]]];
  If[NumericQ[r], r, $Failed]];

(* One instantiation of a definite record: {"good", num, tier} | {"bad", num} |
   {"skip"}.  tier: strict (default strategy, no message, 1e-7), loose (1e-5,
   after a message or from the second strategy), loose-silentwrong (the
   default strategy returned a different number WITHOUT a message: an
   NIntegrate wrong answer, reported). *)
grJudge[f_, x_, lo_, hi_, vi_] := Module[{a, b},
  a = grNIntA[f, x, lo, hi];
  If[NumericQ[a[[1]]] && Abs[a[[1]]] >= 10^-6,
    If[a[[2]] === "strict" && grRelErr[a[[1]], vi] < 10^-7, Return[{"good", a[[1]], "strict"}]];
    If[grRelErr[a[[1]], vi] < 10^-5, Return[{"good", a[[1]], "loose"}]]];
  b = grNIntB[f, x, lo, hi];
  If[! NumericQ[b] || Abs[b] < 10^-6, Return[If[NumericQ[a[[1]]] && a[[2]] === "strict", {"bad", a[[1]]}, {"skip"}]]];
  If[grRelErr[b, vi] < 10^-5,
    Return[{"good", b, If[a[[2]] === "strict", "loose-silentwrong", "loose"]}]];
  If[NumericQ[a[[1]]] && grRelErr[a[[1]], b] < 10^-5, {"bad", b}, {"skip"}]];

grIntDeclared[conds_] := Union[Cases[conds, Element[grS_Symbol, Integers] :> SymbolName[grS], {0, Infinity}]];

grInstances[params_, conds_, seed_] := Module[{decl = grIntDeclared[conds], tries = {}},
  SeedRandom[seed];
  Do[
    AppendTo[tries, Table[
      Symbol[grP] -> Which[
        MemberQ[decl, grP], RandomChoice[{1, 2, 3, 4, 5, 6}],
        MemberQ[grIntNames, grP] && grI <= 24, RandomChoice[grIntPool],
        grI <= 40, RandomChoice[grRealPool],
        (* late draws admit negatives, for conditions such as a b < 0 *)
        True, RandomChoice[Join[grRealPool, -grRealPool]]], {grP, params}]],
    {grI, 64}];
  DeleteDuplicates[tries]];

grCondsHold[conds_, inst_] := AllTrue[conds, TrueQ[# /. grInst] &];

(* {Hold[f, spec], Hold[value]} for lhs == rhs, or $Failed.  The book's value
   stays HELD: it is only ever released after numbers have been substituted
   (grRel).  Evaluating it with symbolic parameters sent finite Sums with
   symbolic limits into closed-form summation / Simplify (G&R 2.673: > 25 min,
   not interruptible by TimeConstrained).  lhs may be c*Integrate[..] + d with
   c, d free of the integral (evaluated: they come from the lhs only). *)
grIsolate[lhsH_, rhsH_] := Module[{ints, l, cc, dd},
  ints = Cases[lhsH, HoldPattern[Integrate[grF_, grS_]] :> Hold[grF, grS], {0, Infinity}];
  If[Length[ints] =!= 1, Return[$Failed]];
  l = ReleaseHold[lhsH /. HoldPattern[Integrate[_, _]] -> grII];
  If[FreeQ[l, grII], Return[$Failed]];
  If[l === grII, Return[{ints[[1]], rhsH}]];
  l = Expand[l];
  cc = Coefficient[l, grII]; dd = l /. grII -> 0;
  If[! FreeQ[{cc, dd}, grII] || cc === 0, Return[$Failed]];
  {ints[[1]], With[{c1 = cc, d1 = dd}, rhsH /. Hold[grZ_] :> Hold[(grZ - d1)/c1]]}];

(* Release a held value after substituting numbers. *)
grRel[valH_, rules_] := ReleaseHold[valH /. rules];

grStr[e_] := ToString[e, InputForm];

(* Which failed attempt to report when no variant/reading verifies. *)
grRank["VERIFIED1"] = 5; grRank["NOHEAD"] = 4; grRank["MISMATCH"] = 3; grRank["BADEXPR"] = 2; grRank["NOCONV"] = 1;

grFF[e_] := ToString[FullForm[e]];
(* FullForm of e /. rules WITHOUT evaluating the result: an antiderivative with
   x still symbolic can send a finite Sum into closed-form summation that
   TimeConstrained cannot interrupt.  Python (pycheck.py) evaluates it. *)
grFFHeld[e_Hold, rules_] := Module[{s = ToString[FullForm[e /. rules]]},
  StringTake[s, {6, StringLength[s] - 1}]];
grFFHeld[e_, rules_] := grFFHeld[Hold @@ {e}, rules];
(* InputForm text of a held expression, without the Hold wrapper. *)
grStrHeld[e_Hold] := Module[{s = ToString[e, InputForm]}, StringTake[s, {6, StringLength[s] - 1}]];
grStrHeld[e_] := grStr[e];

(* One (variant, rhs) attempt.  Returns
   {status, kind, var, f, lim, val, params, good, why, probes}.
   probes: up to three admissible instantiations in FullForm, with what Mathilda
   made of them -- {inst, f, var, lo, hi, value, N[value], NIntegrate verdict} --
   for the independent Python cross-check (tools/gr/pycheck.py, mpmath+scipy),
   which re-decides every instantiation without Mathilda. *)
grTryOne[id_, lhsH_, rhsH_, conds_] := Module[
  {iso, f, spec, kind, var, lim, val, params, extra, nohead, insts, good = {}, bad = {},
   fi, vi, nv, xs, err, cand, fl = {}},
  iso = grIsolate[lhsH, rhsH];
  If[iso === $Failed, Return[{"BADEXPR", "", "", "", "", "", {}, {}, "isolate", {}}]];
  f = Extract[iso[[1]], 1]; spec = Extract[iso[[1]], 2]; val = iso[[2]];
  kind = If[ListQ[spec], "definite", "indefinite"];
  var = If[ListQ[spec], spec[[1]], spec];
  lim = If[ListQ[spec], Rest[spec], {}];
  params = DeleteCases[grSymNames[{f, lim}], SymbolName[var]];
  extra = Complement[grSymNames[val], params, {SymbolName[var]}];
  (* A Global` symbol used as a head is a function Mathilda does not define
     (ParabolicCylinderD, HankelH1, ...): a Mathilda gap, not a bad reading. *)
  nohead = Union[Cases[{f, val}, grH_Symbol[___] /; Context[grH] === "Global`" :> SymbolName[grH], {0, Infinity}, Heads -> True]];
  If[nohead =!= {}, Return[{"NOHEAD", kind, var, f, lim, val, params, {}, StringRiffle[nohead, ","], {}}]];
  If[extra =!= {}, Return[{"BADEXPR", kind, var, f, lim, val, params, {}, "value has extra symbols " <> grStr[extra], {}}]];
  insts = If[params === {}, {{}}, grInstances[params, conds, Total[ToCharacterCode[id]]]];
  Do[
    If[Length[good] >= 2, Break[]];
    If[conds =!= {} && Length[insts] > 1 && ! grCondsHold[conds, grInst], Continue[]];
    If[kind === "definite",
      vi = Quiet[N[grRel[val, grInst], 20]];
      (* the book's value does not numericize in Mathilda: keep the residue so
         the missing evaluation can be read off (a Mathilda numeric gap) *)
      If[Length[fl] < 3 && ! grFiniteQ[vi], AppendTo[fl, {grInst, grFF[f /. grInst], grFF[var],
          grFF[lim[[1]] /. grInst], grFF[lim[[2]] /. grInst], grFFHeld[val, grInst], "novalue", ""}]];
      If[! grFiniteQ[vi], AppendTo[bad, {grInst, "novalue", StringTake[grStr[vi], Min[160, StringLength[grStr[vi]]]]}]; Continue[]];
      If[Abs[vi] < 10^-6, Continue[]];
      nv = grJudge[f /. grInst, var, lim[[1]] /. grInst, lim[[2]] /. grInst, vi];
      If[Length[fl] < 3, AppendTo[fl, {grInst, grFF[f /. grInst], grFF[var],
          grFF[lim[[1]] /. grInst], grFF[lim[[2]] /. grInst], grFFHeld[val, grInst], grStr[vi], grStr[nv]}]];
      Switch[nv[[1]],
        "good", AppendTo[good, {grInst, Quiet[N[grRel[val, grInst], 30]], nv[[2]], nv[[3]]}],
        "bad", AppendTo[bad, {grInst, nv[[2]], vi}],
        _, Null],
      (* indefinite: the derivative test, as a 60-digit central difference
         with the parameters AND x substituted in one replacement, so the
         book's value is only ever evaluated at numbers.  Symbolic D (or
         substituting the parameters while x is still symbolic) sent a finite
         Sum with binomials into closed-form summation / Simplify / FLINT gcd
         that TimeConstrained could not interrupt (G&R 2.673, >25 min), and D
         does not thread over a Sum with a symbolic limit.  h = 10^-15 at 60
         digits leaves ~45 significant digits in the difference and an O(h^2)
         truncation error of ~10^-30. *)
      xs = {31/100, 77/100, 143/100};
      err = Quiet[Max[Table[
          Module[{h = 10^-15, xp, xm, x0, fv, dv},
            x0 = N[grXX, 60]; xp = N[grXX + h, 60]; xm = N[grXX - h, 60];
            dv = grRel[val /. Abs[grU_] :> Sqrt[grU^2], Join[grInst, {var -> xp}]] -
                 grRel[val /. Abs[grU_] :> Sqrt[grU^2], Join[grInst, {var -> xm}]];
            dv = dv/(2 h);
            fv = f /. Join[grInst, {var -> x0}];
            Abs[dv - fv]/Max[1, Abs[fv]]], {grXX, xs}]]];
      If[Length[fl] < 3, AppendTo[fl, {grInst, grFF[f /. grInst], grFF[var], "", "",
          grFFHeld[val, grInst], "", grStr[err]}]];
      If[! NumericQ[err], Continue[]];
      If[err < 10^-8, AppendTo[good, {grInst, "", "", "strict"}], AppendTo[bad, {grInst, err}]]],
    {grInst, insts}];
  Which[
    Length[good] >= Min[2, Max[1, Length[insts]]],
      {"VERIFIED", kind, var, f, lim, val, params, good, "", fl},
    good === {} && bad =!= {},
      {"MISMATCH", kind, var, f, lim, val, params, {}, grStr[Take[bad, Min[2, Length[bad]]]], fl},
    Length[good] == 1,
      {"VERIFIED1", kind, var, f, lim, val, params, good, "", fl},
    True, {"NOCONV", kind, var, f, lim, val, params, {}, "", fl}]];

grEmit[id_, r_, variant_] := Module[{st, kind, var, f, lim, val, params, good, why, fl},
  {st, kind, var, f, lim, val, params, good, why, fl} = r;
  Print["@@R\t", id, "\t", st, "\t", kind, "\t", grStr[var], "\t", grStr[f], "\t",
    If[lim === {}, "", grStr[lim]], "\t", grStrHeld[val], "\t", StringRiffle[params, ","], "\t",
    variant, "\t", why];
  Scan[Print["@@G\t", id, "\t", grStr[#[[1]]], "\t", grStr[#[[2]]], "\t", grStr[#[[3]]], "\t",
      #[[4]], "\t", grStr[f /. #[[1]]], "\t", If[lim === {}, "", grStr[lim /. #[[1]]]]] &, good];
  Scan[Print["@@F\t", id, "\t", grStr[#[[1]]], "\t", #[[2]], "\t", #[[3]], "\t", #[[4]], "\t",
      #[[5]], "\t", #[[6]], "\t", #[[7]], "\t", #[[8]]] &, fl]];

grVerify[id_String, lhs_String, rhs_List, condStrs_List] := Module[
  {lhsH, rhsL, conds, vrs, best = {"NOCONV", "", "", "", "", "", {}, {}, "", {}}, bestV = 0, out, done = False},
  lhsH = Quiet[ToExpression[lhs, InputForm, Hold]];
  If[lhsH === $Failed, Print["@@R\t", id, "\tBADEXPR\t\t\t\t\t\t\t0\tlhs"]; Return[Null]];
  rhsL = DeleteCases[Quiet[ToExpression[#, InputForm, Hold]] & /@ rhs, $Failed];
  If[rhsL === {}, Print["@@R\t", id, "\tBADEXPR\t\t\t\t\t\t\t0\trhs"]; Return[Null]];
  conds = DeleteCases[Quiet[ToExpression[#, InputForm]] & /@ condStrs, $Failed];
  vrs = If[grHasPlaceholder[{lhsH, rhsL}], grVariants, {{}}];
  out = TimeConstrained[
    Catch[
      Do[
        Do[
          out = grTryOne[id, lhsH /. vrs[[grVi]], rhsL[[grRi]] /. vrs[[grVi]], conds];
          If[out[[1]] === "VERIFIED", Throw[{out, grVi}]];
          (* bestV == 0: keep the first real attempt even when it is NOCONV,
             so its FullForm probes still reach the Python cross-check *)
          If[bestV == 0 || grRank[out[[1]]] > grRank[best[[1]]], best = out; bestV = grVi],
          {grRi, Length[rhsL]}],
        {grVi, Length[vrs]}];
      {best, bestV}],
    300, $Failed];
  If[out === $Failed,
    Print["@@R\t", id, "\tTIMEOUT\t\t\t\t\t\t\t0\t"],
    grEmit[id, out[[1]], out[[2]]]];
  Null];
