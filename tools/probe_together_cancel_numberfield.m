(* ::Package:: *)
(* tools/probe_together_cancel_numberfield.m
 *
 * Probe of the rational-function and factoring builtins over number fields:
 * Together, Cancel (with and without Extension -> Automatic), Apart (with and
 * without Extension -> Automatic), Factor (default, Extension -> Automatic,
 * Extension -> {gens}) and FactorList.
 *
 * Fields: Q(I), Q(Sqrt[2], Sqrt[3]), Q(I, Sqrt[3]), Q(2^(1/3)), Q(Root of a
 * cubic), Q((-1)^(1/3)), and a parametric one (symbols a, k with the symbolic
 * radical Sqrt[k - 1]).  1-3 variables; rational inputs from ~100 to ~5000
 * leaves, built as random nested sums / products / quotients of products of a
 * small pool of factors (so that real cancellations occur); polynomial inputs
 * are expanded products of known factors over the field.
 *
 * Usage, from the repository root:
 *
 *   ./Mathilda -file tools/probe_together_cancel_numberfield.m
 *       generates the cases (once, into $probeDir/cases.m), runs every
 *       operation, checks every result, writes $probeDir/mathilda.m, and, if
 *       $probeDir/reference.m exists, prints the comparison table.
 *   wolframscript -file tools/probe_together_cancel_numberfield.m
 *       runs the same cases in the reference system (timing and LeafCount
 *       only -- the reference is never used as a correctness oracle) and
 *       writes $probeDir/reference.m.  Run it after the first Mathilda run,
 *       then run Mathilda again for the table.
 *
 * $probeDir defaults to "probe_nf_out" (relative to the working directory);
 * set probeDir before Get-ing this file to override.
 *
 * Correctness (Mathilda side only): Together / Cancel / Apart results are
 * compared with the input numerically, at two random rational points, to 40
 * digits; Factor / FactorList results must Expand back to the input exactly.
 * The factor count is compared with the reference's.
 *
 * A row is flagged when Mathilda is more than 5x slower (and slower than
 * 10 ms) or produces more than 5x the reference's LeafCount, when a check
 * fails, when an operation times out, or when a factor count differs. *)

If[!ValueQ[probeDir], probeDir = "probe_nf_out"];
If[!ValueQ[probeTimeout], probeTimeout = 20];
(* optional case-id range, to split a run over several processes (each with
   its own probeDir holding a copy of cases.m; merge the result files) *)
If[!ValueQ[probeFirst], probeFirst = 1];
If[!ValueQ[probeLast], probeLast = Infinity];
isMathilda = StringContainsQ[$Version, "Mathilda"];
casesFile = probeDir <> "/cases.m";
outFile = probeDir <> "/" <> If[isMathilda, "mathilda.m", "reference.m"];
refFile = probeDir <> "/reference.m";
If[!DirectoryQ[probeDir], CreateDirectory[probeDir]];

(* ---------------------------------------------------------------- *)
(* Case generation (Mathilda only; the reference reads the file)   *)
(* ---------------------------------------------------------------- *)

fields = {
  {"QI",    {I},                         {2},    {}},
  {"Q23",   {Sqrt[2], Sqrt[3]},          {2, 2}, {}},
  {"QI3",   {I, Sqrt[3]},                {2, 2}, {}},
  {"Q2c",   {2^(1/3)},                   {3},    {}},
  {"QRoot", {Root[#^3 - # - 1 &, 1]},    {3},    {}},
  {"Qw",    {(-1)^(1/3)},                {2},    {}},
  {"Param", {Sqrt[k - 1]},               {2},    {a}}
};

(* the power-product basis of the field (parameters multiply in) *)
fieldBasis[gens_, degs_, params_] := Module[{b = {1}},
  Do[b = Flatten[Table[m gens[[i]]^j, {m, b}, {j, 0, degs[[i]] - 1}]], {i, Length[gens]}];
  Do[b = Join[b, p b], {p, params}];
  b];

randElt[basis_] := Module[{n = RandomInteger[{1, Min[3, Length[basis]]}], picks},
  picks = RandomSample[basis, n];
  Total[Table[RandomChoice[{-3, -2, -1, 1, 2, 3}] m, {m, picks}]]];

randNonzeroElt[basis_] := Module[{c = randElt[basis]}, While[c === 0, c = randElt[basis]]; c];

(* a pool of factors over the field: linear and quadratic in the first
   variable, the other variables entering linearly *)
factorPool[vars_, basis_, n_] := Table[
  Module[{v = First[vars], rest = Rest[vars], lin},
    lin = Total[Table[RandomInteger[{-1, 1}] w, {w, rest}]];
    If[RandomInteger[{0, 2}] > 0,
      v + randElt[basis] + lin,
      v^2 + randElt[basis] v + randNonzeroElt[basis] + lin]],
  {n}];

genPoly[pool_, basis_] := Module[{f},
  f = randNonzeroElt[basis] Times @@ RandomChoice[pool, RandomInteger[{1, 3}]];
  If[RandomInteger[] == 1, Expand[f], f]];

genRat[pool_, basis_, 0] :=
  If[RandomInteger[{0, 3}] > 0, genPoly[pool, basis]/genPoly[pool, basis], genPoly[pool, basis]];
genRat[pool_, basis_, d_] := Module[{c = RandomInteger[{1, 4}]},
  Which[
    c == 1, Total[Table[genRat[pool, basis, d - 1], {RandomInteger[{2, 3}]}]],
    c == 2, genRat[pool, basis, d - 1] genRat[pool, basis, d - 1],
    c == 3, genRat[pool, basis, d - 1]/genRat[pool, basis, d - 1],
    True,   genRat[pool, basis, d - 1]^2 + genRat[pool, basis, d - 1]]];

sizeBuckets = {{100, 400}, {400, 1500}, {1500, 5000}};

genRatCase[basis_, vars_, {lo_, hi_}] := Module[{e = Null, tries = 0, d = 1, lc},
  While[e === Null && tries < 400,
    tries++;
    d = Which[lo < 300, RandomInteger[{1, 2}], lo < 1000, RandomInteger[{2, 3}], True, RandomInteger[{3, 4}]];
    Module[{cand = genRat[factorPool[vars, basis, RandomInteger[{3, 5}]], basis, d]},
      lc = LeafCount[cand];
      If[lo <= lc <= hi && !FreeQ[cand, Power[_, _?Negative]] &&
         Length[Intersection[Variables[cand], vars]] == Length[vars], e = cand]]];
  e];

genPolyCase[basis_, vars_, {lo_, hi_}] := Module[{p = Null, tries = 0, pool, fs, lc},
  While[p === Null && tries < 400,
    tries++;
    pool = factorPool[vars, basis, RandomInteger[{2, 4}]];
    fs = Table[f^RandomInteger[{1, 2}], {f, pool}];
    Module[{cand = Expand[randNonzeroElt[basis] Times @@ fs]},
      lc = LeafCount[cand];
      If[lo <= lc <= hi && Length[Intersection[Variables[cand], vars]] == Length[vars],
        p = {cand, pool}]]];
  p];

makeCases[] := Module[{cases = {}, id = 0},
  SeedRandom[20261127];
  Do[
    Module[{name = fld[[1]], gens = fld[[2]], degs = fld[[3]], params = fld[[4]], basis},
      basis = fieldBasis[gens, degs, params];
      Do[
        Module[{vars = Take[{x, y, z}, nv]},
          Do[
            Module[{e = genRatCase[basis, vars, b]},
              If[e =!= Null,
                id++;
                AppendTo[cases, {id, name, nv, "rat", e, vars, gens}]]],
            {b, sizeBuckets}];
          Do[
            Module[{pp = genPolyCase[basis, vars, b]},
              If[pp =!= Null,
                id++;
                AppendTo[cases, {id, name, nv, "poly", pp[[1]], vars, gens}]]],
            {b, {{60, 400}, {400, 2000}}}]],
        {nv, 1, 3}]],
    {fld, fields}];
  cases];

If[!FileExistsQ[casesFile],
  If[!isMathilda, Print["no cases file; run the probe in Mathilda first"]; Quit[]];
  Print["generating cases ..."];
  Put[makeCases[], casesFile]];
cases = Get[casesFile];
Print[Length[cases], " cases"];

(* ---------------------------------------------------------------- *)
(* Operations                                                       *)
(* ---------------------------------------------------------------- *)

ratOps = {"Together", "Cancel", "TogetherExt", "CancelExt", "Apart", "ApartExt"};
polyOps = {"Factor", "FactorExt", "FactorGens", "FactorList"};

runOp["Together", e_, v_, g_] := Together[e];
runOp["Cancel", e_, v_, g_] := Cancel[e];
runOp["TogetherExt", e_, v_, g_] := Together[e, Extension -> Automatic];
runOp["CancelExt", e_, v_, g_] := Cancel[e, Extension -> Automatic];
runOp["Apart", e_, v_, g_] := Apart[e, First[v]];
runOp["ApartExt", e_, v_, g_] := Apart[e, First[v], Extension -> Automatic];
runOp["Factor", e_, v_, g_] := Factor[e];
runOp["FactorExt", e_, v_, g_] := Factor[e, Extension -> Automatic];
runOp["FactorGens", e_, v_, g_] := Factor[e, Extension -> g];
runOp["FactorList", e_, v_, g_] := FactorList[e];

(* non-constant factor bases of a product *)
factorCount[r_] := Module[{fs = If[Head[r] === Times, List @@ r, {r}]},
  Length[Select[fs /. Power[b_, _Integer?Positive] :> b, !NumericQ[#] &]]];
factorCountList[r_List] := Length[Select[r, !NumericQ[First[#]] &]];
factorCountList[_] := -1;

(* numeric value check at random rational points (Mathilda only) *)
randPoint[vars_, params_] := Join[
  Thread[vars -> Table[RandomInteger[{-40, 40}]/RandomInteger[{7, 29}], {Length[vars]}]],
  Thread[params -> Table[RandomInteger[{30, 90}]/RandomInteger[{7, 11}], {Length[params]}]]];

numVal[e_, pt_] := Quiet[N[e /. pt, 50]];

goodNum[z_] := NumericQ[z] && FreeQ[z, Indeterminate | ComplexInfinity | DirectedInfinity];

checkValue[e_, r_, vars_, params_] := Module[{ok = True, tested = 0, tries = 0, pt, ve, vr},
  While[tested < 2 && tries < 12,
    tries++;
    pt = randPoint[vars, params];
    ve = numVal[e, pt];
    If[goodNum[ve],
      vr = numVal[r, pt];
      tested++;
      If[!goodNum[vr] || Abs[ve - vr] > 10^-25 Max[1, Abs[ve]], ok = False]]];
  If[tested == 0, "untested", ok]];

(* Expand[r - p] === 0 is exact but not decisive: Expand leaves
   (-1 + k)^(3/2) - (-1 + k) Sqrt[-1 + k] and (-1)^(1/3) - (-1)^(2/3) - 1 as
   they are.  When it is not 0 the difference is checked numerically too. *)
checkFactor[p_, r_, vars_] := Quiet[Expand[r - p] === 0] ||
  TrueQ[checkValue[p, r, Join[vars, {a, k}], {}]];
checkFactorList[p_, r_List, vars_] := checkFactor[p, Times @@ (#[[1]]^#[[2]] & /@ r), vars];
checkFactorList[p_, _, _] := False;

paramsOf[c_] := If[c[[2]] === "Param", {a, k}, {}];

(* ---------------------------------------------------------------- *)
(* Run                                                              *)
(* ---------------------------------------------------------------- *)

results = If[FileExistsQ[outFile], Get[outFile], {}];
done = Map[{#[[1]], #[[2]]} &, results];

Do[
  Module[{id = c[[1]], kind = c[[4]], e = c[[5]], vars = c[[6]], gens = c[[7]], ops, row},
    ops = If[kind === "rat", ratOps, polyOps];
    Do[
      If[!MemberQ[done, {id, op}] && probeFirst <= id <= probeLast,
        Module[{t, r, ok = "n/a", nf = -1, lc},
          {t, r} = AbsoluteTiming[TimeConstrained[Quiet[runOp[op, e, vars, gens]], probeTimeout, $probeTimedOut]];
          If[r === $probeTimedOut,
            t = "TO"; lc = -1,
            lc = LeafCount[r];
            (* the reference's Apart has no Extension option: it returns the
               call unevaluated (Apart::optx); such a row is not compared *)
            If[op === "ApartExt" && Head[r] === Apart, lc = -2];
            nf = Which[op === "FactorList", factorCountList[r], kind === "poly", factorCount[r], True, -1];
            If[isMathilda,
              ok = Which[
                op === "FactorList", checkFactorList[e, r, vars],
                kind === "poly", checkFactor[e, r, vars],
                True, checkValue[e, r, Join[vars, Complement[paramsOf[c], vars]], {}]]]];
          row = {id, op, t, lc, ok, nf, LeafCount[e]};
          AppendTo[results, row];
          Put[results, outFile];
          If[ok === False || t === "TO", Print["  ", row]]]],
      {op, ops}]],
  {c, cases}];

Print["wrote ", outFile, " (", Length[results], " rows)"];

(* ---------------------------------------------------------------- *)
(* Comparison (Mathilda only, when the reference has been run)      *)
(* ---------------------------------------------------------------- *)

If[isMathilda && FileExistsQ[refFile],
  Module[{ref = Get[refFile], refT, flagged = {}, byOp},
    refT = Association[Map[({#[[1]], #[[2]]} -> #) &, ref]];
    Do[
      Module[{key = {row[[1]], row[[2]]}, rr, tm, tr, lm, lr, why = {}},
        If[KeyExistsQ[refT, key] && refT[key][[4]] =!= -2,
          rr = refT[key];
          tm = row[[3]]; tr = rr[[3]]; lm = row[[4]]; lr = rr[[4]];
          If[tm === "TO" && tr =!= "TO", AppendTo[why, "timeout"]];
          If[NumericQ[tm] && NumericQ[tr] && tm > 0.01 && tm > 5 Max[tr, 0.0005], AppendTo[why, "slow"]];
          If[lm > 0 && lr > 0 && lm > 5 lr, AppendTo[why, "large"]];
          If[row[[5]] === False, AppendTo[why, "WRONG"]];
          If[row[[6]] >= 0 && rr[[6]] >= 0 && row[[6]] =!= rr[[6]], AppendTo[why, "nfactors"]];
          If[why =!= {}, AppendTo[flagged, Join[row, {tr, lr, rr[[6]], why}]]]]],
      {row, results}];
    Print["\n=== summary by operation (Mathilda vs reference) ==="];
    byOp = GatherBy[results, #[[2]] &];
    Do[
      Module[{op = g[[1, 2]], n = Length[g], tms, trs, ratios},
        tms = Select[g, NumericQ[#[[3]]] &];
        ratios = Select[Map[Module[{k = {#[[1]], #[[2]]}},
              If[KeyExistsQ[refT, k] && refT[k][[4]] =!= -2 && NumericQ[refT[k][[3]]] && NumericQ[#[[3]]],
                (#[[3]] + 0.0005)/(refT[k][[3]] + 0.0005), Null]] &, g], NumericQ];
        Print[op, ": ", n, " cases, ",
          Count[g, r_ /; r[[5]] === False], " wrong, ",
          Count[g, r_ /; r[[3]] === "TO"], " timeouts, ",
          "median time ratio ", If[ratios === {}, "-", N[Median[ratios], 3]],
          ", max ", If[ratios === {}, "-", N[Max[ratios], 3]],
          ", flagged ", Count[flagged, f_ /; f[[2]] === op]]],
      {g, byOp}];
    Print["\n=== flagged rows: {id, op, t, leaf, ok, nfac, inleaf, ref t, ref leaf, ref nfac, why} ==="];
    Do[Print[f], {f, flagged}];
    Put[flagged, probeDir <> "/flagged.m"]]];
