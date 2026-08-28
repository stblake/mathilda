(* reduce_check_prelude.m
 *
 * Sample-point EQUIVALENCE + WITNESS verifier for the Reduce[] stress corpus
 * (tests/reduce_corpus.m).  Loaded once by the C runner
 * (tests/test_reduce_corpus.c) before it forks per case, and usable by a
 * standalone driver.
 *
 * Reduce returns a quantifier-free LOGICAL FORMULA (an And/Or/Not tree of
 * relational atoms), not a list of solution rules -- so verification is by
 * SAMPLING, never by string comparison, and is therefore immune to output
 * ordering, radical/complex spelling, and non-minimal-but-sound forms.
 *
 * The contract a case must satisfy, keyed on its `expected` tag:
 *   "decline"  Reduce must bubble back unevaluated (MatchQ[_Reduce]).  This is
 *              the soundness invariant: an engine that cannot decide an
 *              input must NOT invent a formula for it.
 *   True/False Reduce must decide to exactly that, AND the input must have that
 *              constant truth at every grid point (a soundness cross-check).
 *   "solved"   Reduce must produce a formula that is EQUIVALENT to the input:
 *                (a) grid check   -- at every point of a bounded rational/integer
 *                    grid over the variables and free parameters where BOTH the
 *                    input and the output decide to a Boolean, they must agree
 *                    (catches unsoundness AND incompleteness on the grid);
 *                (b) witness check -- each fully-determined solution branch of
 *                    the output (x == value, with parameters/generated C[k]
 *                    instantiated) must satisfy the input.  This verifies the
 *                    equality branches whose solutions (radicals, complex, off-
 *                    grid) the grid never lands on.
 *              A "solved" case with zero judged points fails (no vacuous pass).
 *
 * reduceVerdict returns {code, why} with code 0 = PASS, 1 = FAIL,
 * 2 = UNEVALUATED (bubbled back when a result was expected).  reduceCheckCode
 * returns just the integer (used by the C runner); reduceReport prints a TSV
 * line (used by a standalone driver).
 *)

(* Symbols that Cases would otherwise mistake for a free variable: domain
 * names / constants, and the relational operators that appear as ELEMENTS
 * (not heads) inside a chained Inequality[a, Less, x, Less, b]. *)
$rdConstants = {Reals, Integers, Rationals, Complexes, GaussianIntegers,
                Booleans, Automatic, Pi, E,
                Less, Greater, LessEqual, GreaterEqual, Equal, Unequal,
                And, Or, Not, Xor, Implies, Inequality, Element};

(* Leaf variables/parameters of a (possibly relational/logical) expression.
 * Variables[] does NOT descend into Greater/Equal/And/..., so use Cases at all
 * levels (heads excluded by default) and drop the domain constants. *)
rdFreeSyms[e_] := Complement[Union[Cases[e, _Symbol, {0, Infinity}]], $rdConstants];

(* Per-symbol sample set, sized by symbol count so the Tuples product stays
 * bounded: 1 -> 11, 2 -> 81, 3 -> 125, >=4 -> 3^n points. *)
rdSamples[dom_, nsyms_] := Which[
  nsyms <= 1, If[dom === Integers, {-4, -3, -2, -1, 0, 1, 2, 3, 4},
                 {-3, -2, -3/2, -1, -1/2, 0, 1/2, 1, 3/2, 2, 3}],
  nsyms == 2, If[dom === Integers, {-4, -2, -1, 0, 1, 2, 4},
                 {-3, -2, -1, -1/2, 0, 1/2, 1, 2, 3}],
  nsyms == 3, If[dom === Integers, {-3, -1, 0, 1, 3}, {-2, -1, 0, 1, 2}],
  True,       {-1, 0, 1}
];

(* Tri-state truth of a substituted relation/formula:
 *   1 = True, 0 = False, -1 = undecided (a free symbol survives, a division by
 *   zero produced ComplexInfinity, a Root inequality did not resolve, ...). *)
rdTruth[b_] := Which[TrueQ[b], 1, b === False, 0, True, -1];

(* A relation operand is OUT OF DOMAIN at pt when it evaluates to an infinite /
 * indeterminate value (a pole -- 1/x at x==0), or, over an ordered domain, to a
 * genuinely non-real number (Sqrt[x] at x<0).  Reduce's domain-restricted output
 * legitimately excludes such points, so the raw input must not be JUDGED there --
 * this is what keeps `1/x != 0 -> x != 0` and `Sqrt[x] != y -> x >= 0 && ...`
 * from spuriously "disagreeing" with their own inputs at the excluded points. *)
rdBadOperand[o_, pt_, dom_] := Module[{v, nv},
  v = o /. pt;
  If[! FreeQ[v, ComplexInfinity | Indeterminate | DirectedInfinity], Return[True]];
  If[! MemberQ[{Reals, Integers, Rationals}, dom], Return[False]];
  nv = N[v];
  NumericQ[nv] && TrueQ[Im[nv] != 0]];

(* Truth of the INPUT at pt, but -1 (skip) at any point where a relation operand
 * leaves the domain (see rdBadOperand). *)
rdInputTruth[e_, pt_, dom_] := Module[{ops},
  ops = Join[
    Cases[e, (Equal | Unequal | Less | Greater | LessEqual | GreaterEqual)[a_, b_]
            :> Sequence[a, b], {0, Infinity}],
    Cases[e, Inequality[a__] :> a, {0, Infinity}]];
  If[AnyTrue[ops, rdBadOperand[#, pt, dom] &], -1, rdTruth[e /. pt]]];

(* Build the full assignment lists (rule lists) for a Tuples grid over `syms`. *)
rdGridPoints[syms_, samp_] := Module[{n = Length[syms]},
  If[n == 0, {{}},
     Map[Function[tup, Thread[syms -> tup]], Tuples[samp, n]]]];

(* Witness points: for each conjunction-clause of `red` that pins every solve
 * variable to a value, instantiate the parameters and generated C[k] over a
 * small sample and keep the fully-numeric points. *)
rdWitness[red_, varsL_, psyms_, dom_] := Module[
  {clauses, pcombo, out = {}},
  clauses = If[Head[red] === Or, List @@ red, {red}];
  pcombo = If[Length[psyms] == 0, {{}},
     Map[Function[tup, Thread[psyms -> tup]],
         Tuples[rdSamples[dom, Length[psyms]], Length[psyms]]]];
  Scan[Function[cl,
     Module[{atoms},
       atoms = If[Head[cl] === And, List @@ cl, {cl}];
       Scan[Function[pa,
          Module[{inst, vrules, cand},
             inst   = atoms /. pa;
             vrules = Cases[inst, Equal[v_, val_] /; MemberQ[varsL, v] :> (v -> val)];
             cand   = Join[pa, vrules];
             (* Keep only points that (i) pin every solve variable to a fully
              * numeric value and (ii) genuinely lie in red -- the other atoms
              * of the clause (parametric conditions, side equations) must hold,
              * else this is not actually a solution red claims. *)
             If[Length[Union[First /@ vrules]] === Length[varsL]
                && AllTrue[Last /@ vrules, FreeQ[#, _Symbol] &]
                && rdTruth[red /. cand] === 1,
                AppendTo[out, cand]]]],
          pcombo]]],
     clauses];
  out];

(* Core equivalence engine -> {code, why}. *)
rdCheckPoints[exprE_, red_, varsL_, dom_] := Module[
  {params, genparams, syms, samp, gridPts, witPts, judged = 0, i, pt, lb, rb, why = ""},
  params    = Complement[rdFreeSyms[exprE], varsL];
  genparams = Union[Cases[red, Element[c_, _] :> c, Infinity]];
  syms      = Join[varsL, params];
  samp      = rdSamples[dom, Max[Length[syms], 1]];
  gridPts   = rdGridPoints[syms, samp];

  (* (a) grid: agree wherever both sides decide. *)
  Do[
     pt = gridPts[[i]];
     lb = rdInputTruth[exprE, pt, dom];
     rb = rdTruth[red   /. pt];
     If[lb =!= -1 && rb =!= -1,
        judged++;
        If[why === "" && lb =!= rb,
           why = "disagree at " <> ToString[pt, InputForm]
                 <> " (in=" <> ToString[lb] <> " out=" <> ToString[rb] <> ")"]],
     {i, Length[gridPts]}];
  If[why =!= "", Return[{1, why}]];

  (* (b) witnesses: each claimed solution must satisfy the input. *)
  witPts = rdWitness[red, varsL, Join[params, genparams], dom];
  Do[
     pt = witPts[[i]];
     lb = rdTruth[exprE /. pt];
     judged++;
     If[why === "" && lb =!= 1,
        why = "witness fails input at " <> ToString[pt, InputForm]
              <> " (in=" <> ToString[lb] <> ")"],
     {i, Length[witPts]}];
  If[why =!= "", Return[{1, why}]];

  If[judged == 0, Return[{1, "no judged points"}]];
  {0, ""}];

(* Core verdict: {code, why}. *)
SetAttributes[reduceVerdict, HoldAll];
reduceVerdict[expr_, vars_, dom_, expected_] := Module[
  {varsL, red, declined, exprN},
  varsL = If[Head[vars] === List, vars, {vars}];
  (* A List in the expr slot is the conjunction of its elements (Reduce accepts
   * it directly); normalise it for the back-sampler, whose truth oracle judges
   * a single relational/logical expression, not a raw List. *)
  exprN = If[Head[expr] === List, And @@ expr, expr];
  red   = If[dom === Automatic, Reduce[expr, vars], Reduce[expr, vars, dom]];
  declined = MatchQ[red, _Reduce] || ! FreeQ[red, _Reduce];

  If[expected === "decline",
     Return[If[declined, {0, ""},
               {1, "expected decline, got: " <> ToString[red, InputForm]}]]];

  If[declined,
     Return[{2, "unevaluated: " <> ToString[red, InputForm]}]];

  If[expected === True && red =!= True,
     Return[{1, "expected True, got: " <> ToString[red, InputForm]}]];
  If[expected === False && red =!= False,
     Return[{1, "expected False, got: " <> ToString[red, InputForm]}]];
  If[! MemberQ[{True, False, "solved"}, expected],
     Return[{1, "bad expected tag: " <> ToString[expected]}]];

  rdCheckPoints[exprN, red, varsL, dom]];

(* Integer-returning entry used by the C runner; prints a diagnostic on any
 * non-PASS so the CTest log names the failing case. *)
SetAttributes[reduceCheckCode, HoldAll];
reduceCheckCode[label_, expr_, vars_, dom_, expected_] := Module[{v},
  v = reduceVerdict[expr, vars, dom, expected];
  If[v[[1]] =!= 0,
     Print["REDUCE-", If[v[[1]] === 2, "UNEVAL", "FAIL"], "\t", label, "\t", v[[2]]]];
  v[[1]]];

(* TSV-printing entry used by a standalone driver. *)
SetAttributes[reduceReport, HoldAll];
reduceReport[label_, expr_, vars_, dom_, expected_] := Module[{v, tag},
  v = reduceVerdict[expr, vars, dom, expected];
  tag = Switch[v[[1]], 0, "PASS", 2, "UNEVAL", _, "FAIL"];
  Print["REDUCE\t", tag, "\t", label, "\t", v[[2]]];
  v[[1]]];
