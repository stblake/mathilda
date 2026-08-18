(* solve_check_prelude.m
 *
 * Form-invariant back-substitution verifier for the Solve[] stress
 * corpus (tests/solve_corpus.m).  Loaded once by the C runner
 * (tests/test_solve_corpus.c) before it forks per case, and by the
 * standalone driver (tests/solve_corpus_run.m).
 *
 * The contract a case must satisfy, in order:
 *   1. Solve[eqn, vars, dom] must NOT bubble back unevaluated.
 *   2. The number of returned solution rule-lists must equal `expected`.
 *   3. Every returned rule-list, substituted into every equation of the
 *      system, must drive the residual lhs-rhs to zero.
 *   4. Concrete solution values must lie in the requested domain.
 *
 * Verification is by BACK-SUBSTITUTION, never by string comparison, so
 * the corpus is immune to solution ordering and radical/formatting
 * spelling.  The residual-is-zero decision mirrors solverad.c's
 * verify_candidate ladder: PossibleZeroQ -> Simplify -> |N[.]| < tol,
 * with Root[] objects accepted by numericization only and periodic
 * ConditionalExpression[value, Element[C[k],Integers]] families checked
 * across a small integer sample.
 *
 * solveVerdict returns {code, why} with code 0 = PASS, 1 = FAIL,
 * 2 = UNEVALUATED.  solveCheckCode returns just the integer (used by the
 * C runner); solveReport prints a TSV line (used by the standalone driver).
 *)

$solveTol = 1*^-9;
$ckSamples = {-2, -1, 0, 1, 2};

(* Numeric zero test: does the residual numericize to ~0? *)
scNumZero[r_] := Module[{n},
  n = N[r];
  NumberQ[n] && Abs[n] < $solveTol
];

(* Residual-is-zero for a plain algebraic / radical / complex residual
 * (no Root and no ConditionalExpression present). *)
scZeroPlain[r_] := Module[{s},
  If[TrueQ[PossibleZeroQ[r]], Return[True]];
  s = Simplify[r];
  If[s === 0 || TrueQ[PossibleZeroQ[s]], Return[True]];
  scNumZero[s]
];

(* Residual-is-zero, shape-aware. *)
scZero[r_] := Module[{params, stripped},
  (* Root[] objects describe the unique root of an irreducible factor;
   * accept via numericization only, never back-substitute. *)
  If[!FreeQ[r, Root], Return[scNumZero[r]]];
  (* Periodic / branch-cut ConditionalExpression: pull the generated
   * integer parameters out of the Element[_,Integers] predicates, strip
   * the wrapper to its value, and require the residual to vanish for
   * every integer sample. *)
  If[!FreeQ[r, ConditionalExpression],
     params  = Union[Cases[r, Element[q_, Integers] :> q, Infinity]];
     stripped = r /. ConditionalExpression[v_, _] :> v;
     If[params === {}, Return[scZero[stripped]]];
     (* NB: Mathilda's named Function does not close over its variable
      * inside a nested pure function, so build the sample substitution
      * with Thread (params -> m) rather than a nested Map. *)
     Return[AllTrue[$ckSamples,
        Function[m, scZero[stripped /. Thread[params -> m]]]]]
  ];
  scZeroPlain[r]
];

(* Verify one solution rule-list `rl` against one equation lhs==rhs. *)
scVerifyOne[lhs_ == rhs_, rl_, dom_] := Module[{resid},
  resid = (lhs - rhs) /. rl;
  If[MatchQ[dom, Modulus -> _],
     Return[TrueQ[Mod[resid, dom[[2]]] == 0] || scNumZero[resid]]];
  scZero[resid]
];
(* Non-equation operands (a system element that evaluated to True) are
 * trivially satisfied. *)
scVerifyOne[True, _, _] := True;
scVerifyOne[_, _, _] := True;

(* Does a concrete solution value lie in the requested domain? *)
scInDomain[val_, dom_] := Which[
  dom === Reals,             TrueQ[PossibleZeroQ[Im[val]]] || scNumZero[Im[val]],
  dom === Integers,          IntegerQ[val],
  dom === Rationals,         MatchQ[val, _Integer | _Rational],
  MatchQ[dom, Modulus -> _], IntegerQ[val] && 0 <= val < dom[[2]],
  True,                      True   (* Complexes / no domain: unconstrained *)
];

(* Core verdict: {code, why}. *)
SetAttributes[solveVerdict, HoldAll];
solveVerdict[eqn_, vars_, dom_, expected_] := Module[
  {sol, eqs, i, rl, ok = True, why = ""},

  sol = If[dom === Automatic, Solve[eqn, vars], Solve[eqn, vars, dom]];

  (* (1) unevaluated bubble-back. *)
  If[MatchQ[sol, _Solve] || !FreeQ[sol, _Solve],
     Return[{2, "unevaluated: " <> ToString[sol, InputForm]}]];

  (* (2) solution count. *)
  If[Length[sol] =!= expected,
     Return[{1, "count " <> ToString[Length[sol]] <> " != "
                 <> ToString[expected] <> "  " <> ToString[sol, InputForm]}]];

  (* Normalise the equation(s) to a flat list of _Equal operands. *)
  eqs = Which[Head[eqn] === List, eqn,
              Head[eqn] === And,  List @@ eqn,
              True,               {eqn}];
  eqs = Select[eqs, MatchQ[#, _Equal] &];

  (* (3) back-substitution + (4) domain membership, per solution. *)
  Do[
     rl = sol[[i]];
     Scan[Function[eq,
        If[!TrueQ[scVerifyOne[eq, rl, dom]],
           ok = False;
           why = "resid!=0 sol#" <> ToString[i] <> "  " <> ToString[sol, InputForm]]],
        eqs];
     Scan[Function[rule,
        Module[{v = rule[[2]]},
           If[FreeQ[v, ConditionalExpression] && FreeQ[v, Root] && FreeQ[v, C]
              && !TrueQ[scInDomain[v, dom]],
              ok = False;
              why = "domain sol#" <> ToString[i] <> "  " <> ToString[sol, InputForm]]]],
        rl],
     {i, Length[sol]}];

  {If[ok, 0, 1], why}
];

(* Integer-returning entry used by the C runner; prints a diagnostic on
 * any non-PASS so the CTest log names the failing case. *)
SetAttributes[solveCheckCode, HoldAll];
solveCheckCode[label_, eqn_, vars_, dom_, expected_] := Module[{v},
  v = solveVerdict[eqn, vars, dom, expected];
  If[v[[1]] =!= 0,
     Print["SOLVE-", If[v[[1]] === 2, "UNEVAL", "FAIL"], "\t", label, "\t", v[[2]]]];
  v[[1]]
];

(* TSV-printing entry used by the standalone driver. *)
SetAttributes[solveReport, HoldAll];
solveReport[label_, eqn_, vars_, dom_, expected_] := Module[{v, tag},
  v = solveVerdict[eqn, vars, dom, expected];
  tag = Switch[v[[1]], 0, "PASS", 2, "UNEVAL", _, "FAIL"];
  Print["SOLVE\t", tag, "\t", label, "\t", v[[2]]];
  v[[1]]
];
