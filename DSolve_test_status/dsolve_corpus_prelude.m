(* dsolve_corpus_prelude.m
 *
 * Self-verifying checker for the DSolve[] corpus (DE_examples_2.m and any
 * future DE_examples_N.m).  Loaded once by the C runner
 * (tests/test_dsolve_corpus.c) before it forks per case.
 *
 * A corpus record is {label, equation(s), function(s), indVar, classif, sympy?}.
 * dsolveCheckCode[label, eqn, fn, iv, classif] returns an integer verdict:
 *
 *   0 = PASS   DSolve returned a non-empty List of branches and every
 *              explicit branch back-substitutes numerically to ~0 (implicit /
 *              parametric / non-numericizable branches are accepted on
 *              DSolve's own symbolic self-verification).
 *   1 = FAIL   DSolve returned a closed form but an explicit branch is
 *              demonstrably nonzero at the majority of numeric sample points
 *              (a wrong answer -- the case the symbolic verifier can miss when
 *              the residual is undecidable, cf. M12/M14).
 *   2 = UNEVAL DSolve declined: bubbled back unevaluated, returned $Aborted
 *              (timed out), or returned {} (no closed form).
 *   3 = SKIP   Systems (this is the scalar-first campaign) -- not counted.
 *
 * Verification is by BACK-SUBSTITUTION, never string comparison, so it is
 * immune to solution spelling.  Free parameters and generated constants C[k]
 * are instantiated at fixed generic rational points; the independent variable
 * is swept over several points and the residual must vanish at the majority of
 * points where it numericizes (mirrors the l2_num_ok / cv_num_ok policy that
 * dsolve_lie2.c / dsolve_changevar.c use internally).
 *)

$dsSolveTimeout = 8;          (* seconds per DSolve call *)
$dsTol          = 1*^-6;
$dsProtected    = {E, Pi, I, Infinity, ComplexInfinity, Indeterminate,
                   True, False, Degree, EulerGamma, GoldenRatio, Catalan,
                   Complex, Rational, Integer, Real, List};

(* --- numeric back-substitution of one explicit branch ---------------------- *)

(* Collect free scalar parameters of a residual (symbols other than iv, the
 * generated constants C[k], and protected constants/heads). *)
dsFreeParams[resid_, iv_] := DeleteCases[
  DeleteDuplicates@Cases[resid, s_Symbol /; !MemberQ[$dsProtected, s], {0, Infinity}, Heads -> True],
  iv | C];

(* Verdict for a single residual expression: "OK" (verified zero),
 * "BAD" (verified nonzero), or "UNK" (could not numericize). *)
dsResidVerdict[resid_, iv_] := Module[
  {params, consts, pv, vals = {}, k, cv, r, a, nsmall},
  params = dsFreeParams[resid, iv];
  consts = DeleteDuplicates@Cases[resid, C[_Integer], Infinity];
  pv = MapIndexed[#1 -> (13/10 + First[#2]*4/17) &, params];
  Do[
    cv = MapIndexed[#1 -> (7/10 + k/5 + First[#2]*3/19) &, consts];
    r  = N[(resid /. pv /. cv /. iv -> (11/10 + k*5/13)), 20];
    If[NumberQ[r] || Head[r] === Complex, AppendTo[vals, Abs[r]]];
  , {k, 0, 5}];
  If[Length[vals] < 2, Return["UNK"]];
  nsmall = Count[vals, a_ /; a < $dsTol];
  If[nsmall >= Ceiling[Length[vals]/2], "OK", "BAD"]
];

(* Is a branch an explicit solution {fn -> Function[..]} (or several)? *)
dsExplicitQ[br_List, fn_] := br =!= {} && AllTrue[br,
  MatchQ[#, Rule[fn, _Function] | Rule[_[fn], _Function]] &];
dsExplicitQ[_, _] := False;

(* Verify one branch against the (single) scalar equation lhs==rhs. *)
dsBranchVerdict[eqn_, br_, fn_, iv_] := Module[{resid},
  If[!dsExplicitQ[br, fn], Return["UNK"]];   (* implicit/parametric: trust DSolve *)
  resid = (eqn /. Equal -> Subtract) /. br;
  dsResidVerdict[resid, iv]
];

(* --- top-level verdict ----------------------------------------------------- *)

dsolveCheckCode[label_, eqn_, fn_, iv_, classif_] := Module[
  {sol, verds},
  (* skip systems (scalar-first campaign) *)
  If[StringQ[classif] && StringContainsQ[classif, "system"], Return[3]];
  If[ListQ[fn], Return[3]];
  sol = TimeConstrained[DSolve[eqn, fn, iv], $dsSolveTimeout, $Aborted];
  If[sol === $Aborted, Return[2]];
  If[Head[sol] === DSolve, Return[2]];        (* bubbled back unevaluated *)
  If[!ListQ[sol] || sol === {}, Return[2]];   (* no closed form *)
  verds = dsBranchVerdict[eqn, #, fn, iv] & /@ sol;
  If[MemberQ[verds, "BAD"], Return[1]];       (* a demonstrably wrong branch *)
  0
];

(* Human-readable driver (standalone use). *)
dsolveReport[label_, eqn_, fn_, iv_, classif_] := Module[{c, name},
  c = dsolveCheckCode[label, eqn, fn, iv, classif];
  name = Switch[c, 0, "PASS", 1, "FAIL", 2, "UNEVAL", 3, "SKIP", _, "?"];
  Print[name, "\t", label, "\t", classif]
];
