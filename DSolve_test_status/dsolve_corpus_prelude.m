(* dsolve_corpus_prelude.m
 *
 * Self-verifying checker for the DSolve[] corpus (DE_examples_2.m and any
 * future DE_examples_N.m).  Loaded once by the C runner
 * (tests/test_dsolve_corpus.c) before it forks per case.
 *
 * A corpus record is {label, equation(s), function(s), indVar, classif, sympy?}.
 * The equation slot is a single ODE for a general-solution problem, or the
 * DSolve-native list {ode, ic1, ...} for an initial-value problem (each ic a
 * point equation y[x0]==v / y'[x0]==v).  For a SYSTEM the equation slot is a
 * List of equations {ode1, ode2, ...(, ic1, ...)} and the function slot is a
 * List of the dependent functions {x, y, ...}; a system is verified exactly like
 * a scalar ODE -- every branch is a rule-list {x->Function[..], y->Function[..]}
 * back-substituted into every equation (see dsExplicitQ / dsBranchVerdict).
 * dsolveCheckCode[label, eqn, fn, iv, classif] returns an integer verdict:
 *
 *   0 = PASS   DSolve returned a non-empty List of branches and every explicit
 *              branch back-substitutes numerically to ~0 -- for an IVP that is
 *              the ODE residual AND every initial condition (implicit /
 *              parametric / non-numericizable branches are accepted on DSolve's
 *              own symbolic self-verification).
 *   1 = FAIL   DSolve returned a closed form but an explicit branch is
 *              demonstrably nonzero at the majority of numeric sample points
 *              (a wrong answer -- the case the symbolic verifier can miss when
 *              the residual is undecidable, cf. M12/M14).
 *   2 = UNEVAL DSolve declined (bubbled back unevaluated, $Aborted, or {}), OR
 *              returned an IVP solution that still carries a generated constant
 *              C[k] -- i.e. the general solution with the initial condition
 *              unfitted, so the IVP is not actually solved (not a wrong answer).
 *   3 = SKIP   Reserved for a record that is not a solvable ODE/system (e.g. a
 *              converter placeholder).  Systems are NO LONGER skipped -- they are
 *              solved and back-substituted like any other record (M27).
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

(* Collect free scalar parameters of a residual: ARGUMENT-POSITION symbols only
 * (other than iv, the generated constants C[k], and protected constants).  NOT
 * Heads -> True -- a free parameter is a value, never an operator/function head,
 * and collecting heads (Plus, Times, Tan, Sec, ...) then substituting numbers for
 * them turns every residual into non-numericizable garbage -> a vacuous "UNK".
 * (Matches the internal l2_num_ok/cv_num_ok policy: instantiate params, never heads.) *)
dsFreeParams[resid_, iv_] := DeleteCases[
  DeleteDuplicates@Cases[resid, s_Symbol /; !MemberQ[$dsProtected, s], {0, Infinity}],
  iv | C];

(* Verdict for a single residual expression: "OK" (verified zero),
 * "BAD" (verified nonzero), or "UNK" (could not numericize). *)
(* The sweep index and the residual-value holder are $-prefixed so they can never
 * collide with an ODE PARAMETER of the same name: the loop variable dynamically
 * rebinds every occurrence of its symbol while the residual is evaluated, so a
 * plain `k` loop over a residual that carries a parameter `k` (e.g. §2.2.4-387,
 * `m x''+k x==F0 Cos[om t]`) would force the spring constant to 0..5 instead of
 * its generic sample value -> a bogus nonzero residual -> a FALSE "BAD". The
 * converter never emits a `$`-prefixed symbol, so `$dsSweep`/`$dsVal` are safe. *)
dsResidVerdict[resid_, iv_] := Module[
  {params, consts, pv, vals = {}, $dsSweep, cv, $dsVal, $dsRex, a, nsmall},
  params = dsFreeParams[resid, iv];
  consts = DeleteDuplicates@Cases[resid, C[_Integer], Infinity];
  pv = MapIndexed[#1 -> (13/10 + First[#2]*4/17) &, params];
  Do[
    cv = MapIndexed[#1 -> (7/10 + $dsSweep/5 + First[#2]*3/19) &, consts];
    $dsRex = (resid /. pv /. cv /. iv -> (11/10 + $dsSweep*5/13));  (* exact *)
    $dsVal = N[$dsRex, 20];
    (* A fast-growing (large-eigenvalue) solution back-substitutes to a residual that
     * is a difference of huge E^{lambda x} terms, so a TRUE zero looks large at
     * 20-digit precision -- catastrophic cancellation (e.g. the eigenvalue-64 system
     * 2.2.11-1001: |resid|@20 = 8.9*^43 at x=3, but @120 = 2*^-56).  When the 20-digit
     * value is NOT already small, re-evaluate the EXACT residual at high precision:
     * a genuine nonzero stays nonzero, a cancellation artifact collapses to ~0.  This
     * only ever turns a spurious "not small" into "small" -- monotone: it can lower a
     * section's non-PASS count, never raise it, and never introduces a FAIL. *)
    If[(NumberQ[$dsVal] || Head[$dsVal] === Complex) && Abs[$dsVal] >= $dsTol,
       $dsVal = N[$dsRex, 200]];
    If[NumberQ[$dsVal] || Head[$dsVal] === Complex, AppendTo[vals, Abs[$dsVal]]];
  , {$dsSweep, 0, 5}];
  If[Length[vals] < 2, Return["UNK"]];
  nsmall = Count[vals, a_ /; a < $dsTol];
  If[nsmall >= Ceiling[Length[vals]/2], "OK", "BAD"]
];

(* The target function symbol of an explicit rule (fn -> Function[..] or the
 * applied form fn[iv] -> Function[..]); $Failed for an implicit/parametric rule. *)
dsRuleFn[Rule[f_Symbol, _Function]]      := f;
dsRuleFn[Rule[(_)[f_Symbol], _Function]] := f;
dsRuleFn[_]                              := $Failed;

(* Normalise the function slot to a list of dependent-function symbols: a scalar
 * problem passes the bare symbol y, a system passes the List {x, y, ...}. *)
dsFnList[fn_List] := fn;
dsFnList[fn_]     := {fn};

(* Is a branch a fully explicit solution?  Every rule must resolve to a Function
 * (no implicit/parametric member) AND every dependent function must be present
 * -- so a single scalar {y->Function[..]} and a system {x->Function[..],
 * y->Function[..]} are both explicit, but a partially-solved system is not. *)
dsExplicitQ[br_List, fn_] := Module[{fns = dsFnList[fn], got},
  got = dsRuleFn /@ br;
  br =!= {} && FreeQ[got, $Failed] && Complement[fns, got] === {}];
dsExplicitQ[_, _] := False;

(* An IVP/BVP problem gives the equation slot as a List {ode, ic1, ...}; a point
 * condition (y[x0]==v) is free of the independent variable iv, the ODE is not. *)
dsIsIVP[eqn_, iv_] := ListQ[eqn] && AnyTrue[eqn, FreeQ[#, iv] &];

(* Verify one branch against the equation slot, which is either a single scalar
 * ODE lhs==rhs, or (for an IVP) a List {ode, ic1, ...}: every member -- the ODE
 * residual AND each initial condition -- must back-substitute to ~0.  Verdicts:
 *   "OK"    every member numerically verified (>=1 member decided, none BAD)
 *   "BAD"   a member is demonstrably nonzero and no generated constant remains
 *           (a genuinely wrong answer, cf. M12/M14)
 *   "UNFIT" a member fails BUT the branch still carries a generated constant
 *           C[k] -- DSolve returned the general solution without fitting the
 *           initial condition, so the IVP is UNSOLVED (scored UNEVAL, not a
 *           wrong answer)
 *   "UNK"   implicit/parametric or non-numericizable: trust DSolve's own verify *)
dsBranchVerdict[eqn_, br_, fn_, iv_] := Module[{eqs, resids, verds, leaked},
  eqs = If[ListQ[eqn], eqn, {eqn}];
  leaked = !FreeQ[br, C[_Integer]];
  If[!dsExplicitQ[br, fn],                     (* implicit/parametric *)
    Return[If[dsIsIVP[eqn, iv] && leaked, "UNFIT", "UNK"]]];
  resids = ((# /. Equal -> Subtract) /. br) & /@ eqs;
  verds  = dsResidVerdict[#, iv] & /@ resids;
  If[MemberQ[verds, "BAD"], Return[If[leaked, "UNFIT", "BAD"]]];
  If[AllTrue[verds, # === "UNK" &], "UNK", "OK"]
];

(* --- top-level verdict ----------------------------------------------------- *)

dsolveCheckCode[label_, eqn_, fn_, iv_, classif_] := Module[
  {sol, verds},
  (* Systems (fn a List) are solved and verified exactly like scalars (M27): a
   * branch is a rule-list {x->Function[..], ...} back-substituted per equation. *)
  sol = TimeConstrained[DSolve[eqn, fn, iv], $dsSolveTimeout, $Aborted];
  If[sol === $Aborted, Return[2]];
  If[Head[sol] === DSolve, Return[2]];        (* bubbled back unevaluated *)
  If[!ListQ[sol] || sol === {}, Return[2]];   (* no closed form *)
  verds = dsBranchVerdict[eqn, #, fn, iv] & /@ sol;
  If[MemberQ[verds, "BAD"], Return[1]];       (* a demonstrably wrong branch *)
  If[MemberQ[verds, "UNFIT"], Return[2]];     (* IVP constant unfitted -> unsolved *)
  0
];

(* Human-readable driver (standalone use). *)
dsolveReport[label_, eqn_, fn_, iv_, classif_] := Module[{c, name},
  c = dsolveCheckCode[label, eqn, fn, iv, classif];
  name = Switch[c, 0, "PASS", 1, "FAIL", 2, "UNEVAL", 3, "SKIP", _, "?"];
  Print[name, "\t", label, "\t", classif]
];
