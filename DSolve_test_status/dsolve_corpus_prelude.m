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
 *              branch's residual is proved zero by PossibleZeroQ -- for an IVP
 *              the ODE residual AND every initial condition (implicit /
 *              parametric branches are accepted on DSolve's own symbolic
 *              self-verification).
 *   1 = FAIL   DSolve returned a closed form and PossibleZeroQ says an explicit
 *              branch's residual is NOT zero (a wrong answer).
 *   2 = UNEVAL DSolve declined (bubbled back unevaluated, $Aborted, or {}), OR
 *              returned an IVP solution that still carries a generated constant
 *              C[k] and fails its conditions -- the general solution with the
 *              initial condition unfitted, so the IVP is not actually solved.
 *   3 = SKIP   Reserved for a record that is not a solvable ODE/system (e.g. a
 *              converter placeholder).  Systems are NO LONGER skipped -- they are
 *              solved and back-substituted like any other record (M27).
 *
 * Verification is by BACK-SUBSTITUTION, never string comparison, so it is
 * immune to solution spelling, and PossibleZeroQ is THE verifier: the residual
 * keeps C[k] and every parameter symbolic, there is no numeric acceptance path,
 * and a PossibleZeroQ false negative is fixed in PossibleZeroQ.  (Until
 * v0.6161 a branch passed when a majority of six real sample points were small,
 * which accepted answers valid only on part of the line.)  See dsResidVerdict
 * for the branch policy and dsMemberResid for conditions at a removable
 * singularity.
 *)

(* The harness runs core_init() only, never src/internal/init.m, so without
 * this the corpus measured DSolve WITHOUT the distribution rules every REPL
 * session loads: DiracDelta[-1] and HeavisideTheta[-1] stayed unevaluated in
 * impulse-response answers (2.2.40-3976 ...), which the old numeric verifier
 * could not numericize and so accepted as "UNK".  init.m loads exactly this
 * module; keep the two in step. *)
LoadModule["distributions.m"];

$dsSolveTimeout = 8;          (* seconds per DSolve call *)
$dsProtected    = {E, Pi, I, Infinity, ComplexInfinity, Indeterminate,
                   True, False, Degree, EulerGamma, GoldenRatio, Catalan,
                   Complex, Rational, Integer, Real, List,
                   (* Relational / logical / structural operators.  These are
                    * never ODE parameters (the converter emits value symbols
                    * only), but a CHAINED inequality in a Piecewise/step
                    * condition -- `0 <= t < 2` -> Inequality[0, LessEqual, t,
                    * Less, 2] -- puts the operators in ARGUMENT positions, where
                    * Cases[..., s_Symbol] below collects them as bogus "free
                    * parameters".  Substituting numbers for LessEqual/Less then
                    * corrupts the Piecewise into garbage that SIGSEGVs when the
                    * residual is differentiated (2.2.27-2690).  Single-sided
                    * conditions (`Less[t, 2]`, operator as head) never hit this,
                    * which is why the piecewise-heavy S2.2.15 stayed green. *)
                   Less, LessEqual, Greater, GreaterEqual, Equal, Unequal,
                   Inequality, And, Or, Not, Xor, Piecewise, UnitStep,
                   HeavisideTheta, DiracDelta};

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

(* Verdict for a single residual expression: "OK" (PossibleZeroQ proves it
 * zero), "BAD" (PossibleZeroQ says it is not zero), or "UNK" (PossibleZeroQ did
 * not answer within $dsPZQTimeout).
 *
 * PossibleZeroQ is THE verifier.  The residual keeps the generated constants C[k]
 * and every ODE parameter SYMBOLIC, so a PASS means the branch satisfies the ODE
 * identically -- for every value of the constants and for all complex iv -- not
 * merely at a handful of real sample points.  The verifier this replaces scored
 * a branch by a majority of six real points (with 20/210-digit retries), which
 * accepted answers valid only on part of the line: 2.2.16-1582,
 * x^2 y y' == (y^2 - 1)^(3/2), whose explicit inversion holds only where
 * Re(x/(1 - C x)) > 0 because Sqrt[w^2] == w is not an identity.  There is no
 * numeric acceptance path: a PossibleZeroQ false negative is a PossibleZeroQ bug,
 * to be fixed there, never rescued here.  Numerics are for triage only.
 *
 * Messages are muted: the residual is a probe, and a leaked diagnostic from it
 * would be attributed to the record. *)
(* Branch policy.  An ODE whose right-hand side carries a branch -- y' = x
 * Sqrt[y], y' == (y^2 - 1)^(3/2) -- has NO explicit solution valid for all
 * complex x: y = (x^2 + 2 C)^2/16 satisfies it only where x^2 + 2 C is on the
 * principal sheet, and the reference system's answers have the same property.
 * Such a branch is accepted when PossibleZeroQ proves the PowerExpand-ed
 * residual zero: PowerExpand takes every power on its principal branch
 * consistently, so this certifies "valid on the principal region" -- the local
 * notion of an ODE solution -- symbolically, with no numeric acceptance.  A
 * residual that is not zero even then is a wrong answer (BAD). *)
$dsPZQTimeout = 4;
(* Every radicand factored (exactly), so PowerExpand can see a square that was
 * expanded inside the radical: Sqrt[C^2 + 2 C Log[x] + Log[x]^2]. *)
dsFactorRadicands[e_] := Quiet@TimeConstrained[
  e /. Power[b_, p_Rational] /; !NumericQ[b] :> Power[Factor[b], p], $dsPZQTimeout, e];
dsPZQ[e_] := Quiet@TimeConstrained[PossibleZeroQ[e], $dsPZQTimeout, $Aborted];
dsResidVerdict[resid_, iv_] := Module[{z, zp},
  z = dsPZQ[resid];
  If[z === True, Return["OK"]];
  (* undecided within the budget: the PowerExpand variants are no cheaper, and
   * three more budgets per member would outrun the harness's per-case alarm *)
  If[z === $Aborted, Return["UNK"]];
  zp = dsPZQ[PowerExpand[resid]];
  If[zp =!= True, zp = dsPZQ[PowerExpand[dsFactorRadicands[resid]]]];
  Which[zp === True, "OK", z === False || zp === False, "BAD", True, "UNK"]
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

(* The residual of one member of the equation slot under branch br.  A point
 * condition y[x0] == v / y'[x0] == v whose direct substitution is Indeterminate
 * or infinite sits on a removable singularity of the closed form (e.g.
 * (x^2/2)^(1/4)/Sqrt[x] at x0 = 0): the condition means the two-sided LIMIT
 * there, so y^(k)[x0] is replaced by Limit[D[y[iv] /. br, {iv, k}], iv -> x0].
 * A branch-restricted answer has no two-sided limit and still fails. *)
dsMemberResid[eq_, br_, fn_, iv_] := Module[{fns = dsFnList[fn], r},
  r = Quiet[(eq /. Equal -> Subtract) /. br];   (* a probe: 1/0 at a removable point *)
  If[!FreeQ[eq, iv] || FreeQ[r, Indeterminate | ComplexInfinity | DirectedInfinity],
    Return[r]];
  Quiet@TimeConstrained[
    (eq /. Equal -> Subtract) /. {
      Derivative[k_][f_Symbol][x0_] /; MemberQ[fns, f] && FreeQ[x0, iv] :>
        Limit[D[f[iv] /. br, {iv, k}], iv -> x0],
      f_Symbol[x0_] /; MemberQ[fns, f] && FreeQ[x0, iv] :> Limit[f[iv] /. br, iv -> x0]},
    $dsPZQTimeout, r]
];

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
  resids = dsMemberResid[#, br, fn, iv] & /@ eqs;
  verds  = dsResidVerdict[#, iv] & /@ resids;
  (* A failing branch that still carries C[k] is UNFIT -- the general solution
   * returned without fitting the conditions -- only for an IVP/BVP.  For a plain
   * ODE every general solution carries C[k], so this would score every wrong
   * answer as UNEVAL instead of FAIL. *)
  If[MemberQ[verds, "BAD"],
    Return[If[leaked && dsIsIVP[eqn, iv], "UNFIT", "BAD"]]];
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
