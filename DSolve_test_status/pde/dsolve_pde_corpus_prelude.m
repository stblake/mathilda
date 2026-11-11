(* dsolve_pde_corpus_prelude.m
 *
 * Self-verifying checker for the DSolve PDE corpus (DE_PDE_*.m, imported from
 * N. Abbasi's "PDE in CAS" by tools/pde_corpus/import_abbasi_pde.py).  Loaded
 * once by the C runner (test_dsolve_pde_corpus.c) before it forks per case.
 *
 * A record is {label, section, eqns, fn, vars, opts, colA, colB}: eqns is the
 * PDE or a List {pde, conditions...} (conditions may themselves be Lists, so
 * the slot is Flatten-ed), fn the unknown as given (u or u[x, t]), vars the
 * independent variables, opts extra DSolve options (Assumptions -> ...).
 *
 * pdeCheckCode[label, section, eqns, fn, vars, opts] returns
 *
 *   0 PASS          every explicit branch back-substitutes into the PDE and every
 *                   condition: PossibleZeroQ proves each residual zero, with the
 *                   arbitrary functions C[k][...] symbolic, or -- when symbolic
 *                   PossibleZeroQ does not decide -- for BOTH of two generic
 *                   pinnings of every C[k] (a residual that vanishes for two
 *                   unrelated nonlinear functions of every arity).
 *   1 FAIL          a residual is proved NONZERO: a wrong answer.
 *   2 UNEVAL        DSolve declined: unevaluated, {}, or the 8 s budget ran out.
 *   3 SKIP          reserved (a record that is not a PDE problem).
 *   4 PASS-IMPLICIT an implicit answer phi == C[1][psi] (first order), verified
 *                   by the implicit-function rule: u_v = -R_v/R_u substituted into
 *                   the PDE vanishes identically in (vars, u).
 *   5 PASS-SERIES   an eigenfunction series Inactive[Sum][term, {K[1], lo, hi}]:
 *                   the generic mode (K[1] a symbolic integer) satisfies the
 *                   homogeneous PDE and every homogeneous condition, the non-sum
 *                   part carries the inhomogeneity, the modes are orthogonal on
 *                   the interval the conditions bound, and the coefficient is the
 *                   orthogonal projection of the expanded data.
 *   6 UNK           an answer the checks can neither prove nor refute (a
 *                   distributional condition such as the heat kernel at t = 0, a
 *                   higher-order implicit relation, a Bessel series).
 *   7 UNFIT         an initial/boundary-value answer that still carries C[k] and
 *                   fails a condition: the general solution, not the IVP's.
 *
 * PossibleZeroQ is THE verifier (Simplify under the K[1]-integer assumption only
 * where a series needs Sin[K[1] Pi] == 0); there is no numeric acceptance path.
 * A PossibleZeroQ false negative is a PossibleZeroQ bug, fixed there.
 *)

LoadModule["distributions.m"];

$pdSolveTimeout = 8;    (* seconds per DSolve call *)
$pdPZQTimeout   = 4;    (* seconds per PossibleZeroQ / Integrate probe *)

(* first u-occurrence argument list in e, or $Failed.  (Cases with a count,
 * not FirstCase: the 4-argument FirstCase[e, p, default, levelspec] form is not
 * implemented yet -- tasks/todo.md.) *)
pdFirstArgs[e_, u_] := Replace[
  Cases[e, (u | Derivative[__][u])[a__] :> {a}, {0, Infinity}, 1], {{a_} :> a, _ -> $Failed}];

pdTC[e_, d_] := Quiet@TimeConstrained[e, $pdPZQTimeout, d];
SetAttributes[pdTC, HoldFirst];

(* --- one residual --------------------------------------------------------- *)

pdFactorRadicands[e_] := pdTC[
  e /. Power[b_, p_Rational] /; !NumericQ[b] :> Power[Factor[b], p], e];
pdPZQ[e_] := pdTC[PossibleZeroQ[e], $Aborted];

(* "OK" / "BAD" / "UNK" for a residual; the PowerExpand rungs certify the
 * principal-branch region, exactly as the ODE corpus prelude does. *)
pdResidVerdict[resid_] := Module[{z, zp},
  If[!FreeQ[resid, Indeterminate | ComplexInfinity | DirectedInfinity], Return["UNK"]];
  z = pdPZQ[resid];
  If[z === True, Return["OK"]];
  If[z === $Aborted, Return["UNK"]];
  zp = pdPZQ[PowerExpand[resid]];
  If[zp =!= True, zp = pdPZQ[PowerExpand[pdFactorRadicands[resid]]]];
  Which[zp === True, "OK", z === False || zp === False, "BAD", True, "UNK"]
];

(* --- generic pinnings of the arbitrary functions C[k][a1, ..., am] --------- *)

pdPinA[k_, a_List] := Sin[k + Range[Length[a]] . a] + Total[a]^2 + k;
pdPinB[k_, a_List] := Exp[Total[a]/(k + 1)] + Times @@ a + 1;
pdPin[e_, Identity] := e;               (* symbolic: C[k] left as is *)
pdPin[e_, pin_] := e /. C[k_Integer][a___] :> pin[k, {a}];
pdHasArbFun[e_] := !FreeQ[e, C[_Integer][___]];

(* Verdict of a residual MAKER res[pin] (pin = Identity for symbolic C[k]). *)
pdPinnedVerdict[res_] := Module[{v, va, vb},
  v = pdResidVerdict[res[Identity]];
  If[v === "OK" || !pdHasArbFun[res[Identity]], Return[v]];
  va = pdResidVerdict[res[pdPinA]];
  vb = pdResidVerdict[res[pdPinB]];
  Which[va === "OK" && vb === "OK", "OK",
        v === "BAD" || va === "BAD" || vb === "BAD", "BAD",
        True, "UNK"]
];

(* --- problem structure ----------------------------------------------------- *)

pdUSym[fn_Symbol] := fn;
pdUSym[fn_] := Head[fn];

(* A member is a CONDITION when no occurrence of u (or a derivative of u) is
 * evaluated at the full variable list -- u[x, 0], Derivative[1, 0][u][0, t]. *)
pdIsCondition[eq_, u_, vars_] :=
  FreeQ[eq, u @@ vars] && FreeQ[eq, Derivative[__][u] @@ vars];

(* An explicit rule -> a Function of the variables; $Failed otherwise. *)
pdRuleFun[Rule[_Symbol, f_Function], _] := f;
pdRuleFun[Rule[_[a__], body_], _] := Function @@ {{a}, body};
pdRuleFun[_, _] := $Failed;

(* Residual of member eq with u replaced by the Function F. *)
pdSubst[eq_, u_, F_] := Quiet[(eq /. Equal -> Subtract) /. u -> F];

(* A condition whose direct substitution is singular (a removable singularity
 * of the closed form on the data line): take the one-sided limit onto the line
 * of each fixed slot, as the ODE prelude does at a point. *)
pdCondResid[eq_, u_, F_, vars_] := Module[{r, occ, a, fixed, rr},
  r = pdSubst[eq, u, F];
  If[FreeQ[r, Indeterminate | ComplexInfinity | DirectedInfinity], Return[r]];
  occ = pdFirstArgs[eq, u];
  If[occ === $Failed || Length[occ] =!= Length[vars], Return[r]];
  fixed = Select[Range[Length[vars]], occ[[#]] =!= vars[[#]] &];
  rr = (eq /. Equal -> Subtract) /.
         (h : (u | Derivative[__][u]))[b__] /; Length[{b}] == Length[vars] :>
           h @@ ReplacePart[{b}, Thread[fixed -> vars[[fixed]]]];
  rr = Quiet[rr /. u -> F];
  pdTC[Fold[Limit[#1, vars[[#2]] -> occ[[#2]]] &, rr, fixed], r]
];

(* --- explicit branch ----------------------------------------------------- *)

pdExplicitVerdict[eqs_, F_, u_, vars_, ivp_] := Module[{verds, leaked},
  verds = Table[
    With[{eq = eq},
      If[pdIsCondition[eq, u, vars],
        pdPinnedVerdict[Function[pin, pdCondResid[eq, u, pdPin[F, pin], vars]]],
        pdPinnedVerdict[Function[pin, pdSubst[eq, u, pdPin[F, pin]]]]]],
    {eq, eqs}];
  leaked = pdHasArbFun[F] || !FreeQ[F, C[_Integer]];
  Which[
    MemberQ[verds, "BAD"], If[ivp && leaked, "UNFIT", "BAD"],
    MemberQ[verds, "UNK"], If[ivp && leaked, "UNFIT", "UNK"],
    True, "OK"]
];

(* --- implicit branch: phi(vars, u) == rhs ----------------------------------- *)

pdImplicitVerdict[eqs_, rel_Equal, u_, vars_] := Module[
  {U, R, pdes, order, check},
  pdes = Select[eqs, !pdIsCondition[#, u, vars] &];
  If[Length[pdes] =!= Length[eqs], Return["UNK"]];       (* conditions: not checked here *)
  order = Max[0, Cases[pdes, Derivative[n__][u][___] :> Total[{n}], {0, Infinity}]];
  If[order =!= 1, Return["UNK"]];
  check = Function[pin, Module[{Rp},
    Rp = pdPin[(rel /. Equal -> Subtract) /. u @@ vars -> U, pin];
    Total[Table[Together[(pde /. Equal -> Subtract) /. {
        Derivative[n__][u][Sequence @@ vars] :>
          With[{j = First@Flatten@Position[{n}, 1]},
            -D[Rp, vars[[j]]]/D[Rp, U]],
        u @@ vars -> U}]^2, {pde, pdes}]]]];
  Switch[pdPinnedVerdict[check], "OK", "IMPL", "BAD", "UNK", _, "UNK"]
];

(* --- series branch ----------------------------------------------------------- *)

pdSeriesVerdict[eqs_, body_, u_, vars_] := Module[
  {sums, s0, term, kk, lo, assum, simp, R, modeOK, pdes, conds, cverd, F0, Ft},
  sums = Cases[body, Inactive[Sum][_, {_, _, _}], {0, Infinity}];
  If[Length[sums] =!= 1 || !FreeQ[sums[[1, 1]], Inactive[Sum]] ||
     !FreeQ[body, BesselJ | BesselY | LegendreP | SphericalBesselJ], Return["UNK"]];
  term = sums[[1, 1]];  kk = sums[[1, 2, 1]];  lo = sums[[1, 2, 2]];
  s0 = body /. Inactive[Sum][__] -> 0;
  assum = Element[kk, Integers] && kk >= lo;
  simp[e_] := pdTC[Simplify[e, assum], e];
  R[w_, eq_] := pdSubst[eq, u, Function @@ {vars, w}];
  pdes  = Select[eqs, !pdIsCondition[#, u, vars] &];
  conds = Select[eqs, pdIsCondition[#, u, vars] &];
  (* PDE: the full ansatz carries the inhomogeneity, each mode is homogeneous *)
  modeOK = AllTrue[pdes, pdResidVerdict[simp[R[s0 + term, #]]] === "OK" &&
                         pdResidVerdict[simp[R[term, #] - R[0, #]]] === "OK" &];
  If[!modeOK, Return[If[AnyTrue[pdes, pdResidVerdict[simp[R[term, #] - R[0, #]]] === "BAD" &],
                        "BAD", "UNK"]]];
  cverd = Table[Module[{full, mode, line, j, ends, x, G, T, Tm, m, ip, orth},
      full = pdResidVerdict[simp[R[s0 + term, c]]];
      mode = pdResidVerdict[simp[R[term, c] - R[0, c]]];
      If[full === "OK", "OK",
        (* the expanded condition: project onto the modes *)
        line = pdFirstArgs[c, u];
        If[line === $Failed, "UNK",
          j = Select[Range[Length[vars]], line[[#]] === vars[[#]] &];
          If[Length[j] =!= 1, "UNK",
            x = vars[[j[[1]]]];
            ends = Union@Cases[conds, (u | Derivative[__][u])[a__] /;
                       Length[{a}] == Length[vars] && FreeQ[{a}[[j[[1]]]], x] :>
                       {a}[[j[[1]]]], {0, Infinity}];
            If[Length[ends] =!= 2, "UNK",
              G = R[s0, c];  T = R[term, c] - R[0, c];
              Tm = T /. kk -> m;
              orth = pdTC[Simplify[Integrate[T Tm, {x, ends[[1]], ends[[2]]}],
                       assum && Element[m, Integers] && m >= lo && m != kk], $Failed];
              ip = pdTC[Simplify[Integrate[G T + T^2, {x, ends[[1]], ends[[2]]}], assum],
                       $Failed];
              Which[
                !FreeQ[{orth, ip}, $Failed | Integrate], "UNK",
                pdResidVerdict[orth] =!= "OK", "UNK",
                True, Switch[pdResidVerdict[ip], "OK", "OK", "BAD", "BAD", _, "UNK"]]]]]]],
    {c, conds}];
  Which[MemberQ[cverd, "BAD"], "BAD", MemberQ[cverd, "UNK"], "UNK", True, "SER"]
];

(* --- one branch, top level ------------------------------------------------ *)

pdBranchVerdict[eqs_, br_, u_, vars_, ivp_] := Module[{rules, F},
  rules = If[ListQ[br], br, {br}];
  If[Length[rules] =!= 1, Return["UNK"]];
  If[Head[rules[[1]]] === Equal, Return[pdImplicitVerdict[eqs, rules[[1]], u, vars]]];
  F = pdRuleFun[rules[[1]], vars];
  If[F === $Failed, Return["UNK"]];
  If[!FreeQ[F, Inactive[Sum]], Return[pdSeriesVerdict[eqs, Last[F], u, vars]]];
  pdExplicitVerdict[eqs, F, u, vars, ivp]
];

pdeCheckCode[label_, section_, eqns_, fn_, vars_, opts_] := Module[
  {eqs, u, sol, ivp, verds},
  eqs = Flatten[{eqns}];
  u = pdUSym[fn];
  ivp = AnyTrue[eqs, pdIsCondition[#, u, vars] &];
  sol = Quiet@TimeConstrained[DSolve[eqns, fn, vars, Sequence @@ opts],
                              $pdSolveTimeout, $Aborted];
  If[sol === $Aborted || Head[sol] === DSolve || !ListQ[sol] || sol === {}, Return[2]];
  verds = Quiet[pdBranchVerdict[eqs, #, u, vars, ivp] & /@ sol];   (* probes *)
  Which[
    MemberQ[verds, "BAD"],   1,
    MemberQ[verds, "UNFIT"], 7,
    MemberQ[verds, "UNK"],   6,
    MemberQ[verds, "SER"],   5,
    MemberQ[verds, "IMPL"],  4,
    True,                    0]
];

(* Verify a GIVEN answer (no DSolve call): the prelude's own self-test hook. *)
pdeCheckAnswer[eqns_, fn_, vars_, sol_] := Module[{eqs, u, ivp, verds},
  eqs = Flatten[{eqns}]; u = pdUSym[fn];
  ivp = AnyTrue[eqs, pdIsCondition[#, u, vars] &];
  verds = Quiet[pdBranchVerdict[eqs, #, u, vars, ivp] & /@ sol];
  Which[MemberQ[verds, "BAD"], 1, MemberQ[verds, "UNFIT"], 7, MemberQ[verds, "UNK"], 6,
        MemberQ[verds, "SER"], 5, MemberQ[verds, "IMPL"], 4, True, 0]
];

pdeVerdictName[c_] := Switch[c, 0, "PASS", 1, "FAIL", 2, "UNEVAL", 3, "SKIP",
  4, "PASS-IMPLICIT", 5, "PASS-SERIES", 6, "UNK", 7, "UNFIT", _, "?"];
