(* NDSolveBVP.m -- boundary-value problems for NDSolve by shooting.

   Loaded lazily by NDSolve (src/numerical_calculus/ndsolve.c) when the
   conditions of an ODE system are not all plain initial values
   f^(k)[a] == v at one point a: they sit at two or more points (u[0] == 0,
   u'[1] == 0), or they are not of that plain form (u[0] + u'[0] == 1,
   u[0] == u[1]).  All symbols are qualified in NDSolveBVP`.

   METHOD.  Let a be the start of the range (or the condition point equal to
   it, so u[0] == 0 is recognised against {x, 0., 1.}) and let the state at a
   be every f^(k)[a], k < order of f.  Plain conditions at a fix part of that
   state; the rest is unknown, s.  For a trial s, NDSolve integrates the
   initial-value problem from a over the whole range, and the remaining
   conditions, evaluated on that solution, are the residual r(s).  Newton's
   method drives r to zero: a forward-difference Jacobian (one extra IVP per
   unknown), a linear solve, and a step halved until |r| decreases.  The
   difference step is 10^(-wp/4) relative, well above the IVP's own error so
   the Jacobian is not noise; a linear problem's residual is affine in s, so
   any step gives its Jacobian exactly and Newton ends in one step.

   The IVPs are solved by NDSolve itself with the caller's options, so
   WorkingPrecision, PrecisionGoal, Method, ... apply to every shot and an
   MPFR problem is solved to MPFR precision.  The answer is that of the final
   shot, requested in the caller's form (u or u[x]).

   RESULT.  {"ok", solution}, {"berr", solution, rel} when Newton stops with
   a relative residual above Sqrt[tolerance] (the solution is still returned, as NDSolve
   does, with NDSolve::berr from the C side), or {"bvdim", nunknown, ncond}
   when the counts disagree, or {"fail"}. *)

BeginPackage["NDSolveBVP`"];
Shoot::usage = "NDSolveBVP`Shoot[odes, conds, funcs, {x, a, b}, heads, orders, opts] solves the ODE boundary-value problem by shooting from a.";
Begin["`Private`"];

(* f^(k)[p] for a function head f *)
lit[f_, 0, p_] := f[p];
lit[f_, k_, p_] := Derivative[k][f][p];

(* the points at which a condition evaluates the dependent functions *)
condPoints[c_, heads_] := Union[Cases[c,
   (h_Symbol | Derivative[_][h_Symbol])[p_] /; MemberQ[heads, h] && NumericQ[p] :> p,
   {0, Infinity}, Heads -> True]];

(* plain condition  lit == v  with lit in the state list and v numeric *)
plainValue[Equal[l_, v_], state_] /; MemberQ[state, l] && NumericQ[v] := {l, v};
plainValue[Equal[v_, l_], state_] /; MemberQ[state, l] && NumericQ[v] := {l, v};
plainValue[_, _] := None;

Shoot[odes_List, conds_List, funcs_, {x_, xa_, xb_}, heads_List, orders_List,
      opts___] :=
  Module[{pts, a, state, plain, fixed, rest, unk, m, wp, pg, tol, h, ivp, resid,
          s, r, rel, nr, tolN, J, ds, lam, sTry, rTry, iter = 0, stepN, sol},
    pts = Union @@ (condPoints[#, heads] & /@ conds);
    a = SelectFirst[pts, TrueQ[N[#] == N[xa]] &, xa];
    state = Flatten[Table[lit[heads[[i]], k, a], {i, Length[heads]},
                          {k, 0, orders[[i]] - 1}]];
    plain = plainValue[#, state] & /@ conds;
    (* a plain condition fixes its state entry only once; a second plain
       condition on the same entry stays a residual (and is then over-determined) *)
    fixed = {}; rest = {};
    Do[If[plain[[i]] =!= None && !MemberQ[fixed[[All, 1]], plain[[i, 1]]],
          AppendTo[fixed, plain[[i]]], AppendTo[rest, conds[[i]]]],
       {i, Length[conds]}];
    unk = Complement[state, fixed[[All, 1]]];
    unk = Select[state, MemberQ[unk, #] &];        (* state order *)
    m = Length[unk];
    If[m != Length[rest], Return[{"bvdim", m, Length[rest]}]];
    wp = Lookup[Cases[{opts}, _Rule | _RuleDelayed], WorkingPrecision, MachinePrecision];
    (* kept exact: a machine-real wp would make the difference step a machine
       number, and one machine number in the Jacobian drops every shot to
       machine precision *)
    wp = If[wp === MachinePrecision || !NumericQ[wp], MachinePrecision, Rationalize[wp]];
    (* Newton tolerance = the IVP's PrecisionGoal (Automatic: wp/2, as in
       NDSolve), since no residual is better than the shots it is read from;
       NDSolve::berr only when it is missed by more than half its digits. *)
    pg = Lookup[Cases[{opts}, _Rule | _RuleDelayed], PrecisionGoal, Automatic];
    pg = If[NumericQ[pg] && pg > 0, Round[pg],
            If[wp === MachinePrecision, 8, Round[wp/2]]];
    {tol, h} = If[wp === MachinePrecision, {10.^-pg, 10.^-4},
                  {10^-pg, 10^-Round[wp/4]}];
    (* Newton itself aims past the goal, at the working precision: a shot is
       usually far better than its PrecisionGoal, and a final step costs only
       m + 2 shots.  Once the goal is met an unproductive step just ends the
       iteration (no step halving). *)
    tolN = If[wp === MachinePrecision, 10.^-13, 10^-Round[0.9 wp]];
    (* one shot; a shot NDSolve declines is $Failed INSIDE the Quiet -- an
       unevaluated NDSolve[...] held in a variable would re-run (and re-warn)
       every time the variable is read *)
    ivp[sv_] := Quiet[Replace[NDSolve[Join[odes, Equal @@@ fixed, Thread[unk == sv]],
                                      heads, {x, xa, xb}, opts],
                              Except[{{__Rule}, ___}] -> $Failed]];
    (* {residual vector, its size relative to the condition values}: the
       Newton iteration uses the raw residual (affine in s for a linear
       problem), convergence is judged on the relative one *)
    resid[sv_] := Module[{so = ivp[sv], lhs, rhs, rv},
      If[so === $Failed, Return[$Failed, Module]];
      lhs = (First /@ rest) /. First[so]; rhs = (Last /@ rest) /. First[so];
      rv = lhs - rhs;
      If[!VectorQ[rv, NumericQ], Return[$Failed, Module]];
      {rv, Max[Abs[rv]/(1 + Abs[lhs] + Abs[rhs])]}];
    If[m == 0, Return[{"ok", NDSolve[Join[odes, Equal @@@ fixed], funcs, {x, xa, xb}, opts]}]];
    s = If[wp === MachinePrecision, ConstantArray[0., m], ConstantArray[N[0, wp], m]];
    r = resid[s];
    If[r === $Failed, Return[{"fail"}]];
    {r, rel} = r; nr = Norm[r];
    While[iter < 40 && rel > tolN,
      iter++;
      J = Transpose[Table[
            Module[{hj = h (1 + Abs[s[[j]]]), sj = s, rj},
              sj[[j]] += hj; rj = resid[sj];
              If[rj === $Failed, Return[{"fail"}, Module]];
              (First[rj] - r)/hj], {j, m}]];
      If[!MatrixQ[J, NumericQ], Return[{"fail"}]];
      ds = Quiet[LinearSolve[J, -r]];
      If[!VectorQ[ds, NumericQ], Return[{"fail"}]];
      (* damped step: halve until the residual decreases *)
      lam = 1; rTry = $Failed;
      While[lam >= If[rel <= tol, 1, 2^-12],
        sTry = s + lam ds; rTry = resid[sTry];
        If[rTry =!= $Failed && Norm[First[rTry]] < nr, Break[]];
        lam /= 2];
      If[rTry === $Failed || Norm[First[rTry]] >= nr, Break[]];   (* stalled *)
      stepN = Norm[lam ds];
      s = sTry; {r, rel} = rTry; nr = Norm[r];
      (* at the shots' noise floor: the step no longer moves s *)
      If[stepN <= tolN (1 + Norm[s]), Break[]]];
    sol = NDSolve[Join[odes, Equal @@@ fixed, Thread[unk == s]], funcs, {x, xa, xb}, opts];
    If[!MatchQ[sol, {{__Rule}, ___}], Return[{"fail"}]];
    If[rel > Sqrt[tol], {"berr", sol, rel}, {"ok", sol}]];

End[];
EndPackage[];
