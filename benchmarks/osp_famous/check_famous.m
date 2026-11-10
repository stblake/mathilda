(* check_famous.m -- the local grader of benchmarks/osp_famous (see osp_famous.py).
   Every check returns a string "grade" or "grade: reason". Numeric comparisons
   run at 60 digits against a relative tolerance of 10^-30. *)

ospTol = 10^-30;
ospUnevalHeads = {Integrate, Sum, Product, Limit, DSolve, RSolve, Solve, Reduce, Resolve,
   MinimalPolynomial, LaplaceTransform, FunctionExpand, FullSimplify,
   integrate, sum, product, harmonic};

ospUneval[e_, head_] := Or @@ (Not[FreeQ[e, #]] & /@ ospUnevalHeads) ||
   (head =!= None && Not[FreeQ[e, head]]);

(* the free (non-numeric) symbols of an expression *)
ospFree[e_] := Complement[Union[Cases[{e}, s_Symbol /; Not[NumericQ[s]] &&
     Not[MemberQ[{True, False, List, Rule, Integers, Reals, Complexes, Null}, s]], {0, Infinity}]],
   (* a bound iterator (Sum[n^-n, {n, 1, Infinity}]) is not a free parameter *)
   Union[Cases[{e}, (Sum | Product | Table | Integrate)[_, {b_Symbol, ___}] :> b, {0, Infinity}]]];

ospN[e_] := Quiet[N[e, 60]];

(* |a - b| small relative to b *)
ospClose[a_, b_] := Module[{d = ospN[a - b], s = ospN[b]},
   If[Not[NumericQ[d]] || Not[NumericQ[s]], False,
     Abs[d] <= ospTol Max[1, Abs[s]] || (Abs[d] <= 10^-12 Max[1, Abs[s]] && Precision[d] < 40)]];

(* the points at which a symbol left free after the parameters is sampled *)
ospPts = {41/100, 83/100, 137/100};

(* Do two closed forms agree, once the parameters p are given values?  Free
   symbols left over (the transform variable, the indefinite-integral variable)
   are sampled; `indef` allows a constant difference. *)
ospAgree[ans_, ref_, p_, indef_] := Module[{a = ans /. p, r = ref /. p, fr, pts, ds},
   If[Head[a] === List && Head[r] === List,
     Return[Length[a] == Length[r] && And @@ MapThread[ospAgree[#1, #2, {}, indef] &, {a, r}]]];
   fr = ospFree[{a, r}];
   If[fr === {}, Return[ospClose[a, r]]];
   pts = Table[Thread[fr -> RotateLeft[ospPts, j][[1 ;; Length[fr]]]], {j, 0, 2}];
   If[indef =!= False,
     (* antiderivatives: the derivatives agree (a constant, or a piecewise
        constant from a branch cut, may separate them) *)
     ds = D[a - r, indef];
     And @@ Table[With[{v = ospN[ds /. q]}, NumericQ[v] && Abs[v] < 10^-25], {q, pts}],
     And @@ Table[ospClose[a /. q, r /. q], {q, pts}]]];

(* the heads of an expression that are not arithmetic *)
ospHeads[e_] := Complement[Union[Cases[e, h_Symbol[___] :> h, {0, Infinity}]],
   {Plus, Times, Power, Rational, Complex, List, Hold, Sqrt}];

(* heads that never count as worked out, whatever the input held *)
ospNeverWorkedOut = {Root, RootOf, Integrate, Sum, Limit};

(* unevaluated: a task head left in the answer, other than one the reference
   itself is written with (Sophomore's dream is Sum[n^-n, {n, 1, Infinity}]) *)
ospUnevalBut[e_, head_, ref_] := Module[{hs = Complement[ospUnevalHeads, ospHeads[ref]]},
   Or @@ (Not[FreeQ[e, #]] & /@ hs) ||
   (head =!= None && Not[MemberQ[ospHeads[ref], head]] && Not[FreeQ[e, head]])];

ospCheckValue[ans_, ref_, head_, p_, simp_, Hold[inp_]] := Module[{bad},
   If[ospUnevalBut[ans, head, ref], Return["unsolved: left unevaluated"]];
   If[ans === $Failed || ans === Null, Return["error: no answer"]];
   If[Not[ospAgree[ans, ref, p, mrIndef]],
     Return[If[NumericQ[ospN[ans /. p]] || Head[ans] === List, "wrong: differs from the reference",
       "unclear: not numeric at the sample point"]]];
   If[simp,
     bad = Complement[Intersection[ospHeads[ans], Union[ospHeads[Hold[inp]], ospNeverWorkedOut]],
       ospHeads[ref]];
     If[bad =!= {}, Return["unsolved: equal, but keeps " <> ToString[bad]]];
     If[ospNest[ans] > ospNest[ref],
       Return["unsolved: equal, but radicals nested " <> ToString[ospNest[ans]] <>
         " deep against " <> ToString[ospNest[ref]]]];
     If[LeafCount[ans] > 3/2 LeafCount[ref] + 4,
       Return["unsolved: equal, but " <> ToString[LeafCount[ans]] <> " leaves against " <>
         ToString[LeafCount[ref]]]]];
   "solved"];

ospExpectUneval[ans_, head_] := If[ospUneval[ans, head] || (head === None),
   "solved: expected unevaluated", "wrong: claims an answer to an open problem"];

ospCheckMinPoly[ans_, ref_, x_] := Module[{q},
   If[ospUneval[ans, MinimalPolynomial], Return["unsolved: left unevaluated"]];
   If[Not[PolynomialQ[ans, x]], Return["wrong: not a polynomial"]];
   q = Together[ans/ref];
   If[FreeQ[q, x] && q =!= 0, "solved", "wrong: a different polynomial"]];

ospCheckMultiset[ans_, ref_, p_] := Module[{key = Sort[{Re[#], Im[#]} & /@ ospN[#]] &},
   If[ospUneval[ans, Eigenvalues], Return["unsolved"]];
   If[Max[Abs[Flatten[key[ans /. p] - key[ref /. p]]]] < 10^-25, "solved", "wrong"]];

(* ---- differential and difference equations ---- *)

(* arbitrary constants and functions, given values *)
ospConst[e_] := e /. {C[i_][z__] :> Plus @@ ((i/3 + Sin[#] + #^2/(5 i)) & /@ {z}),
     C[i_] :> (i + 2)/7 + 3/(10 i)};

(* a solution {f[v..] -> e} or {f -> Function[..]}, as a pure function *)
ospSolFun[ans_, f_, vars_] := Module[{r},
   r = Cases[ans, (f[___] -> _) | (f -> _), {0, Infinity}];
   If[r === {}, Return[$Failed]];
   r = First[r];
   If[MatchQ[r, f -> _], r[[2]],
     Function @@ {vars, r[[2]] /. Thread[List @@ r[[1]] -> vars]}]];

ospResid[eq_] := Which[eq === True, 0, eq === False, 1,
   Head[eq] === Equal, eq[[1]] - eq[[2]], True, eq];

(* radical nesting depth: Sqrt[a + Sqrt[b]] is 2 *)
ospNest[e_?AtomQ] := 0;
ospNest[Power[b_, q_Rational]] := 1 + ospNest[b];
ospNest[e_] := Max[0, Max[ospNest /@ (List @@ e)]];

(* is every residual zero at the sample points? *)
ospZeroAt[res_List, pts_] := And @@ Flatten[Table[
     (Function[v, NumericQ[v] && Abs[v] < 10^-25]) /@ ospN[res /. q], {q, pts}]];

(* approached from inside: a condition at a removable singularity (y(0) of sin x/x) *)
ospZeroNear[res_List, vars_] := And @@ ((Function[v, NumericQ[v] && Abs[v] < 10^-8]) /@
     ospN[res]);

ospCheckDE[op_, ans_, f_, vars_, eqs_Hold, conds_Hold, ref_Hold, p_, order_] :=
  Module[{fun, eq2, cd2, pts, nconst, implicit, n0, vals, refsol, reffun},
   If[ospUneval[ans, None] && FreeQ[ans, Rule], Return["unsolved: left unevaluated"]];
   fun = ospSolFun[ans, f, vars];
   If[fun === $Failed,
     implicit = Cases[ans, _Equal, {0, Infinity}];
     If[implicit =!= {} && order <= 2 && op === "solve_ode" && conds === Hold[{}],
       Return[ospCheckImplicit[First[implicit], f, vars[[1]], eqs, p, order]]];
     Return["unclear: no explicit solution to substitute"]];
   nconst = Length[Union[Cases[{fun}, C[_], {0, Infinity}]]];
   If[conds === Hold[{}] && nconst < order && op =!= "solve_pde",
     Return["wrong: " <> ToString[nconst] <> " constants, the general solution needs " <>
       ToString[order]]];
   fun = ospConst[fun] /. p;
   If[op === "solve_recurrence",
     eq2 = ReleaseHold[eqs /. f -> fun /. p];
     cd2 = ReleaseHold[conds /. f -> fun /. p];
     n0 = Max[0, Min[Cases[conds, f[k_Integer] :> k, {0, Infinity}] /. {} -> {0}]];
     If[Not[ospZeroAt[ospResid /@ eq2, Table[{vars[[1]] -> n}, {n, n0, n0 + 10}]]],
       Return["wrong: does not satisfy the recurrence"]];
     If[Not[ospZeroAt[ospResid /@ cd2, {{}}]], Return["wrong: does not meet the conditions"]];
     Return["solved"]];
   pts = Table[Thread[vars -> (RotateLeft[{41/100, 63/100, 7/10, 29/100}, j][[1 ;; Length[vars]]])],
     {j, 0, 2}];
   eq2 = ReleaseHold[eqs /. f -> fun /. p];
   If[Not[ospZeroAt[ospResid /@ eq2, pts]], Return["wrong: does not satisfy the equation"]];
   cd2 = ReleaseHold[conds /. f -> fun /. p];
   If[cd2 =!= {} && Not[ospZeroAt[ospResid /@ cd2, pts]],
     (* a condition at a removable singularity (y'(0) of Sin[x]/x): approached
      * from inside, as the one-sided limit of the shifted solution *)
     cd2 = ReleaseHold[conds /. f -> (Function @@ {vars, (fun @@ vars) /.
             Thread[vars -> vars + ospEps]}) /. p];
     cd2 = Quiet[Limit[#, ospEps -> 0, Direction -> "FromAbove"]] & /@ (ospResid /@ cd2);
     If[Not[ospZeroNear[cd2, vars]],
       Return["wrong: does not meet the conditions"]]];
   "solved"];

(* F(x, y) == 0, first order: y' by implicit differentiation, C solved for *)
ospCheckImplicit[rel_, f_, x_, eqs_Hold, p_, order_: 1] := Module[{Y, F, y1, y2, cs, res, pts},
   (* F(x, Y) == 0: Y' = -F_x/F_Y, Y'' = d/dx Y' along the curve.  The equations
    * are released FIRST: they spell y' as D[y[x], x], and replacing y[x] -> Y
    * inside the hold would turn it into D[Y, x] = 0. *)
   F = ospResid[rel] /. f[x] -> Y /. p;
   cs = Union[Cases[{F}, C[_], {0, Infinity}]];
   F = F /. Thread[cs -> Table[23/10 + i/7, {i, Length[cs]}]];
   y1 = -D[F, x]/D[F, Y];
   y2 = D[y1, x] + D[y1, Y] y1;
   res = ReleaseHold[eqs] /. p;
   res = res /. {Derivative[2][f][x] -> y2, Derivative[1][f][x] -> y1} /. f[x] -> Y;
   res = ospResid /@ res;
   pts = {{x -> 31/100, Y -> 47/100}, {x -> 13/10, Y -> 21/100}, {x -> 7/10, Y -> 9/10}};
   If[ospZeroAt[res, pts], "solved", "wrong: the implicit solution does not satisfy the equation"]];

(* ---- equations ---- *)

ospUnwrap[e_] := e /. ConditionalExpression[a_, _] :> a;

ospCheckSolve[ans0_, eqs_Hold, ref_, vars_, p_] := Module[{ans = ospUnwrap[ans0], sols, refs,
    par, inst, rinst, ok, tup},
   If[ospUneval[ans, None], Return["unsolved: left unevaluated"]];
   If[Not[MatchQ[ans, {___List}]], Return["unclear: not a list of solutions"]];
   sols = Select[#, MatchQ[#, _Rule] &] & /@ ans;
   refs = Cases[#, Equal[a_, b_] :> (a -> b)] & /@ ref;
   par = Union[Cases[{ans}, C[_], {0, Infinity}]];
   inst[s_, c_] := s /. Thread[par -> c] /. p;
   (* every answer satisfies the equations *)
   ok = And @@ Flatten[Table[
        With[{r = ReleaseHold[eqs /. inst[s, Table[c, {Length[par]}]]]},
          (Function[v, NumericQ[v] && Abs[v] < 10^-20]) /@ ospN[ospResid /@ r]],
        {s, sols}, {c, {0, 1, 2}}]];
   If[Not[ok], Return["wrong: an answer does not satisfy the equations"]];
   tup[s_] := ospN[vars /. s];
   rinst = Table[tup[r /. k -> 0 /. p], {r, refs}];
   If[par === {} && FreeQ[ref, k],
     If[Length[sols] =!= Length[refs], Return["wrong: " <> ToString[Length[sols]] <>
         " solutions, the reference has " <> ToString[Length[refs]]]]];
   inst = Flatten[Table[tup[inst[s, Table[c, {Length[par]}]]], {s, sols}, {c, -2, 3}], 1];
   If[And @@ Table[Or @@ (Function[t, Max[Abs[Flatten[{t - r}]]] < 10^-20] /@ inst), {r, rinst}],
     "solved", "wrong: a reference solution is missing"]];
