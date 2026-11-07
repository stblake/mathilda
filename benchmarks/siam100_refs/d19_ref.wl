(* Reference B for the finite proxy of SIAM100 App. D problem 19 (exp 194):
   an independent DUALITY CERTIFICATE, in exact/80-digit Wolfram arithmetic.
   With e_e(v) = v_e^T Kq v_e the element energies,
       J(c) = max_v [2 F.v - Sum_e c_e e_e(v)]           (compliance)
   so for any v,
       min_c J(c) >= LB(v) = 2 F.v - max_{c feasible} Sum_e c_e e_e(v),
   where the inner max is a fractional knapsack (c = 100 on the largest
   energies, 1 on the smallest, one cell fractional) -- no optimiser needed.
   UB = J(c_star) for the design c_star from reference A (siam100_refs/d19_ref.py,
   file d19_cstar_n<n>.txt), recomputed here by LinearSolve; LB at
   v = U(c_star).  LB == UB to the printed digits proves min J_n = UB (weak
   duality), independently of how c* was found.
   Run: wolframscript -file d19_ref.wl 8      (or 4) *)
n = ToExpression[Last[$ScriptCommandLine]];
prec = 90;
cstar = SetPrecision[ToExpression /@ StringSplit[Import[
     FileNameJoin[{DirectoryName[$InputFileName], "d19_cstar_n" <> ToString[n] <> ".txt"}], "Text"]], prec];
kq = {{4, -1, -2, -1}, {-1, 4, -1, -2}, {-2, -1, 4, -1}, {-1, -2, -1, 4}}/6;
nodeIdx[i_, j_] := If[1 <= i <= n - 1 && 1 <= j <= n - 1, (i - 1) (n - 1) + j, 0];
elems = Flatten[Table[{nodeIdx[a, b], nodeIdx[a + 1, b], nodeIdx[a + 1, b + 1], nodeIdx[a, b + 1]},
    {a, 0, n - 1}, {b, 0, n - 1}], 1];
nInt = (n - 1)^2;
F = ConstantArray[1/n^2, nInt];
stiff[c_] := Module[{K = ConstantArray[0, {nInt, nInt}]},
  Do[With[{idx = elems[[e]]}, Do[If[idx[[r]] > 0 && idx[[s]] > 0,
       K[[idx[[r]], idx[[s]]]] += c[[e]] kq[[r, s]]], {r, 4}, {s, 4}]], {e, Length[elems]}]; K];
U = LinearSolve[stiff[cstar], F];
UB = F . U;
energy[v_] := Table[With[{ve = Table[If[k > 0, v[[k]], 0], {k, elems[[e]]}]}, ve . kq . ve], {e, Length[elems]}];
en = energy[U];
(* fractional knapsack: budget Sum c = 101/2 n^2, 1 <= c <= 100 *)
ord = Ordering[en, All, Greater];
budget = 101/2 n^2 - n^2;                 (* extra above c = 1 everywhere *)
cmax = ConstantArray[1, n^2];
Do[With[{add = Min[99, budget]}, cmax[[e]] += add; budget -= add], {e, ord}];
LB = 2 F . U - cmax . en;
Print["n = ", n, "  feasible c*: ", {Min[cstar] >= 1, Max[cstar] <= 100, Total[cstar] - 101/2 n^2}];
Print["UB = J(c*) = ", N[UB, 75]];
Print["LB         = ", N[LB, 75]];
Print["UB - LB    = ", N[UB - LB, 5]];
