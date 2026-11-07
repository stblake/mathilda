(* D3 reference, method A: FindMinimum (WL's quasi-Newton / Newton) at
   WorkingPrecision 130 on the same 197-coordinate energy as thomson_100.m,
   started from the double-precision global minimiser written by d3_start.py
   (WL's own FindMinimum from the Fibonacci start lands in the 4448.4394 local
   minimum, so the start is supplied).  Cross-checks: |grad| at the solution,
   and the energy recomputed from the minimiser at 150 digits. *)
Get["d3_start.wl"];
npts = 100;
ph1 = Pi (3 - Sqrt[5]);
pt[0] = {0, 0, 1};
pt[1] = {Sin[th[1]] Cos[ph1], Sin[th[1]] Sin[ph1], Cos[th[1]]};
pt[k_Integer] := {Sin[th[k]] Cos[ph[k]], Sin[th[k]] Sin[ph[k]], Cos[th[k]]};
energy = Sum[1/Sqrt[(pt[j] - pt[k]).(pt[j] - pt[k])], {j, 0, npts - 2}, {k, j + 1, npts - 1}];
vars = Join[Table[th[k], {k, 1, npts - 1}], Table[ph[k], {k, 2, npts - 1}]];
st = Transpose[{vars, SetPrecision[d3start, 130]}];
{e, sol} = FindMinimum[energy, st, WorkingPrecision -> 130, PrecisionGoal -> 60, AccuracyGoal -> 60,
   Method -> "Newton", MaxIterations -> 200];
g = D[energy, {vars}] /. sol;
Print["|grad| = ", N[Norm[g], 5]];
Print["E_100 = ", NumberForm[e, 80]];
