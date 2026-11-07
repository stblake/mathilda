(* Holdout A, family hd3, method A: closed forms where they exist
   (N = 5 triangular bipyramid, N = 12 icosahedron, N = 30 regular polygon),
   otherwise FindMinimum Method -> "Newton" at WorkingPrecision 120 on the
   benchmark's symbolic energy, started from the double minimiser that
   hd3_ref.py wrote (hd3_start_<k>.wl), with the gradient norm and the energy
   at the minimiser printed.  Run: wolframscript -file hd3_ref.wl *)
ph1 = Pi (3 - Sqrt[5]);
pt[0] = {0, 0, 1};
pt[1] = {Sin[th[1]] Cos[ph1], Sin[th[1]] Sin[ph1], Cos[th[1]]};
pt[k_Integer] := {Sin[th[k]] Cos[ph[k]], Sin[th[k]] Sin[ph[k]], Cos[th[k]]};
gC[r2_] := 1/Sqrt[r2]; gR2[r2_] := 1/r2; gLog[r2_] := -Log[r2]/2; gDist[r2_] := -Sqrt[r2];
energyS[n_, g_] := Sum[g[(pt[j] - pt[k]).(pt[j] - pt[k])], {j, 0, n - 2}, {k, j + 1, n - 1}];
Print["1 closed ", NumberForm[N[1/2 + 3 Sqrt[2] + Sqrt[3], 80], 75]];
Print["2 closed ", NumberForm[N[3 + 30/Sqrt[2 - 2/Sqrt[5]] + 30/Sqrt[2 + 2/Sqrt[5]], 80], 75]];
Print["8 closed ", NumberForm[N[30/4 Sum[Csc[Pi k/30], {k, 1, 29}], 80], 75]];
probs = <|3 -> {8, gC}, 4 -> {32, gC}, 5 -> {20, gR2}, 6 -> {30, gLog}, 7 -> {16, gDist},
          9 -> {50, gC}, 10 -> {72, gC}|>;
Do[With[{n = probs[k][[1]], g = probs[k][[2]]},
  Get["hd3_start_" <> ToString[k] <> ".wl"];
  Module[{en = energyS[n, g], vars, e, sol, gr},
    vars = Join[Table[th[j], {j, 1, n - 1}], Table[ph[j], {j, 2, n - 1}]];
    {e, sol} = FindMinimum[en, Transpose[{vars, SetPrecision[hstart, 120]}],
       WorkingPrecision -> 120, PrecisionGoal -> 60, AccuracyGoal -> 60,
       Method -> "Newton", MaxIterations -> 200];
    gr = D[en, {vars}] /. sol;
    Print[k, " FindMinimum ", NumberForm[e, 75], "   |grad| = ", N[Norm[gr], 3],
          "  E recomputed at 150: ", NumberForm[N[en /. SetPrecision[sol, 150], 80], 75]]]],
  {k, {3, 4, 5, 6, 7, 9, 10}}];
