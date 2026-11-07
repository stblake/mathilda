(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 208 -- SIAM 100-digit HOLDOUT A, family d3 (minimal point
   energies; original: experiment 178, the Thomson problem for N = 100).

   min over N points x_j of Sum_{j<k} g(|x_j - x_k|):
     hd3.1   sphere, N = 5,  g = 1/r            (= 1/2 + 3 Sqrt[2] + Sqrt[3])
     hd3.2   sphere, N = 12, g = 1/r            (icosahedron, closed form)
     hd3.3   sphere, N = 8,  g = 1/r            (square antiprism)
     hd3.4   sphere, N = 32, g = 1/r
     hd3.5   sphere, N = 20, g = 1/r^2          (Riesz s = 2)
     hd3.6   sphere, N = 30, g = -log r         (logarithmic energy, Whyte)
     hd3.7   sphere, N = 16, g = -r             (max sum of distances, Fejes Toth)
     hd3.8   circle, N = 30, g = 1/r            (= (N/4) Sum_k Csc[Pi k/N])
     hd3.9   sphere, N = 50, g = 1/r
     hd3.10  sphere, N = 72, g = 1/r            (the harder member)
   The answer is the GLOBAL minimum value (confirmed by random-start BFGS
   searches in siam100_refs/hd3_ref.py); a run that stops in another local
   minimum has given a wrong answer.

   Idiomatic statement, exactly as experiment 178: FindMinimum of the
   symbolic energy over spherical coordinates, point 0 held at the north pole
   and phi_1 held at Pi (3 - Sqrt[5]) (2N - 3 unknowns), started from the
   Fibonacci lattice theta_k = ArcCos[1 - 2k/N], phi_k = k Pi (3 - Sqrt[5]).
   The circle: angles t_1..t_{N-1} (t_0 = 0 held), start t_k = 2 Pi (k/N)^(3/2).
   Rows: machine precision; WorkingPrecision -> 60 from the same start.
   True values: closed forms (hd3.1, .2, .8) / mpmath Newton at 120 digits
   == wolframscript FindMinimum Method -> "Newton" at WP 120 (hd3_ref.py and .wl).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum"}];

ph1 = Pi (3 - Sqrt[5]);
pt[0] = {0, 0, 1};
pt[1] = {Sin[th[1]] Cos[ph1], Sin[th[1]] Sin[ph1], Cos[th[1]]};
pt[k_Integer] := {Sin[th[k]] Cos[ph[k]], Sin[th[k]] Sin[ph[k]], Cos[th[k]]};
pc[0] = {1, 0};
pc[k_Integer] := {Cos[t[k]], Sin[t[k]]};
(* g as a function of the squared distance *)
gC[r2_] := 1/Sqrt[r2];  gR2[r2_] := 1/r2;  gLog[r2_] := -Log[r2]/2;  gDist[r2_] := -Sqrt[r2];
energyS[n_, g_] := Sum[g[(pt[j] - pt[k]).(pt[j] - pt[k])], {j, 0, n - 2}, {k, j + 1, n - 1}];
startS[n_, pr_] := Join[Table[{th[k], N[ArcCos[1 - 2 k/n], pr]}, {k, 1, n - 1}],
                        Table[{ph[k], N[k ph1, pr]}, {k, 2, n - 1}]];
energyC[n_] := Sum[gC[(pc[j] - pc[k]).(pc[j] - pc[k])], {j, 0, n - 2}, {k, j + 1, n - 1}];
startC[n_, pr_] := Table[{t[k], N[2 Pi (k/n)^(3/2), pr]}, {k, 1, n - 1}];

probs = {{5, gC, 60}, {12, gC, 60}, {8, gC, 60}, {32, gC, 60}, {20, gR2, 60},
         {30, gLog, 60}, {16, gDist, 60}, {30, "circle", 60}, {50, gC, 120}, {72, gC, 120}};
Do[With[{n = probs[[i, 1]], g = probs[[i, 2]], cap = probs[[i, 3]]},
   With[{en = If[g === "circle", energyC[n], energyS[n, g]],
         st = If[g === "circle", startC[n, MachinePrecision], startS[n, MachinePrecision]],
         st60 = If[g === "circle", startC[n, 60], startS[n, 60]]},
     vrun["hd3." <> ToString[i] <> " machine", First[FindMinimum[en, st]], 10, 1, cap];
     vrun["hd3." <> ToString[i] <> " wp50",
          First[FindMinimum[en, st60, WorkingPrecision -> 60]], 48, 1, cap]]],
  {i, 1, 10}];
