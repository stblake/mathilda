(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 238 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   problem 3 (Thomson problem, N = 100): ten other minimal-energy
   configurations, each a FindMinimum (or FindMaximum) over many coordinates
   from a deterministic start, where landing in a non-global local optimum is
   a WRONG answer.

     bd3.1   Thomson N = 5   (triangular bipyramid)       E = 6.474691494...
     bd3.2   Thomson N = 12  (icosahedron)                E = 49.16525305...
     bd3.3   30 points on the unit CIRCLE, Coulomb        E = 505.1594202...
     bd3.4   Thomson N = 24  (snub cube)                  E = 223.3470740...
     bd3.5   Thomson N = 32                               E = 412.2612746...
     bd3.6   logarithmic energy Sum -Log|x_j-x_k|, N = 20 E = -54.01112997...
     bd3.7   Riesz s = 2 energy Sum |x_j-x_k|^-2, N = 20  E = 133.9369785...
     bd3.8   Lennard-Jones cluster LJ13 in R^3            E = -44.32680141...
     bd3.9   Thomson N = 60                               E = 1543.830400...
     bd3.10  MAXIMAL Sum |x_j-x_k| over N = 30 sphere points (FindMaximum)
                                                          S = 597.7729699...

   Sphere rows: the parametrisation of experiment 178 -- point 0 held at the
   north pole, phi_1 held at Pi (3 - Sqrt[5]), unknowns theta_1..theta_{N-1},
   phi_2..phi_{N-1}; start = the Fibonacci lattice theta_k = ArcCos[1 - 2k/N],
   phi_k = k Pi (3 - Sqrt[5]).  bd3.3: angles phi_1..phi_29 (phi_0 = 0) from
   phi_k = 2 Pi (k/30)^2.  bd3.8: 39 Cartesian coordinates, start = an atom at
   the origin plus 12 Fibonacci points (theta_k = ArcCos[1 - (2k+1)/12]) at
   radius 1.1.  From each of these starts scipy BFGS reaches the global
   optimum (siam100_refs/bd3_ref.py, with random-start checks).
   Rows: machine precision; WorkingPrecision -> 60 (default goals: the optimum
   VALUE carries about twice the digits of the optimiser).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "FindMaximum"}];

ph1 = Pi (3 - Sqrt[5]);
spt[0] = {0, 0, 1};
spt[1] = {Sin[th[1]] Cos[ph1], Sin[th[1]] Sin[ph1], Cos[th[1]]};
spt[k_Integer] := {Sin[th[k]] Cos[ph[k]], Sin[th[k]] Sin[ph[k]], Cos[th[k]]};
kern["coul"][q_] := 1/Sqrt[q];
kern["log"][q_] := -Log[q]/2;
kern["r2"][q_] := 1/q;
kern["dist"][q_] := Sqrt[q];
senergy[n_, kd_] := Sum[kern[kd][(spt[j] - spt[k]).(spt[j] - spt[k])],
                        {j, 0, n - 2}, {k, j + 1, n - 1}];
sstart[n_, prec_] := Join[Table[{th[k], N[ArcCos[1 - 2 k/n], prec]}, {k, 1, n - 1}],
                          Table[{ph[k], N[k ph1, prec]}, {k, 2, n - 1}]];
sphereCase[p_, n_, kd_] := Module[{e = senergy[n, kd], s1 = sstart[n, MachinePrecision],
                                   s60 = sstart[n, 60], lab = "bd3." <> ToString[p]},
  If[kd === "dist",
    vrun[lab <> " machine", First[FindMaximum[e, s1]], 10, 1, 60];
    vrun[lab <> " wp50", First[FindMaximum[e, s60, WorkingPrecision -> 60]], 48, 1, 120],
    vrun[lab <> " machine", First[FindMinimum[e, s1]], 10, 1, 60];
    vrun[lab <> " wp50", First[FindMinimum[e, s60, WorkingPrecision -> 60]], 48, 1, 120]]];

sphereCase[1, 5, "coul"];
sphereCase[2, 12, "coul"];

(* bd3.3: 30 points on the unit circle *)
cpt[0] = {1, 0}; cpt[k_Integer] := {Cos[al[k]], Sin[al[k]]};
cenergy = Sum[1/Sqrt[(cpt[j] - cpt[k]).(cpt[j] - cpt[k])], {j, 0, 28}, {k, j + 1, 29}];
vrun["bd3.3 machine", First[FindMinimum[cenergy, Table[{al[k], N[2 Pi (k/30)^2]}, {k, 1, 29}]]],
     10, 1, 60];
vrun["bd3.3 wp50", First[FindMinimum[cenergy, Table[{al[k], N[2 Pi (k/30)^2, 60]}, {k, 1, 29}],
                                     WorkingPrecision -> 60]], 48, 1, 120];

sphereCase[4, 24, "coul"];
sphereCase[5, 32, "coul"];
sphereCase[6, 20, "log"];
sphereCase[7, 20, "r2"];

(* bd3.8: Lennard-Jones 13 *)
lj[i_] := {x[i], y[i], z[i]};
ljenergy = Sum[With[{q = (lj[i] - lj[j]).(lj[i] - lj[j])}, 4 (q^-6 - q^-3)], {i, 0, 11}, {j, i + 1, 12}];
ljp[0] = {0, 0, 0};
ljp[k_] := With[{t = ArcCos[1 - (2 (k - 1) + 1)/12], f = (k - 1) ph1},
                 11/10 {Sin[t] Cos[f], Sin[t] Sin[f], Cos[t]}];
ljstart[prec_] := Flatten[Table[Transpose[{lj[i], N[ljp[i], prec]}], {i, 0, 12}], 1];
vrun["bd3.8 machine", First[FindMinimum[ljenergy, ljstart[MachinePrecision]]], 10, 1, 60];
vrun["bd3.8 wp50", First[FindMinimum[ljenergy, ljstart[60], WorkingPrecision -> 60]], 48, 1, 120];

sphereCase[9, 60, "coul"];
sphereCase[10, 30, "dist"];
