(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 230 -- SIAM 100-Digit Challenge, sealed holdout B for Problem 5
   ("Hard Work"): ten more best uniform (minimax) polynomial approximations,
   stated -- as in experiment 170 -- as the DISCRETE minimax problem a user
   writes: NMinimize over the polynomial's coefficients of the max of
   |f - p| over sample points.  Each row is scored against the discrete
   minimax value E itself (truth.tsv).

   Real, on M equispaced points x_k = a + (b-a) k/(M-1), k = 0..M-1,
   p(x) = c0 + c1 x + ... + cn x^n (n+1 real unknowns):
     bc5.1   e^x           on [-1,1],  n = 2, M = 101
     bc5.2   |x|           on [-1,1],  n = 4, M = 201   (non-smooth f)
     bc5.3   sqrt(x)       on [0,1],   n = 3, M = 101   (endpoint branch)
     bc5.4   1/(1+25x^2)   on [-1,1],  n = 6, M = 201   (Runge)
     bc5.5   Gamma(x)      on [1,2],   n = 3, M = 101
     bc5.6   arctan(x)     on [-1,1],  n = 5, M = 201
     bc5.7   J0(10x)       on [0,1],   n = 8, M = 201   (9 unknowns)
   Complex, on the M-th roots of unity z_k = exp(2 pi i k/M),
   p(z) = Sum (a_j + i b_j) z^j (2(n+1) real unknowns, no symmetry reduction):
     bc5.8   e^z           n = 3, M = 64
     bc5.9   sqrt(z+2)     n = 2, M = 96
     bc5.10  Gamma(z+2)    n = 4, M = 128  (10 unknowns; pole at -2)

   The objective is NON-SMOOTH (a max of M moduli); at the optimum n+2 (real)
   or more (complex) of them are equal.  Rows: NMinimize at machine precision
   ("bc5.k machine") and with WorkingPrecision -> 60 ("bc5.k wp50").
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize", "Gamma", "BesselJ"}];

cv = {c0, c1, c2, c3, c4, c5, c6, c7, c8};
av = {a0, a1, a2, a3, a4};  bv = {b0, b1, b2, b3, b4};

realObj[f_, {a_, b_}, n_, m_, prec_] := Module[{xs, xk, fk},
  xs = a + (b - a) Range[0, m - 1]/(m - 1);
  xk = N[xs, prec]; fk = N[f /@ xs, prec];
  Max[Abs[fk - (c0 + Sum[cv[[j + 1]] xk^j, {j, 1, n}])]]];
cplxObj[f_, n_, m_, prec_] := Module[{zs, zk, fk},
  zs = Exp[2 Pi I Range[0, m - 1]/m];
  zk = N[zs, prec]; fk = N[f /@ zs, prec];
  Max[Abs[fk - ((a0 + I b0) + Sum[(av[[j + 1]] + I bv[[j + 1]]) zk^j, {j, 1, n}])]]];

probs = {
  {"r", Exp, {-1, 1}, 2, 101}, {"r", Abs, {-1, 1}, 4, 201},
  {"r", Sqrt, {0, 1}, 3, 101}, {"r", 1/(1 + 25 #^2) &, {-1, 1}, 6, 201},
  {"r", Gamma, {1, 2}, 3, 101}, {"r", ArcTan, {-1, 1}, 5, 201},
  {"r", BesselJ[0, 10 #] &, {0, 1}, 8, 201},
  {"c", Exp, 3, 64}, {"c", Sqrt[# + 2] &, 2, 96}, {"c", Gamma[# + 2] &, 4, 128}};

Do[With[{pr = probs[[k]]},
    With[{F = If[pr[[1]] === "r", realObj[pr[[2]], pr[[3]], pr[[4]], pr[[5]], #] &,
                                  cplxObj[pr[[2]], pr[[3]], pr[[4]], #] &],
          vars = If[pr[[1]] === "r", Take[cv, pr[[4]] + 1],
                    Join[Take[av, pr[[3]] + 1], Take[bv, pr[[3]] + 1]]]},
      vrun["bc5." <> ToString[k] <> " machine",
           First[NMinimize[F[MachinePrecision], vars]], 10, 1, 120];
      vrun["bc5." <> ToString[k] <> " wp50",
           First[NMinimize[F[60], vars, WorkingPrecision -> 60]], 48, 1, 120]]],
  {k, 10}];
