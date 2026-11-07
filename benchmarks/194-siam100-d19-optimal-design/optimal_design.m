(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 194 -- SIAM 100-Digit Challenge, Appendix D, Problem 19
   ("optimal design").

     -div(c(x) grad u) = 1 on the unit square, u = 0 on the boundary.  What
     is the infimum of Int u over measurable c that equal 1 on half the area
     of the square and 100 on the rest?

   A RESEARCH PROBLEM: no digit of the answer is established.  (The book says
   "infimum"; the authors' web version of Appendix D says "supremum" -- the
   benchmark follows the printed book.)  For one load in a scalar problem the
   relaxation (homogenisation) replaces c by a laminate along the flux, whose
   conductivity there is the arithmetic mean; the relaxed problem is the
   convex "variable thickness" problem
        min { J(c) = Int u : 1 <= c <= 100, mean(c) = 101/2 }.

   WHAT IS ASKED HERE is a well-defined finite version of that relaxed
   problem (the "proxy"): n x n uniform squares, n = 8, c constant on each
   (64 unknowns), bilinear (Q1) finite elements for u (49 interior nodes),
       J_8(c) = F . U,   K(c) U = F,   F_i = h^2,
   element stiffness c_e (1/6)[[4,-1,-2,-1],[-1,4,-1,-2],[-2,-1,4,-1],
   [-1,-2,-1,4]] (corners counter-clockwise), minimised over 1 <= c_e <= 100,
   Sum c_e = 50.5 n^2.  J_8 is convex in c, so the minimum VALUE is unique:
       min J_8 = 5.38244792399511293370464160567...e-4     (truth.tsv)
   The optimum mixes c = 1 (12 cells), c = 100 (8 cells) and six interior
   values (44 cells) -- the discrete trace of the relaxed microstructure.
   It is the infimum of the proxy, not of the problem.

   Idiomatic statement: FindMinimum with box and equality constraints over the
   n^2 cell values, the objective a black box built with LinearSolve (glue:
   assembling K).  Rows:
     d19.machine          n = 8 (64 unknowns), machine precision
     d19.machine-reduced  n = 4 (16 unknowns; min J_4 = 5.0626006269...e-4)
     d19.wp50             n = 4 with WorkingPrecision -> 60 (n = 8 at 60
                          digits is out of reach of every system's built-in
                          optimiser within the cap)  Python: scipy.optimize.minimize SLSQP (machine);
   mpmath has no constrained optimiser -> wp50 SKIP.  Maxima: fmin_cobyla
   (derivative-free, constrained) for the machine row; no bigfloat
   constrained optimiser -> wp50 SKIP.  FriCAS: no constrained optimiser ->
   SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "LinearSolve"}];

kq = {{4, -1, -2, -1}, {-1, 4, -1, -2}, {-2, -1, 4, -1}, {-1, -2, -1, 4}}/6;
nodeIdx[n_, i_, j_] := If[1 <= i <= n - 1 && 1 <= j <= n - 1, (i - 1) (n - 1) + j, 0];
elems[n_] := elems[n] = Flatten[Table[{nodeIdx[n, a, b], nodeIdx[n, a + 1, b],
     nodeIdx[n, a + 1, b + 1], nodeIdx[n, a, b + 1]}, {a, 0, n - 1}, {b, 0, n - 1}], 1];

(* compliance Int u_h = F.U of the n x n design cv (a numeric black box) *)
compliance[n_, cv_ /; VectorQ[cv, NumericQ]] := Module[{K, F, e, r, s, idx, m = (n - 1)^2},
  K = ConstantArray[0 First[cv], {m, m}];
  Do[idx = elems[n][[e]];
     Do[If[idx[[r]] > 0 && idx[[s]] > 0,
          K[[idx[[r]], idx[[s]]]] += cv[[e]] kq[[r, s]]], {r, 4}, {s, 4}],
     {e, n^2}];
  F = ConstantArray[1/n^2, m];
  F . LinearSolve[K, F]];

designMin[n_, opts___] := Module[{cv = Array[cc, n^2]},
  First[FindMinimum[{compliance[n, cv], And @@ Thread[1 <= cv <= 100],
      Total[cv] == 101/2 n^2},
    Transpose[{cv, ConstantArray[101/2, n^2]}], opts]]];

vrun["d19.machine", designMin[8], 10];
vrun["d19.machine-reduced", designMin[4], 10];
vrun["d19.wp50", designMin[4, WorkingPrecision -> 60, PrecisionGoal -> 50], 48];
