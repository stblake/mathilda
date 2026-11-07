(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 180 -- SIAM 100-Digit Challenge, Appendix D Problem 5
   (Trefethen): a Chebyshev polynomial of a matrix.

       A = the 48 x 48 Toeplitz matrix with -1 on the first subdiagonal, +1 on
           the main diagonal and the first three superdiagonals;
       min over monic p of degree 8 of ||p(A)||_2  =  1766.3135...

   The objective c -> ||A^8 + Sum_{j<8} c_j A^j||_2 is CONVEX but NONSMOOTH:
   at the minimiser the largest singular value is double (1766.31...,
   the next is 1756.80), so a smooth optimiser stalls on the crease and
   reports a value that is too high (scipy's Nelder-Mead stops at 1766.33,
   BFGS from 0 at 2049).  The reference is a KKT/Newton solve of the
   crease conditions (siam100_refs/d5_ref.wl and .py).

   Idiomatic statement: FindMinimum of Norm[p(A)] (the matrix 2-norm) over
   the 8 lower coefficients, from c = 0 (p = x^8), at machine precision and
   at WorkingPrecision -> 60.  The powers A^j are precomputed exactly (glue).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "Norm"}];

nT = 48;
AT = Table[Which[i == j, 1, j == i - 1, -1, 1 <= j - i <= 3, 1, True, 0], {i, nT}, {j, nT}];
pw = Table[MatrixPower[AT, j], {j, 0, 8}];
obj[c0_?NumericQ, c1_?NumericQ, c2_?NumericQ, c3_?NumericQ,
    c4_?NumericQ, c5_?NumericQ, c6_?NumericQ, c7_?NumericQ] :=
  Norm[pw[[9]] + c0 pw[[1]] + c1 pw[[2]] + c2 pw[[3]] + c3 pw[[4]] +
       c4 pw[[5]] + c5 pw[[6]] + c6 pw[[7]] + c7 pw[[8]]];

vrun["d5.machine",
  First[FindMinimum[obj[a0, a1, a2, a3, a4, a5, a6, a7],
    {{a0, 0.}, {a1, 0.}, {a2, 0.}, {a3, 0.}, {a4, 0.}, {a5, 0.}, {a6, 0.}, {a7, 0.}}]], 10];
vrun["d5.wp50",
  First[FindMinimum[obj[a0, a1, a2, a3, a4, a5, a6, a7],
    {{a0, 0}, {a1, 0}, {a2, 0}, {a3, 0}, {a4, 0}, {a5, 0}, {a6, 0}, {a7, 0}},
    WorkingPrecision -> 60]], 48, 1, 600];
