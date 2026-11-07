(* bench-timeout: 420 *)
(* ==========================================================================
   Experiment 183 -- SIAM 100-Digit Challenge, Appendix D Problem 8
   ("A Sobolev constant"; problem 10 on the book's web page).

     c_m(R^N)^2 = (2 Pi)^-N Int_{R^N} ( Sum_{|k| <= m} x^(2k) )^-1 dx,
     k a multi-index (k_1..k_N) >= 0, |k| = Sum k_j, x^k = Prod x_j^k_j.
     E.g. c_3(R^1) = 0.5.  What is c_10(R^10)?     Answer: 5.5489 e-4
     (only ~5 digits are established -- see truth.tsv; no published value).

   Idiomatic statement: the 10-dimensional NIntegrate over R^10, as written.
   Sum_{|k|<=10} y^k (y_j = x_j^2) has C(20,10) = 184756 monomials, so the
   integrand evaluates it as Sum_{j<=10} h_j(y) with the complete homogeneous
   polynomials h_j from the power sums (Newton's identity) -- the same
   polynomial, evaluated in O(100) operations.  No dimension reduction is
   applied: the spherical / radial reduction that makes the problem tractable
   is the reference method (siam100_refs/d8_ref.py), not the case.
     d8.machine        the 10-D integral at machine precision;
     d8.wp50           the same with WorkingPrecision -> 60 (no system has
                       a 50-digit 10-D cubature: expected timeout / SKIP);
     d8.machine c3R1   sanity row from the problem statement, c_3(R^1) = 1/2
                       (1-D NIntegrate over R);
     d8.machine c2R2   reduced instance c_2(R^2) = 0.31875906098...
                       (2-D NIntegrate over R^2), where 10 digits are known.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

(* Sum_{|k|<=10} y^k = Sum_{j=0}^{10} h_j(y), h_j by Newton's identity *)
hsum[y_List] := Module[{p = Table[Total[y^j], {j, 1, 10}], h = Table[0, {11}]},
  h[[1]] = 1;
  Do[h[[j + 1]] = Sum[p[[i]] h[[j - i + 1]], {i, 1, j}]/j, {j, 1, 10}];
  Total[h]];
fS[x1_?NumericQ, x2_?NumericQ, x3_?NumericQ, x4_?NumericQ, x5_?NumericQ,
   x6_?NumericQ, x7_?NumericQ, x8_?NumericQ, x9_?NumericQ, x10_?NumericQ] :=
  1/hsum[{x1, x2, x3, x4, x5, x6, x7, x8, x9, x10}^2];

vrun["d8.machine", Sqrt[NIntegrate[fS[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10],
    {x1, -Infinity, Infinity}, {x2, -Infinity, Infinity}, {x3, -Infinity, Infinity},
    {x4, -Infinity, Infinity}, {x5, -Infinity, Infinity}, {x6, -Infinity, Infinity},
    {x7, -Infinity, Infinity}, {x8, -Infinity, Infinity}, {x9, -Infinity, Infinity},
    {x10, -Infinity, Infinity}]/(2 Pi)^10], 3];
vrun["d8.wp50", Sqrt[NIntegrate[fS[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10],
    {x1, -Infinity, Infinity}, {x2, -Infinity, Infinity}, {x3, -Infinity, Infinity},
    {x4, -Infinity, Infinity}, {x5, -Infinity, Infinity}, {x6, -Infinity, Infinity},
    {x7, -Infinity, Infinity}, {x8, -Infinity, Infinity}, {x9, -Infinity, Infinity},
    {x10, -Infinity, Infinity}, WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^10], 3];
vrun["d8.machine c3R1", Sqrt[NIntegrate[1/(1 + x^2 + x^4 + x^6), {x, -Infinity, Infinity}]/(2 Pi)], 10, 3];
vrun["d8.machine c2R2", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + x^4 + x^2 y^2 + y^4),
    {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2], 10, 3];
