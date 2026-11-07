(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 206 -- SIAM 100-digit HOLDOUT A, family d1 (Kempner-type
   series; original: experiment 176, Sum 1/n over n avoiding "42").

   Ten new digit-restricted harmonic series, each Sum f(n) over the positive
   integers n whose representation passes a digit test:
     hd1.1   1/n,   decimal digits contain no 9          (Kempner 1914)
     hd1.2   1/n,   no digit 0
     hd1.3   1/n,   no substring "00"
     hd1.4   1/n,   no substring "314"                    (sum ~ 2300)
     hd1.5   1/n,   binary digits contain no "11"         (Fibbinary n)
     hd1.6   1/n,   ternary digits contain no 2           (Cantor-set n)
     hd1.7   1/n^2, no digit 9
     hd1.8   1/n,   every decimal digit odd
     hd1.9   1/n,   decimal digits non-decreasing (left to right)
     hd1.10  1/n,   no substring "2026"                   (sum ~ 23021)
   Generalised along the original's axis -- the summand is an indicator on
   the digit string times a smooth term, so the series is not the sampling of
   any smooth function and its tail thins out only geometrically per DIGIT
   (ratio (b-1)/b, 9.9/10, ..., per decade).  Varied in kind: forbidden
   digit / forbidden substring of length 2, 3, 4 / allowed digit set /
   monotone digit strings / bases 2, 3, 10 / exponent 1 and 2.

   Idiomatic statement only: NSum of the indicator-weighted term to Infinity
   (a ?NumericQ helper, as in experiment 176), at machine precision and at
   WorkingPrecision 60.  No built-in knows the Schmelzer-Baillie power-sum
   recursion that sums these (the references use it: siam100_refs/hd1_ref.py and .wl).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "IntegerString", "StringFreeQ", "IntegerDigits", "OrderedQ"}];

free[n_, s_String, b_: 10] := StringFreeQ[IntegerString[Round[n], b], s];
t1[n_?NumericQ] := Boole[free[n, "9"]]/n;
t2[n_?NumericQ] := Boole[free[n, "0"]]/n;
t3[n_?NumericQ] := Boole[free[n, "00"]]/n;
t4[n_?NumericQ] := Boole[free[n, "314"]]/n;
t5[n_?NumericQ] := Boole[free[n, "11", 2]]/n;
t6[n_?NumericQ] := Boole[free[n, "2", 3]]/n;
t7[n_?NumericQ] := Boole[free[n, "9"]]/n^2;
t8[n_?NumericQ] := Boole[And @@ OddQ[IntegerDigits[Round[n]]]]/n;
t9[n_?NumericQ] := Boole[OrderedQ[IntegerDigits[Round[n]]]]/n;
t10[n_?NumericQ] := Boole[free[n, "2026"]]/n;

vrun["hd1.1 machine", NSum[t1[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.1 wp50", NSum[t1[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.2 machine", NSum[t2[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.2 wp50", NSum[t2[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.3 machine", NSum[t3[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.3 wp50", NSum[t3[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.4 machine", NSum[t4[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.4 wp50", NSum[t4[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.5 machine", NSum[t5[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.5 wp50", NSum[t5[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.6 machine", NSum[t6[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.6 wp50", NSum[t6[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.7 machine", NSum[t7[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.7 wp50", NSum[t7[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.8 machine", NSum[t8[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.8 wp50", NSum[t8[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.9 machine", NSum[t9[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.9 wp50", NSum[t9[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd1.10 machine", NSum[t10[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd1.10 wp50", NSum[t10[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
