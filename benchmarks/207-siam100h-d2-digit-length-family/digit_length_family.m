(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 207 -- SIAM 100-digit HOLDOUT A, family d2 (self-similar
   digit-length series; original: experiment 177, Sum 1/f(n) with
   f(n) = n f(#binary digits of n)).

   L_b(n) = number of base-b digits of n.  Ten new series:
     hd2.1   Sum 1/f(n),        f(n) = n^2 f(L_2(n)),  f(1)=1, f(2)=2
     hd2.2   Sum 1/f(n),        f(n) = n f(L_2(n)),    f(1)=1, f(2)=3
     hd2.3   Sum (-1)^(n+1)/f(n), f the original (f(n) = n f(L_2 n), f(1)=1, f(2)=2)
     hd2.4   Sum 1/g(n),        g(n) = n L g(L), L = L_10(n);  g(n) = n for n <= 9
     hd2.5   Sum 1/f(n),        f(n) = n f(L_3(n))^2,  f(1)=1, f(2)=2
     hd2.6   Sum 1/f(n)^2,      f the original
     hd2.7   Sum 1/f(n),        f(n) = n f(L_2(n) - 1) = n f(Floor[Log2 n]),  f(1)=1
     hd2.8   Sum (-1)^(n+1)/g(n), g as in hd2.4
     hd2.9   Sum 1/f(n),        f(n) = n Sqrt[L] f(L), L = L_2(n);  f(1)=1, f(2)=2
     hd2.10  Sum over ODD n of 1/f(n), f the original (stated as Sum_k 1/f(2k-1))
   Generalised along the original's axis: the term is a smooth n^-a times a
   factor constant on each digit-length block, defined recursively through
   the block index, so the series is self-similar.  The non-alternating ones
   with a linear n converge like 1/log n, 1/log^2 n, ... (hd2.2, hd2.4, hd2.5,
   hd2.7, hd2.9, hd2.10): no extrapolation that sees only the terms can
   reach 10 digits.  hd2.1/hd2.6 (squared terms) converge like 1/N and the
   alternating ones (hd2.3, hd2.8) like their first omitted term -- the
   "easier" members; but the term is only piecewise smooth (it jumps at every
   power of b), so an Euler-Maclaurin / Wynn extrapolation still stalls at a
   few digits (verification run: 3-8 digits in every system).

   Idiomatic statement only: NSum of the direct term, machine precision and
   WorkingPrecision -> 60; f as written (glue).  References: exact block sums
   (digamma / Hurwitz zeta) with the self-similar linear equation, see
   siam100_refs/hd2_ref.py and .wl; truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "BitLength", "IntegerLength"}];

f0[n_Integer] := If[n <= 2, n, n f0[BitLength[n]]];
f1[n_Integer] := If[n <= 2, n, n^2 f1[BitLength[n]]];
f2[n_Integer] := Which[n == 1, 1, n == 2, 3, True, n f2[BitLength[n]]];
g4[n_Integer] := If[n <= 9, n, n IntegerLength[n] g4[IntegerLength[n]]];
f5[n_Integer] := If[n <= 2, n, n f5[IntegerLength[n, 3]]^2];
f7[n_Integer] := If[n == 1, 1, n f7[BitLength[n] - 1]];
f9[n_Integer] := If[n <= 2, n, n Sqrt[BitLength[n]] f9[BitLength[n]]];

vrun["hd2.1 machine", NSum[1/f1[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.1 wp50", NSum[1/f1[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.2 machine", NSum[1/f2[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.2 wp50", NSum[1/f2[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.3 machine", NSum[(-1)^(n + 1)/f0[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.3 wp50", NSum[(-1)^(n + 1)/f0[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.4 machine", NSum[1/g4[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.4 wp50", NSum[1/g4[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.5 machine", NSum[1/f5[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.5 wp50", NSum[1/f5[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.6 machine", NSum[1/f0[n]^2, {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.6 wp50", NSum[1/f0[n]^2, {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.7 machine", NSum[1/f7[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.7 wp50", NSum[1/f7[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.8 machine", NSum[(-1)^(n + 1)/g4[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.8 wp50", NSum[(-1)^(n + 1)/g4[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.9 machine", NSum[1/f9[n], {n, 1, Infinity}], 10, 1, 60];
vrun["hd2.9 wp50", NSum[1/f9[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["hd2.10 machine", NSum[1/f0[2 k - 1], {k, 1, Infinity}], 10, 1, 60];
vrun["hd2.10 wp50", NSum[1/f0[2 k - 1], {k, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
