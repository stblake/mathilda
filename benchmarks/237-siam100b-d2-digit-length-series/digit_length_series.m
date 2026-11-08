(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 237 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   problem 2 (Sum 1/f(n), f defined through the digit length of n): ten
   other series whose terms depend on the NUMBER OF DIGITS of n, or that decay
   only logarithmically.  d(n) = IntegerLength[n], b(n) = BitLength[n],
   t(n) = IntegerLength[n, 3].

     bd2.1   Sum_{n>=1} 1/(n d(n)^2)                         = 4.325970846...
     bd2.2   Sum_{n>=1} 1/(n b(n)^2)                         = 1.492454526...
     bd2.3   Sum_{n>=2} 1/(n Log[n]^2)          (smooth)     = 2.109742801...
     bd2.4   Sum_{n>=3} 1/(n Log[n] Log[Log[n]]^2) (smooth)  = 38.40676809...
     bd2.5   Sum 1/g(n), g(1)=1, g(n) = n g(b(n)-1)          = 3.869696587...
     bd2.6   Sum 1/f(n), f(1)=1, f(n) = n f(d(n))^2          = 4.097770429...
     bd2.7   Sum_{n>=1} 1/(n d(n) d(d(n))^2)                 = 10.51183632...
     bd2.8   Sum_{n>=1} 1/(n d(n)!)                          = 4.506555879...
     bd2.9   Sum_{n>=1} (-1)^b(n)/(n b(n))                   = -0.7335342757...
     bd2.10  Sum 1/f(n), f(1)=1, f(2)=2, f(n) = 2n f(t(n))   = 2.180233390...

   All converge, most like 1/log N or slower (bd2.4: like 1/loglog N, so the
   partial sum to n = 10^100 is still wrong in the first digit).  The digit-
   length rows are constant on blocks b^(k-1) <= n < b^k, so no
   acceleration that sees only the terms can get them; grouping by digit
   length (psi differences) does, and is how truth.tsv was computed.  bd2.3
   and bd2.4 are smooth: Euler-Maclaurin with an accurate tail integral gets
   them in milliseconds, which is what a correct NSum should do.

   Idiomatic statement only: NSum of the direct term, at machine precision
   and with WorkingPrecision 60.  The recursive f/g are the definitions as
   written (glue); digit-length summands are ?NumericQ helpers on Round[n].
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "IntegerLength", "BitLength"}];

g5[1] = 1; g5[n_Integer] := n g5[BitLength[n] - 1];
f6[1] = 1; f6[n_Integer] := n f6[IntegerLength[n]]^2;
f10[n_Integer] := If[n <= 2, n, 2 n f10[IntegerLength[n, 3]]];
tm[1][m_] := 1/(m IntegerLength[m]^2);
tm[2][m_] := 1/(m BitLength[m]^2);
tm[5][m_] := 1/g5[m];
tm[6][m_] := 1/f6[m];
tm[7][m_] := 1/(m IntegerLength[m] IntegerLength[IntegerLength[m]]^2);
tm[8][m_] := 1/(m IntegerLength[m]!);
tm[9][m_] := (-1)^BitLength[m]/(m BitLength[m]);
tm[10][m_] := 1/f10[m];
t[p_][n_?NumericQ] := tm[p][Round[n]];
(* NHoldAll: NSum applies N to the summand, which would turn t[1][n] into
   t[1.][n] -- a head with no tm definition, so the sum came back unevaluated. *)
SetAttributes[t, NHoldAll];

Do[
  vrun["bd2." <> ToString[p] <> " machine", NSum[t[p][n], {n, 1, Infinity}], 10, 1, 60];
  vrun["bd2." <> ToString[p] <> " wp50",
       NSum[t[p][n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60],
  {p, {1, 2}}];
vrun["bd2.3 machine", NSum[1/(n Log[n]^2), {n, 2, Infinity}], 10, 1, 60];
vrun["bd2.3 wp50", NSum[1/(n Log[n]^2), {n, 2, Infinity}, WorkingPrecision -> 60], 48, 1, 60];
vrun["bd2.4 machine", NSum[1/(n Log[n] Log[Log[n]]^2), {n, 3, Infinity}], 10, 1, 60];
vrun["bd2.4 wp50", NSum[1/(n Log[n] Log[Log[n]]^2), {n, 3, Infinity},
                        WorkingPrecision -> 60], 48, 1, 60];
Do[
  vrun["bd2." <> ToString[p] <> " machine", NSum[t[p][n], {n, 1, Infinity}], 10, 1, 60];
  vrun["bd2." <> ToString[p] <> " wp50",
       NSum[t[p][n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60],
  {p, {5, 6, 7, 8, 9, 10}}];
