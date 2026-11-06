(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 148 -- ISSAC '97 holdout B, family 3: slowly convergent series
   (the counterpart of Problem 3, Sum[1/(n^Pi + n^2 + n^Sqrt[2] + 1)^(1/3)])
   ==========================================================================

       cd benchmarks/148-issac97b-slow-series
       ../../Mathilda -file slow_series.m

   The original's terms are smooth, positive and decay like n^-1.047.  These
   ten are slow in other ways: logarithmic decay (alternating, and positive
   with a 1/log N tail), conditional convergence of a non-alternating
   oscillatory series, log-weighted and harmonic-number terms, an alternating
   1/sqrt(k) series, arithmetic (non-smooth) terms DivisorSigma and EulerPhi,
   and an oscillatory Bessel series whose terms decay only like n^-1/2.
   True values (closed form where one exists), each established twice
   without either benchmarked system:

   b3.1  sum_{n>=1} (-1)^(n+1)/log(n+1) = 0.924299897222938855959570
         (a) Cohen-Villegas-Zagier acceleration; (b) mpmath nsum (Euler/Levin).
   b3.2  sum_{n>=2} 1/(n log^2 n) = 2.109742801236891974479257  (tail ~ 1/log N)
         (a) Euler-Maclaurin by hand from N = 2000 (B_2..B_14 corrections);
         (b) the published value 2.10974280123689197447 (OEIS A115563).
   b3.3  sum sin(n)/n = (pi - 1)/2 = 1.070796326794896619231322
         (a) the sawtooth Fourier series (pi - x)/2 at x = 1;
         (b) Im(-log(1 - e^i)); (c) Riesz means of order 4 to N = 10^6.
   b3.4  sum log(n)/n^2 = -zeta'(2) = 0.937548254315843753702574
         (a) mpmath zeta derivative; (b) Euler-Maclaurin by hand.
   b3.5  sum H_n/n^2 = 2 zeta(3) = 2.404113806319188570799476  (Euler)
         (a) the identity; (b) Euler-Maclaurin on (psi(x+1) + gamma)/x^2.
   b3.6  sum (-1)^(k+1)/sqrt(k) = (1 - sqrt 2) zeta(1/2)
         = 0.604898643421630370247266   (a) eta(1/2); (b) CVZ acceleration.
   b3.7  sum_{n>=1} (-1)^n log(n)/n = eta'(1) = gamma log 2 - (log 2)^2/2
         = 0.159868903742430971756948   (a) the identity; (b) CVZ.
   b3.8  sum d(n)/n^2 = zeta(2)^2 = pi^4/36 = 2.705808084277845478790009
         (a) Dirichlet-series identity; (b) a sieve to 10^6 plus the tail
         (log N + 1 + 2 gamma)/N: 2.7058080844, 9 digits.
   b3.9  sum phi(n)/n^3 = zeta(2)/zeta(3) = 1.368432777620205875736766
         (a) identity; (b) sieve to 10^6 plus tail 6/(pi^2 N): 1.3684327776202.
   b3.10 sum_{n>=1} J_0(n) = 1/2 (Poisson summation: sum_{n in Z} J_0(n x) =
         2/x for 0 < x < 2 pi, the Fourier transform of J_0 being supported on
         [-1, 1]); Riesz means of order 4: 0.50000342 (N = 10^5), 0.500000342
         (N = 10^6), converging as 1/N.

   METHODS.  NSum with default options at machine precision and with
   WorkingPrecision -> 30.
   CHECKS.  Machine 5-7 digits, 12 at 30 digits; one fewer where the true value
   sits within 0.05 units of a rounding edge (b3.3 at 30 digits, b3.4 machine).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NSum", "HarmonicNumber", "DivisorSigma", "EulerPhi", "BesselJ"}];

mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];

run["b3.1 NSum machine", NSum[(-1)^(n + 1)/Log[n + 1], {n, 1, Infinity}], mc[7], 3];
run["b3.1 NSum wp30", NSum[(-1)^(n + 1)/Log[n + 1], {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.2 NSum machine", NSum[1/(n Log[n]^2), {n, 2, Infinity}], mc[6], 3];
run["b3.2 NSum wp30", NSum[1/(n Log[n]^2), {n, 2, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.3 NSum machine", NSum[Sin[n]/n, {n, 1, Infinity}], mc[6], 3];
run["b3.3 NSum wp30", NSum[Sin[n]/n, {n, 1, Infinity}, WorkingPrecision -> 30], mc[11]];

run["b3.4 NSum machine", NSum[Log[n]/n^2, {n, 1, Infinity}], mc[6], 3];
run["b3.4 NSum wp30", NSum[Log[n]/n^2, {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.5 NSum machine", NSum[HarmonicNumber[n]/n^2, {n, 1, Infinity}], mc[7], 3];
run["b3.5 NSum wp30", NSum[HarmonicNumber[n]/n^2, {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.6 NSum machine", NSum[(-1)^(k + 1)/Sqrt[k], {k, 1, Infinity}], mc[7], 3];
run["b3.6 NSum wp30", NSum[(-1)^(k + 1)/Sqrt[k], {k, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.7 NSum machine", NSum[(-1)^n Log[n]/n, {n, 1, Infinity}], mc[7], 3];
run["b3.7 NSum wp30", NSum[(-1)^n Log[n]/n, {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.8 NSum machine", NSum[DivisorSigma[0, n]/n^2, {n, 1, Infinity}], mc[6], 3];
run["b3.8 NSum wp30", NSum[DivisorSigma[0, n]/n^2, {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.9 NSum machine", NSum[EulerPhi[n]/n^3, {n, 1, Infinity}], mc[6], 3];
run["b3.9 NSum wp30", NSum[EulerPhi[n]/n^3, {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b3.10 NSum machine", NSum[BesselJ[0, n], {n, 1, Infinity}], mc[5], 3];
run["b3.10 NSum wp30", NSum[BesselJ[0, n], {n, 1, Infinity}, WorkingPrecision -> 30], mc[12]];
