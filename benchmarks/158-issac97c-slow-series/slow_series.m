(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 158 -- ISSAC '97 holdout C (SEALED), family c3:
   slowly convergent series (counterpart of Problem 3)
   ==========================================================================

       cd benchmarks/158-issac97c-slow-series
       ../../Mathilda -file slow_series.m

   Ten series whose terms are arithmetic (primes, Liouville), nested
   (harmonic numbers, an alternating Euler sum), q-series (a Lambert series at
   q = 0.999, needing ~10^5 terms), logarithmically convergent (a tail of
   1/Log[N]), oscillating with an irrational phase (Cos[n]/Sqrt[n]), or barely
   convergent (Zeta at 1.01, whose tail is 100 N^(-0.01)).  Where the terms are
   analytic in n the rows are NSum (machine) and NSum with WorkingPrecision ->
   40 (wp40).  For the arithmetic series NSum has nothing to extrapolate, so
   the machine row is a long direct Sum whose truncation error is below the
   digits checked, and the wp row asks for the identity that makes the series
   computable (the symbolic Dirichlet series, the Moebius log-zeta formula,
   the polylogarithm).

   TRUE ANSWERS (A = closed form or mpmath with its own acceleration, 50
   digits; B = PARI/GP at 60 digits: suminf/sumalt/zeta/zetahurwitz, or a
   separate route as stated).  A and B agree to >= 30 digits.  Digits: machine
   6-8 (5 for the 1e5-prime partial sum, whose tail is 5e-8), wp 20-25 (15 for
   the two logarithmically convergent sums).

   c3.1  sum_p 1/p^2 = PrimeZeta(2) = 0.452247420041065498506543364832.
         A: mpmath primezeta; B: PARI suminf of MoebiusMu[k]/k Log[Zeta[2k]].
         machine: the first 10^5 primes (0.4522473688..., 5 digits); wp50: the
         Moebius log-zeta series to k = 80 (terms ~ 4^-k/k).
   c3.2  sum LiouvilleLambda[n]/n^2 = Zeta[4]/Zeta[2] = Pi^2/15
         = 0.657973626739290574588966066658.  B: PARI partial sums to 10^6
         (0.6579736270731, tail ~ N^-1.5).  machine: Sum to 10^5 (7 digits);
         wp30: N of the symbolic Sum to Infinity.
   c3.3  sum (-1)^(n+1) HarmonicNumber[n]/n = Pi^2/12 - Log[2]^2/2
         = 0.582240526465012505902656320160.  B: mpmath nsum.
   c3.4  sum q^n/(1-q^n), q = 999/1000 (Lambert; = sum DivisorSigma[0,n] q^n)
         = 7480.97786909167252764630707088.  A: 2e5-term direct sum at 50
         digits; B: PARI suminf.
   c3.5  sum_{n>=2} 1/(n Log[n]^2) = 2.10974280123689197447925719762.
         A: Euler-Maclaurin at N = 1000 with the exact tail 1/Log[N] and 11
         Bernoulli corrections; B: the published value (OEIS A115563).
         (mpmath's own Euler-Maclaurin nsum returns 1.9513: wrong, not used.)
   c3.6  sum_{n>=2} (-1)^n/Log[n] = 0.924299897222938855959570181360.
         A: mpmath nsum; B: PARI sumalt.
   c3.7  sum Cos[n]/Sqrt[n] = Re PolyLog[1/2, E^I]
         = -0.194108935092182049739149244928.  A: mpmath polylog; B: PARI,
         Jonquiere's formula through zetahurwitz at x = 1/(2 Pi).
         machine: NSum; wp40: the polylogarithm.
   c3.8  sum HarmonicNumber[n]^2/n^2 = 17 Pi^4/360
         = 4.59987374327233731394301571030.  B: partial sum to 2e4 plus the
         asymptotic tail (agrees to 1e-12; mpmath nsum fails here).
   c3.9  sum n^(-101/100) = Zeta[101/100] = 100.577943338496872490282154286.
         A: mpmath zeta; B: PARI zeta.
   c3.10 sum Log[n]/n^2 = -Zeta'[2] = 0.937548254315843753702574094568.
         A: mpmath zeta derivative; B: PARI zeta'.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NSum", "Prime", "MoebiusMu", "LiouvilleLambda", "HarmonicNumber", "PolyLog", "Zeta"}];

run["c3.1 Sum primes machine", Sum[1./Prime[k]^2, {k, 1, 10^5}], sig[#, 5] &];
run["c3.1 MoebiusMu log-zeta wp50",
    Sum[MoebiusMu[k]/k Log[Zeta[N[2 k, 50]]], {k, 1, 80}], sig[#, 25] &];

run["c3.2 Sum machine", Sum[LiouvilleLambda[n]/n^2., {n, 1, 10^5}], sig[#, 7] &];
run["c3.2 Sum Infinity wp30", N[Sum[LiouvilleLambda[n]/n^2, {n, 1, Infinity}], 30], sig[#, 20] &];

run["c3.3 NSum machine", NSum[(-1)^(n + 1) HarmonicNumber[n]/n, {n, 1, Infinity}], sig[#, 8] &, 3];
run["c3.3 NSum wp40", NSum[(-1)^(n + 1) HarmonicNumber[n]/n, {n, 1, Infinity}, WorkingPrecision -> 40],
    sig[#, 20] &];

run["c3.4 NSum machine", NSum[(999/1000)^n/(1 - (999/1000)^n), {n, 1, Infinity}], sig[#, 8] &, 3];
run["c3.4 NSum wp40", NSum[(999/1000)^n/(1 - (999/1000)^n), {n, 1, Infinity}, WorkingPrecision -> 40],
    sig[#, 20] &];

run["c3.5 NSum machine", NSum[1/(n Log[n]^2), {n, 2, Infinity}], sig[#, 6] &, 3];
run["c3.5 NSum wp40", NSum[1/(n Log[n]^2), {n, 2, Infinity}, WorkingPrecision -> 40], sig[#, 15] &];

run["c3.6 NSum machine", NSum[(-1)^n/Log[n], {n, 2, Infinity}], sig[#, 8] &, 3];
run["c3.6 NSum wp40", NSum[(-1)^n/Log[n], {n, 2, Infinity}, WorkingPrecision -> 40], sig[#, 20] &];

run["c3.7 NSum machine", NSum[Cos[n]/Sqrt[n], {n, 1, Infinity}], sig[#, 6] &, 3];
run["c3.7 PolyLog wp40", N[Re[PolyLog[1/2, Exp[I]]], 40], sig[#, 25] &];

run["c3.8 NSum machine", NSum[HarmonicNumber[n]^2/n^2, {n, 1, Infinity}], sig[#, 6] &, 3];
run["c3.8 NSum wp40", NSum[HarmonicNumber[n]^2/n^2, {n, 1, Infinity}, WorkingPrecision -> 40],
    sig[#, 15] &];

run["c3.9 NSum machine", NSum[n^(-101/100), {n, 1, Infinity}], sig[#, 8] &, 3];
run["c3.9 NSum wp40", NSum[n^(-101/100), {n, 1, Infinity}, WorkingPrecision -> 40], sig[#, 20] &];

run["c3.10 NSum machine", NSum[Log[n]/n^2, {n, 1, Infinity}], sig[#, 8] &, 3];
run["c3.10 NSum wp40", NSum[Log[n]/n^2, {n, 1, Infinity}, WorkingPrecision -> 40], sig[#, 20] &];
