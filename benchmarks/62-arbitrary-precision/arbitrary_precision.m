(* ==========================================================================
   Experiment 62 -- Arbitrary-precision numerics (group E)
   ==========================================================================
   WHAT IT MEASURES.  The MPFR-backed high-precision path: N[expr, p] and the
   WorkingPrecision option honoured by NSum, NIntegrate and FindRoot
   (src/precision.c and the *_mpfr.c modules).  This is the one group-E
   experiment whose baseline is NOT scipy -- scipy is double-only -- but mpmath,
   the standard Python arbitrary-precision library.  A README backlog item.

   Single-run (benchOnce): N[Pi, 10000] and friends are expensive on their own,
   and a warm-up would let Mathilda's N[] result cache serve every later rep for
   free -- measuring the cache, not the arithmetic.

   Checks use a 6-digit window pulled from DEEP in the expansion
   (win[v, s] = digits s-5..s of v), so the check verifies real precision, not
   just the leading figures.  FindRoot at WorkingPrecision -> 200 now delivers
   ~wp/2 = 100 correct digits (win50); it previously stopped at ~19 digits while
   labelling the result Precision 100, a correctness bug fixed by scaling
   AccuracyGoal with WorkingPrecision (findroot.c).  NSum at WorkingPrecision ->
   100 now reaches ~90+ digits (win60); it previously floored at ~44 because its
   Euler–Maclaurin base did not scale with precision (nsum.c).
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"N", "NSum", "NIntegrate", "FindRoot", "WorkingPrecision"}];

win[v_, start_] := Mod[IntegerPart[v 10^start], 10^6];

benchOnce["N[Pi, 10000]", N[Pi, 10000];];
check["N[Pi, 10000]", win[N[Pi, 10000], 5000]];

benchOnce["N[Gamma[1/3], 1000]", N[Gamma[1/3], 1000];];
check["N[Gamma[1/3], 1000]", win[N[Gamma[1/3], 1000], 500]];

benchOnce["N[Log[2], 5000]", N[Log[2], 5000];];
check["N[Log[2], 5000]", win[N[Log[2], 5000], 2000]];

benchOnce["NSum 1/k^2 WorkingPrecision 100",
          NSum[1/k^2, {k, 1, Infinity}, WorkingPrecision -> 100];];
check["NSum 1/k^2 WorkingPrecision 100",
      win[NSum[1/k^2, {k, 1, Infinity}, WorkingPrecision -> 100], 60]];

benchOnce["NIntegrate Exp[-x^2] WorkingPrecision 200",
          NIntegrate[Exp[-x^2], {x, 0, 1}, WorkingPrecision -> 200];];
check["NIntegrate Exp[-x^2] WorkingPrecision 200",
      win[NIntegrate[Exp[-x^2], {x, 0, 1}, WorkingPrecision -> 200], 50]];

benchOnce["FindRoot Cos[x]=x WorkingPrecision 200",
          FindRoot[Cos[x] == x, {x, 1}, WorkingPrecision -> 200];];
check["FindRoot Cos[x]=x WorkingPrecision 200",
      win[x /. FindRoot[Cos[x] == x, {x, 1}, WorkingPrecision -> 200], 50]];
