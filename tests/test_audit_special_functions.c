/* Builtin-audit regressions: special functions (review/issues/<Head>.md;
 * tests/test_audit.h).  Every probe runs under TimeConstrained so a
 * regression to the old hangs fails (as $Aborted) instead of wedging. */
#include "test_audit.h"

/* BetaRegularized#1 (high): BetaRegularized hangs for a large second shape
 * parameter (b ~ 1e9 to 1e15).  The machine kernel declined a + b >= 1e7 and
 * arb's 2F1(a, 1 - b; a + 1; z) with |1 - b| ~ 1e10 never returned.  For
 * a = 2 the closed form I_x(2, b) = 1 - (1 - x)^b (1 + b x) is the
 * independent check. */
static void test_betaregularized_large_b(void) {
    const char* id = "BetaRegularized#1";
    audit_true(id, "TimeConstrained[BetaRegularized[0.5, 2, 1.*^10], 10] === 1.");
    audit_true(id, "TimeConstrained[BetaRegularized[0.5, 2, 1.*^9], 10] === 1.");
    audit_true(id, "TimeConstrained[BetaRegularized[0.5, 2, 1.*^15], 10] === 1.");
    audit_true(id, "TimeConstrained[BetaRegularized[0.5, 2, 10^10], 10] === 1.");
    audit_true(id, "TimeConstrained[BetaRegularized[0.5, 0.3, 1.*^9], 10] === 1.");
    /* the closed form, evaluated in log space at 40 digits */
    audit_run("audCl[x_, b_] := Module[{xx = N[x, 40], bb = N[b, 40]}, "
              "1 - Exp[bb Log[1 - xx]] (1 + bb xx)]");
    audit_true(id, "Abs[TimeConstrained[BetaRegularized[1.*^-5, 2, 1.*^10], 10] - "
                   "audCl[1/100000, 10^10]] < 1.*^-15");
    /* below the cutoff the complement is not negligible: still the old route */
    audit_true(id, "Abs[TimeConstrained[BetaRegularized[1.*^-6, 2, 1.*^7], 10] - "
                   "audCl[1/1000000, 10^7]] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[BetaRegularized[1.*^-9, 2, 1.*^10], 10] - "
                   "audCl[1/1000000000, 10^10]] < 1.*^-15");
    /* the machine surfaces agree */
    audit_true(id, "With[{bs = {1.*^8, 1.*^10, 1.*^15}}, "
                   "(BetaRegularized[0.5, 2, #] & /@ bs) === BetaRegularized[0.5, 2, bs] === "
                   "Normal[BetaRegularized[0.5, 2, NDArray[bs]]]]");
    /* sanity: ordinary parameters unchanged */
    audit_true(id, "Abs[BetaRegularized[0.3, 2., 3.] - N[BetaRegularized[3/10, 2, 3], 20]] < 1.*^-15");
    audit_true(id, "BetaRegularized[0.5, 2, 1.*^7] === 1.");
}

/* ExpIntegralE#1 (high): ExpIntegralE[10^15, 1.] hangs (huge integer order
 * with machine argument).  The machine kernel declined n > 1e6 and arb's cost
 * grows with n.  Independent check: the continued fraction
 * E_n(x) = e^-x/(x+n - n/(x+n+2 - 2(n+1)/(x+n+4 - ...))) in exact rationals,
 * depth 12 (converges fast for large n). */
static void test_expintegrale_large_order(void) {
    const char* id = "ExpIntegralE#1";
    audit_run("audCf[n_, x_] := Module[{t = x + n + 24}, "
              "Do[t = (x + n + 2 (i - 1)) - i (n - 1 + i)/t, {i, 12, 1, -1}]; N[Exp[-x]/t, 30]]");
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[10^15, 1.], 10]/audCf[10^15, 1] - 1] < 4.*^-16");
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[10^18, 1.], 10]/audCf[10^18, 1] - 1] < 4.*^-16");
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[10^6 + 1, 0.001], 10]/"
                   "audCf[10^6 + 1, 1/1000] - 1] < 4.*^-16");
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[3*10^9, 10.], 10]/audCf[3*10^9, 10] - 1] < 4.*^-16");
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[10^7, 1.*^6], 10]/audCf[10^7, 10^6] - 1] < 4.*^-16");
    /* a double underflow is recomputed in MPFR: that route must not hang either */
    audit_true(id, "Abs[TimeConstrained[ExpIntegralE[10^15, 1000.], 10]/"
                   "audCf[10^15, 1000] - 1] < 4.*^-16");
    /* arbitrary precision at a large order */
    audit_true(id, "Abs[TimeConstrained[N[ExpIntegralE[10^15, 1], 30], 10] - audCf[10^15, 1]] < "
                   "10^-44");
    /* the recurrence E_{n+1} = (e^-x - x E_n)/n across the old n = 1e6 cutoff */
    audit_true(id, "Abs[ExpIntegralE[10^6 + 1, 2.]/((Exp[-2.] - 2. ExpIntegralE[10^6, 2.])/10^6) - 1] < 1.*^-14");
    /* scalar, packed, NDArray and Compile agree bit for bit */
    audit_true(id, "With[{xs = {0.5, 1., 10., 300.}}, "
                   "(ExpIntegralE[10^15, #] & /@ xs) === ExpIntegralE[10^15, xs] === "
                   "Normal[ExpIntegralE[10^15, NDArray[xs]]] === "
                   "Compile[{{v, _Real, 1}}, ExpIntegralE[1000000000000000, v]][xs]]");
}

/* StruveH#1 (high): StruveH at large machine z (about 1e7 to 1e17) hangs.
 * The large-z asymptotic of the machine kernel turned NaN (a_k overflowed
 * while z^-k underflowed), so every large z fell to the convergent 1F2 at
 * ~1.44 z bits.  Independent check: H = Y + K with the reference's BesselY at
 * 40 digits and the K asymptotic sum in exact arithmetic; and the arb
 * asymptotic route against the convergent 1F2 at z = 3000, 1100 digits. */
static void test_struveh_large_z(void) {
    const char* id = "StruveH#1";
    audit_run("audSH[nu_, z_] := N[BesselY[nu, z] + (1/Pi) Sum[Gamma[k + 1/2] "
              "(z/2)^(nu - 2 k - 1)/Gamma[nu + 1/2 - k], {k, 0, 8}], 40]");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 1.*^8], 10]/audSH[0, 10^8] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 1.*^8], 10]] <= 8.*^-5");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 1.*^7], 10]/audSH[0, 10^7] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 1.*^12], 10]/audSH[0, 10^12] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 1.*^17], 10]/audSH[0, 10^17] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[-0.3, 1.*^8], 10]/audSH[-3/10, 10^8] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[7, 1.*^5], 10]/audSH[7, 10^5] - 1] < 1.*^-15");
    audit_true(id, "Abs[TimeConstrained[StruveH[0, 3000.], 10]/audSH[0, 3000] - 1] < 1.*^-15");
    /* arbitrary precision at large z (the 1F2 needed 1.4e8 bits here) */
    audit_true(id, "Abs[TimeConstrained[N[StruveH[0, 10^8], 30], 10] - audSH[0, 10^8]] < 10^-34");
    audit_true(id, "Abs[TimeConstrained[N[StruveH[1/2, 10^8], 30], 10] - "
                   "N[Sqrt[2/(Pi 10^8)] (1 - Cos[10^8]), 40]] < 10^-34");
    audit_true(id, "With[{lo = N[StruveH[-3/10, 3000], 30], hi = N[StruveH[-3/10, 3000], 1100]}, "
                   "Abs[SetPrecision[lo, 1100] - hi] < 10^-32 Abs[hi]]");
    /* scalar, packed, NDArray and Compile agree bit for bit */
    audit_true(id, "With[{xs = {1.*^3, 1.*^5, 1.*^8, 1.*^12, 1.*^17}}, "
                   "(StruveH[0, #] & /@ xs) === StruveH[0, xs] === Normal[StruveH[0, NDArray[xs]]] === "
                   "Compile[{{v, _Real, 1}}, StruveH[0, v]][xs]]");
    /* sanity: moderate z unchanged */
    audit_true(id, "Abs[StruveH[0, 1000.] - N[StruveH[0, 1000], 20]] < 1.*^-17");
    audit_true(id, "Abs[StruveH[1., 30.] - N[StruveH[1, 30], 20]] < 1.*^-14");
}

/* MathieuCPrime#4 (medium): MathieuCPrime runs without terminating for
 * Im z ~ 20.  The Fourier sum needs ~sqrt q e^{|Im z|} terms; far from the
 * axis the Bessel-series continuation now answers.  Independent check: a
 * Taylor continuation of the Mathieu equation y'' = -(a - 2q cos 2w) y from
 * a point the Fourier sum handles (Im w = 6.4) to one the far path handles
 * (Im w = 7), and the Wronskian ce se' - se ce' on the imaginary axis. */
static void test_mathieu_far_imaginary(void) {
    const char* id = "MathieuCPrime#4";
    /* the reported calls: finite (beyond $MaxNumber: Overflow[]) and fast */
    audit_true(id, "TimeConstrained[MathieuCPrime[1.5, 2.0, 0.7 + 30 I], 10] === Overflow[]");
    audit_true(id, "NumberQ[TimeConstrained[MathieuCPrime[1.5, 2.0, 20. I], 10]]");
    audit_true(id, "NumberQ[TimeConstrained[MathieuCPrime[1.5, 2.0, 0.7 + 20 I], 10]]");
    audit_true(id, "NumberQ[TimeConstrained[MathieuC[1.5, 2.0, 0.7 + 12 I], 10]]");
    /* Taylor continuation across the Fourier / Bessel switch (Im z ~ 6.6 for q = 2) */
    audit_run("audStep[a_, q_, w1_, y0_, dy0_, d_, nt_] := Module[{f, y, c2 = Cos[2 w1], s2 = Sin[2 w1]}, "
              "f = Table[If[j == 0, a - 2 q c2, -2 q (c2 If[EvenQ[j], (-1)^(j/2) 2^j/j!, 0] - "
              "s2 If[OddQ[j], (-1)^((j - 1)/2) 2^j/j!, 0])], {j, 0, nt}]; "
              "y = Table[0, {nt + 1}]; y[[1]] = y0; y[[2]] = dy0; "
              "Do[y[[n + 3]] = -Sum[f[[j + 1]] y[[n - j + 1]], {j, 0, n}]/((n + 2) (n + 1)), {n, 0, nt - 2}]; "
              "{Sum[y[[n + 1]] d^n, {n, 0, nt}], Sum[n y[[n + 1]] d^(n - 1), {n, 1, nt}]}]");
    audit_run("audCont[f_, fp_, a_, q_, w1_, w2_, ns_] := Module[{w = w1, v, dv, d = (w2 - w1)/ns}, "
              "v = N[f[a, q, w1], 60]; dv = N[fp[a, q, w1], 60]; "
              "Do[{v, dv} = audStep[N[a, 60], N[q, 60], N[w, 60], v, dv, N[d, 60], 45]; w = w + d, {ns}]; "
              "{v, dv}]");
    audit_true(id, "With[{t = audCont[MathieuC, MathieuCPrime, 3/2, 2, 7/10 + 64/10 I, 7/10 + 7 I, 300]}, "
                   "Abs[t[[1]]/N[MathieuC[3/2, 2, 7/10 + 7 I], 40] - 1] < 10^-30 && "
                   "Abs[t[[2]]/N[MathieuCPrime[3/2, 2, 7/10 + 7 I], 40] - 1] < 10^-30]");
    /* an odd function at a characteristic value (integer nu, the tan-form series) */
    audit_run("audB3 = N[MathieuCharacteristicB[3, 2], 100]");
    audit_true(id, "With[{t = audCont[MathieuS, MathieuSPrime, audB3, 2, -1/2 + 64/10 I, -1/2 + 7 I, 300]}, "
                   "Abs[t[[1]]/N[MathieuS[audB3, 2, -1/2 + 7 I], 40] - 1] < 10^-30]");
    /* Wronskian ce se' - se ce' is constant: at 20 I it equals its value at 0 */
    audit_true(id, "With[{w = N[MathieuC[3/2, 2, 20 I] MathieuSPrime[3/2, 2, 20 I] - "
                   "MathieuS[3/2, 2, 20 I] MathieuCPrime[3/2, 2, 20 I], 30], "
                   "w0 = N[MathieuC[3/2, 2, 0] MathieuSPrime[3/2, 2, 0], 30]}, Abs[w/w0 - 1] < 10^-25]");
    /* sanity: near the real axis the Fourier path is unchanged */
    audit_true(id, "Abs[MathieuCPrime[1.5, 2.0, 0.7 + 5 I] - N[MathieuCPrime[3/2, 2, 7/10 + 5 I], 30]] < "
                   "1.*^-12 Abs[N[MathieuCPrime[3/2, 2, 7/10 + 5 I], 30]]");
}

/* BernoulliB#1 (medium): BernoulliB[5000] took longer than 20 s (O(n^2)
 * big-rational recurrence over every smaller index, odd zeros included).
 * Independent checks: B_{2n} = (-1)^(n+1) 2 (2n)! zeta(2n)/(2 pi)^(2n), the
 * von Staudt-Clausen denominator, and the polynomial route's own recurrence. */
static void test_bernoullib_large_index(void) {
    const char* id = "BernoulliB#1";
    audit_true(id, "Head[TimeConstrained[BernoulliB[5000], 10]] === Rational");
    audit_true(id, "Abs[N[BernoulliB[5000]/(-2 5000! Zeta[5000]/(2 Pi)^5000), 40] - 1] < 10^-35");
    audit_true(id, "Denominator[BernoulliB[5000]] == Times @@ Select[Divisors[5000] + 1, PrimeQ]");
    audit_true(id, "Abs[N[BernoulliB[2002]/(2 2002! Zeta[2002]/(2 Pi)^2002), 40] - 1] < 10^-35");
    /* the direct number agrees with the polynomial route's recurrence */
    audit_true(id, "Union[Table[BernoulliB[n] === BernoulliB[n, 0], {n, 0, 120}]] === {True}");
    audit_chk(id, "BernoulliB[12]", "-691/2730");
    audit_chk(id, "BernoulliB[1]", "-1/2");
    audit_chk(id, "BernoulliB[4999]", "0");
}

/* Builtin-audit regressions: special functions (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* QPochhammer#6 (low, performance): a large positive n expanded an unbounded
 * symbolic product (10^7 factors: minutes and gigabytes), unlike negative n.
 * Both signs now share one expansion bound; past it an exact/symbolic call
 * stays in closed form and an inexact one goes through (a;q)_inf/(a q^n;q)_inf. */
static void test_qpochhammer_large_n(void) {
    audit_chk("QPochhammer#6", "QPochhammer[a, q, 10^7]", "QPochhammer[a, q, 10000000]");
    audit_chk("QPochhammer#6", "QPochhammer[a, q, -10^7]", "QPochhammer[a, q, -10000000]");
    /* (1/2; 1/2)_n converges to (1/2; 1/2)_inf. */
    audit_true("QPochhammer#6",
               "Abs[QPochhammer[0.5, 0.5, 10^7] - QPochhammer[0.5, 0.5]] < 10^-14");
    audit_true("QPochhammer#6",
               "Abs[QPochhammer[0.5, 0.25, 200000] - QPochhammer[0.5, 0.25, 200]] < 10^-14");
    audit_chk("QPochhammer#6", "Expand[QPochhammer[a, q, 2]]", "1 - a - a q + a^2 q");
    audit_true("QPochhammer#6", "QPochhammer[2, 3, 3] == (1 - 2) (1 - 6) (1 - 18)");
}

/* LegendreQ#1 (high): type 3 with an integer degree used the type-1 (Ferrers)
 * log, Log[1 - x], so it was complex for real x > 1.  The type-3 function is
 * real there (DLMF 14.3.7): Q_0 = ArcCoth[x], Q_1 = x ArcCoth[x] - 1,
 * Q_2 = P_2(x) ArcCoth[x] - 3x/2.  Type 1 keeps its cut value. */
static void test_legendreq_type3(void) {
    const char* id = "LegendreQ#1";
    audit_true(id, "LegendreQ[0, 0, 3, x] === Log[1 + x]/2 - Log[-1 + x]/2");
    audit_true(id, "Head[LegendreQ[0, 0, 3, 2.]] === Real && "
                   "Abs[LegendreQ[0, 0, 3, 2.] - N[ArcCoth[2], 20]] < 1.*^-15");
    audit_true(id, "Abs[LegendreQ[1, 0, 3, 2.] - N[2 ArcCoth[2] - 1, 20]] < 1.*^-15");
    audit_true(id, "Abs[LegendreQ[2, 0, 3, 3.] - N[(3 3^2 - 1)/2 ArcCoth[3] - 9/2, 20]] < 1.*^-15");
    audit_true(id, "Head[LegendreQ[3, 0, 3, 2.]] === Real");
    /* m = 1: (x^2 - 1)^(1/2) d/dx Q_2(x), and continuity in the degree */
    audit_true(id, "Abs[LegendreQ[2, 1, 3, 2.] - N[Sqrt[3] D[(3 x^2 - 1)/2 ArcCoth[x] - 3 x/2, x] /. x -> 2, 20]] < 1.*^-14");
    audit_true(id, "Abs[LegendreQ[2, 1, 3, 2.] - LegendreQ[2.0000001, 1, 3, 2.]] < 1.*^-6");
    /* sanity: type 1 (and type 2 on the cut) unchanged */
    audit_true(id, "Abs[LegendreQ[0, 2.] - (N[ArcCoth[2], 20] - N[Pi/2, 20] I)] < 1.*^-15");
    audit_true(id, "LegendreQ[0, x] === -Log[1 - x]/2 + Log[1 + x]/2");
    audit_true(id, "Abs[LegendreQ[2, 1, 2, 0.5] - LegendreQ[2, 1, 0.5]] < 1.*^-15");
}

/* LegendreQ#2/#3/#4/#5: numeric LegendreQ at a large degree.  The closed form
 * P_n L + v_n evaluated at an inexact x cancelled ~n digits (machine 10^112 at
 * n = 1025; a 30-digit tag on -7.5*10^26 at n = 500), the Frobenius series
 * did the same for a non-integer degree, and n > 2000 stayed unevaluated.
 * arb's ball now answers (with enough first-rung bits).  Independent checks:
 * the exact closed form at x = 3/10 under N (precision-tracked), the bound
 * |Q_nu(cos t)| < Sqrt[Pi/(2 nu sin t)]-ish, and the three-term recurrence
 * (nu+1) Q_{nu+1} = (2nu+1) x Q_nu - nu Q_{nu-1} at a non-integer degree. */
static void test_legendreq_large_degree(void) {
    audit_true("LegendreQ#2", "Abs[LegendreQ[1025, 0.3]/N[LegendreQ[1025, 3/10], 20] - 1] < 1.*^-12");
    audit_true("LegendreQ#2", "Abs[LegendreQ[2000, 0.3]/N[LegendreQ[2000, 3/10], 20] - 1] < 1.*^-12");
    audit_true("LegendreQ#2", "Abs[LegendreQ[2000, 0.3]] < 0.03");
    audit_true("LegendreQ#2", "With[{xs = {0.3, -0.6, 0.95}}, "
                              "(LegendreQ[1500, #] & /@ xs) === LegendreQ[1500, xs]]");
    /* sanity: the machine kernel's range and its value at n = 1000 */
    audit_true("LegendreQ#2", "Abs[LegendreQ[1000, 0.3]/N[LegendreQ[1000, 3/10], 20] - 1] < 1.*^-12");
    audit_true("LegendreQ#3", "Abs[LegendreQ[500, 0.3`30] - N[LegendreQ[500, 3/10], 40]] < 10^-28");
    audit_true("LegendreQ#3", "Abs[LegendreQ[200, 0.3`30] - N[LegendreQ[200, 3/10], 40]] < 10^-28");
    audit_true("LegendreQ#3", "Precision[LegendreQ[200, 0.3`30]] == 30.");
    audit_true("LegendreQ#3", "Abs[LegendreQ[5, 2.`25] - N[LegendreQ[5, 2], 25]] < 10^-20");
    audit_true("LegendreQ#4", "Abs[LegendreQ[2001.5, 0.3]] < 0.03 && Abs[LegendreQ[2000.0000001, 0.3]] < 0.03");
    audit_true("LegendreQ#4", "With[{v = 4001/2, x = 3/10}, Abs[(v + 1) LegendreQ[N[v + 1, 30], N[x, 30]] - "
                              "(2 v + 1) x LegendreQ[N[v, 30], N[x, 30]] + v LegendreQ[N[v - 1, 30], N[x, 30]]] < 10^-20]");
    audit_true("LegendreQ#4", "Abs[LegendreQ[2001.5, 0.3] - LegendreQ[2001.5`30, 0.3`30]] < 1.*^-12");
    audit_true("LegendreQ#5", "NumberQ[LegendreQ[2001, 0.3]] && "
                              "Abs[LegendreQ[2001, 0.3] - LegendreQ[2001, 0.3`30]] < 1.*^-13");
}

/* Hypergeometric0F1#1 (high): 0F1(; b; -10^10) came back as 10^86774: the
 * machine kernel declined and the MPFR series, promoted to 53 bits, kept
 * 290000 bits of cancellation as its answer.  arb's 0F1 now answers, and the
 * MPFR series declines (or re-sums) a loss beyond its guard.  Independent
 * check: 0F1(; b; -x) = Gamma[b] x^((1-b)/2) BesselJ[b-1, 2 Sqrt[x]]. */
static void test_hypergeometric0f1_large_negative(void) {
    const char* id = "Hypergeometric0F1#1";
    audit_true(id, "Abs[Hypergeometric0F1[2, -1.*^10]/N[BesselJ[1, 2*10^5]/10^5, 20] - 1] < 1.*^-9");
    audit_true(id, "Abs[Hypergeometric0F1[1, -1.*^10]/N[BesselJ[0, 2*10^5], 20] - 1] < 1.*^-9");
    audit_true(id, "Abs[Hypergeometric0F1[2, -1.*^8]/N[BesselJ[1, 2*10^4]/10^4, 20] - 1] < 1.*^-10");
    audit_true(id, "Abs[N[Hypergeometric0F1[2, -10^10], 20] - N[BesselJ[1, 2*10^5]/10^5, 25]] < 10^-27");
    /* sanity: moderate and large positive z unchanged */
    audit_true(id, "Abs[Hypergeometric0F1[2, -1.*^6]/N[BesselJ[1, 2000]/1000, 20] - 1] < 1.*^-11");
    audit_true(id, "Abs[Hypergeometric0F1[2.5, 3.] - N[Hypergeometric0F1[5/2, 3], 20]] < 1.*^-14");
    audit_true(id, "Abs[Hypergeometric0F1[2, 1.*^7]/N[BesselI[1, 2 Sqrt[10^7]]/Sqrt[10^7], 20] - 1] < 1.*^-9");
    audit_true(id, "Abs[N[Hypergeometric0F1[1/3, -50], 30] - N[Gamma[1/3] 50^(1/3) BesselJ[-2/3, 2 Sqrt[50]], 35]] < 10^-27");
    audit_true(id, "Abs[N[Hypergeometric1F1[1, 2, -40], 30] - N[(1 - Exp[-40])/40, 35]] < 10^-28");
}

/* HypergeometricU#1/#2: z = 0.  The terminating polynomial built its constant
 * term as (-z)^0, i.e. 0^0 at z = 0; and at an inexact zero the series gave a
 * finite -0.006 at the pole of U(1.5, 7, z).  DLMF 13.2.16-22: U(-m, b, 0) =
 * (-1)^m (b)_m; U(a, b, 0) = Gamma(1-b)/Gamma(a-b+1) for Re b < 1, unbounded
 * for Re b >= 1.  Independent check of the finite value: U at z = 10^-45. */
static void test_hypergeometricu_zero(void) {
    audit_chk("HypergeometricU#1", "HypergeometricU[-1, 2, 0]", "-2");
    audit_true("HypergeometricU#1", "HypergeometricU[-2, b, 0] === b + b^2");
    audit_true("HypergeometricU#1", "HypergeometricU[-1, 2, 0.] === -2.");
    audit_true("HypergeometricU#1", "HypergeometricU[-2, b, z] === b + b^2 - 2 z - 2 b z + z^2");
    audit_chk("HypergeometricU#2", "HypergeometricU[1.5, 7, 0.]", "ComplexInfinity");
    audit_chk("HypergeometricU#2", "HypergeometricU[3/2, 7, 0]", "ComplexInfinity");
    audit_chk("HypergeometricU#2", "HypergeometricU[1.5, 7, 0.`20]", "ComplexInfinity");
    audit_chk("HypergeometricU#2", "HypergeometricU[1, 1, 0]", "ComplexInfinity");
    audit_true("HypergeometricU#2", "HypergeometricU[1/2, 1/3, 0] === Gamma[2/3]/Gamma[7/6]");
    audit_true("HypergeometricU#2", "Abs[N[HypergeometricU[1/2, 1/3, 10^-45], 30] - "
                                    "N[Gamma[2/3]/Gamma[7/6], 30]] < 10^-25");
    audit_true("HypergeometricU#2", "Abs[HypergeometricU[0.5, 1/3, 0.] - N[Gamma[2/3]/Gamma[7/6]]] < 1.*^-15");
    /* sanity: a nonzero z unchanged */
    audit_true("HypergeometricU#2", "Abs[HypergeometricU[1.5, 7, 2.] - N[HypergeometricU[3/2, 7, 2], 20]] < 1.*^-13");
}

/* WhittakerW#1/#2: inexact z = 0.  Only an exact 0 took the special value;
 * 0. reached z^(m+1/2) U(.., 0.) and gave 0.0 (or, for m < 0, leaked
 * Power::infy).  W_{k,m}(z) ~ Gamma(2|m|)/Gamma(1/2+|m|-k) z^(1/2-|m|) (DLMF
 * 13.14.18): ComplexInfinity for |Re m| > 1/2, 1/Gamma(1-k) at m = 1/2, 0
 * below.  Independent check: W at z = 10^-12 against the limit. */
static void test_whittakerw_zero(void) {
    const char* id = "WhittakerW#1";
    audit_chk(id, "WhittakerW[2, 3, 0.]", "ComplexInfinity");
    audit_chk(id, "WhittakerW[1., 1., 0.]", "ComplexInfinity");
    audit_chk(id, "WhittakerW[2, 3, 0.`20]", "ComplexInfinity");
    audit_true(id, "Abs[WhittakerW[1/3, 0.5, 0.] - N[1/Gamma[2/3]]] < 1.*^-15");
    audit_true(id, "Abs[WhittakerW[1/3, 0.5, 1.*^-12] - WhittakerW[1/3, 0.5, 0.]] < 1.*^-9");
    audit_true(id, "WhittakerW[1/3, 1/4, 0.] == 0");
    audit_chk(id, "WhittakerW[2, 3, 0]", "ComplexInfinity");
    audit_chk("WhittakerW#2", "Check[WhittakerW[2, -3, 0.], \"fired\"]", "ComplexInfinity");
    /* sanity: a nonzero z unchanged, and W even in m */
    audit_true(id, "Abs[WhittakerW[2, 3, 4.5] - N[WhittakerW[2, 3, 9/2], 20]] < 1.*^-13");
    audit_true(id, "WhittakerW[2, -3, 4.5] == WhittakerW[2, 3, 4.5]");
}

/* WeberE#1 (high): WeberE[n, 0] for odd n built the constant term of the
 * integer-order polynomial as z^0, so z = 0 gave 0^0.  E_nu(0) =
 * (1 - Cos[nu Pi])/(nu Pi), checked against the machine value. */
static void test_webere_zero(void) {
    const char* id = "WeberE#1";
    audit_chk(id, "WeberE[1, 0]", "2/Pi");
    audit_chk(id, "WeberE[3, 0]", "2/(3 Pi)");
    audit_chk(id, "WeberE[-1, 0]", "-2/Pi");
    audit_chk(id, "WeberE[2, 0]", "0");
    audit_true(id, "Abs[N[WeberE[5, 0]] - WeberE[5, 0.]] < 1.*^-15");
    /* sanity: the symbolic polynomial at nonzero z unchanged */
    audit_true(id, "WeberE[3, x] === 2/(3 Pi) + 2 x^2/(15 Pi) - StruveH[3, x]");
}

/* MathieuC#1-#3, MathieuS#1-#3: Mathieu functions at a characteristic value
 * for larger q.  Hill's determinant (cancelling ~(q^2/16)^k/(k!)^4) no longer
 * snapped nu to the integer, so the DLMF 1/Sqrt[2] was lost (or the plan
 * declined); and the sign Re(c_0) > 0 flips where c_0 passes through zero.
 * Now a Sturm count against the characteristic-value matrices identifies
 * a_r / b_r, and the sign follows DLMF 28.2(vi): ce_{2j}(Pi/2)(-1)^j > 0,
 * ce_{2j+1}'(Pi/2)(-1)^(j+1) > 0, se_{2j+1}(Pi/2)(-1)^j > 0,
 * se_{2j+2}'(Pi/2)(-1)^(j+1) > 0.  Independent checks: the normalization
 * integral (= Pi) by NIntegrate, the DLMF sign conditions, continuity in q,
 * and machine == 30-digit. */
static void test_mathieu_characteristic(void) {
    audit_run("aA[n_, q_] := MathieuCharacteristicA[n, q]; aB[n_, q_] := MathieuCharacteristicB[n, q]");
    audit_true("MathieuC#1", "Abs[NIntegrate[MathieuC[aA[2, 35.], 35., z]^2, {z, 0, 2 Pi}] - Pi] < 1.*^-8");
    audit_true("MathieuC#1", "Abs[NIntegrate[MathieuC[aA[0, 40.], 40., z]^2, {z, 0, 2 Pi}] - Pi] < 1.*^-8");
    audit_true("MathieuC#1", "Abs[MathieuC[aA[2, N[35, 30]], N[35, 30], N[3/10, 30]] - "
                             "MathieuC[aA[2, 35.], 35., 0.3]] < 1.*^-13");
    audit_true("MathieuC#2", "And @@ Table[MathieuC[aA[4, q], q, N[Pi/2]] > 0, {q, {10., 14., 14.5, 20., 25., 40.}}]");
    audit_true("MathieuC#2", "Abs[MathieuC[aA[4, 14.5], 14.5, N[Pi/2]] - MathieuC[aA[4, 14.4], 14.4, N[Pi/2]]] < 0.01");
    audit_true("MathieuC#3", "MathieuC[aA[3, 20.], 20., 0.1] > 0");
    audit_true("MathieuC#3", "And @@ Table[MathieuCPrime[aA[3, q], q, N[Pi/2]] > 0, {q, {5., 15., 20., 30., 35., 40.}}]");
    audit_true("MathieuC#3", "Abs[NIntegrate[MathieuC[aA[3, 20.], 20., z]^2, {z, 0, 2 Pi}] - Pi] < 1.*^-8");
    audit_true("MathieuS#1", "NumberQ[MathieuS[aB[2, 40.], 40., 0.4]] && MathieuS[aB[2, 40.], 40., 0.4] > 0");
    audit_true("MathieuS#1", "MathieuSPrime[aB[2, 40.], 40., N[Pi/2]] < 0");
    audit_true("MathieuS#1", "Abs[NIntegrate[MathieuS[aB[2, 40.], 40., z]^2, {z, 0, 2 Pi}] - Pi] < 1.*^-8");
    audit_true("MathieuS#1", "Abs[MathieuS[aB[2, N[40, 30]], N[40, 30], N[2/5, 30]] - "
                             "MathieuS[aB[2, 40.], 40., 0.4]] < 1.*^-13");
    audit_true("MathieuS#2", "Abs[NIntegrate[MathieuS[aB[1, 15.], 15., z]^2, {z, 0, 2 Pi}] - Pi] < 1.*^-8");
    audit_true("MathieuS#3", "And @@ Table[MathieuS[aB[5, q], q, N[Pi/2]] > 0, {q, {10., 20., 30., 35., 40.}}]");
    audit_true("MathieuS#3", "And @@ Table[MathieuS[aB[1, q], q, N[Pi/2]] > 0, {q, {5., 15., 20., 25., 30.}}]");
    audit_true("MathieuS#3", "MathieuS[aB[5, 35.], 35., 0.4] > 0 && MathieuS[aB[2, 30.], 30., 0.4] > 0");
    /* the wrong parity is still identically 0 */
    audit_true("MathieuS#1", "MathieuS[aA[2, 35.], 35., 0.4] == 0 && MathieuC[aB[2, 40.], 40., 0.4] == 0");
    /* sanity: small q, negative q, generic and complex a unchanged */
    audit_true("MathieuC#1", "MathieuC[aA[2, 1.], 1., 0.3] == 0.9693734655943552");
    audit_true("MathieuC#1", "MathieuC[aA[2, -5.], -5., 0.3] == 0.232723423823903");
    audit_true("MathieuC#1", "MathieuC[1.5, 2., 0.7] == 1.6751744989286563");
    audit_true("MathieuC#1", "MathieuC[2 + I, 1. + 0.5 I, 0.3] == 1.1561237520754863 - 0.1593777733640514 I");
}

/* ================================================================== */
/* Group A2: Gamma, LogGamma, Pochhammer, HarmonicNumber, LerchPhi,    */
/* InverseErf(c), AiryAi, DiracDelta, HaarWavelet, JacobiP.            */
/* ================================================================== */

/* Gamma#1 (high): Gamma[a, 0] = Gamma[a] was applied for Re a <= 0, where
 * Int_0^Inf t^(a-1) e^-t dt diverges at 0 (Gamma(a, z) ~ -z^a/a). */
static void test_a2_gamma_incomplete_at_zero(void) {
    const char* id = "Gamma#1";
    audit_chk(id, "Gamma[-1/2, 0]", "ComplexInfinity");
    audit_chk(id, "Gamma[-1.5, 0]", "ComplexInfinity");
    audit_chk(id, "Gamma[0, 0]", "ComplexInfinity");
    audit_chk(id, "Gamma[-1/2 + I, 0]", "ComplexInfinity");
    audit_chk(id, "Gamma[-3/2, 0, 1]", "ComplexInfinity");
    /* Re a = 0, a != 0: z^a oscillates boundedly, no limit */
    audit_chk(id, "Gamma[I, 0]", "Gamma[I, 0]");
    /* the blow-up the limit follows: ~ 2 z^(-1/2) */
    audit_true(id, "Gamma[-1/2, 1.*^-12] > 10^6");
    /* sanity: Re a > 0 and symbolic a keep Gamma[a] */
    audit_chk(id, "Gamma[1/2, 0]", "Sqrt[Pi]");
    audit_chk(id, "Gamma[3, 0]", "2");
    audit_chk(id, "Gamma[a, 0]", "Gamma[a]");
    audit_chk(id, "Gamma[2, 3]", "4/E^3");
}

/* LogGamma#1 (high): machine complex LogGamma lost up to 8 digits for small |z|
 * with Re z > 0: Log z went through log1p(z - 1), whose argument |z|^2 - 1
 * cancels to nothing as z -> 0. Reference values: N of the exact argument. */
static void test_a2_loggamma_small_complex(void) {
    const char* id = "LogGamma#1";
    audit_true(id, "Abs[LogGamma[0.00001 + 0.00001 I] - "
                   "N[LogGamma[1/100000 + I/100000], 30]] < 1.*^-14");
    audit_true(id, "Abs[LogGamma[0.001 + 0.001 I] - "
                   "N[LogGamma[1/1000 + I/1000], 30]] < 1.*^-14");
    audit_true(id, "Abs[LogGamma[0.2 - 0.3 I] - N[LogGamma[2/10 - 3/10 I], 30]] < 1.*^-14");
    /* scalar == packed == NDArray == Compile */
    audit_true(id, "With[{w = {0.00001 + 0.00001 I, 0.001 + 0.001 I, 0.3 + 2 I}}, "
                   "(LogGamma /@ w) === LogGamma[w] === Normal[LogGamma[NDArray[w]]] === "
                   "(Compile[{{z, _Complex}}, LogGamma[z]] /@ w)]");
    /* sanity: near z = 1 (where log1p is kept) */
    audit_true(id, "Abs[LogGamma[0.7 + 0.3 I] - N[LogGamma[7/10 + 3/10 I], 30]] < 1.*^-15");
    audit_true(id, "Abs[LogGamma[1.00001 + 0.00001 I]/"
                   "LogGamma[SetPrecision[1.00001 + 0.00001 I, 40]] - 1] < 1.*^-14");
}

/* Pochhammer#1 (high): an arbitrary-precision a at a pole of the negative-n
 * product multiplied out an MPFR zero and returned Overflow[]. */
static void test_a2_pochhammer_mpfr_pole(void) {
    const char* id = "Pochhammer#1";
    audit_chk(id, "Pochhammer[N[2, 30], -3]", "ComplexInfinity");
    audit_chk(id, "Pochhammer[N[1, 30], -1]", "ComplexInfinity");
    audit_chk(id, "Pochhammer[2., -3]", "ComplexInfinity");
    audit_chk(id, "Pochhammer[2, -3]", "ComplexInfinity");
    /* sanity: off the pole, 1/((a-1)(a-2)(a-3)) */
    audit_true(id, "Abs[Pochhammer[N[5, 30], -3] - 1/24] < 10^-28");
    audit_chk(id, "Pochhammer[10, -3]", "1/504");
    audit_true(id, "Abs[Pochhammer[N[1/2, 30], -3] + 8/15] < 10^-28");
}

/* HarmonicNumber#1 (high): HarmonicNumber[Infinity, r] answered Zeta[r] for
 * r <= 1, where the sum of positive terms diverges. */
static void test_a2_harmonicnumber_infinity(void) {
    const char* id = "HarmonicNumber#1";
    audit_chk(id, "HarmonicNumber[Infinity, -1]", "Infinity");
    audit_chk(id, "HarmonicNumber[Infinity, 0]", "Infinity");
    audit_chk(id, "HarmonicNumber[Infinity, 1/2]", "Infinity");
    audit_chk(id, "HarmonicNumber[Infinity]", "Infinity");
    audit_chk(id, "HarmonicNumber[Infinity, 0.5]", "Infinity");
    audit_chk(id, "HarmonicNumber[Infinity, 1/2 + I]", "ComplexInfinity");
    audit_chk(id, "HarmonicNumber[Infinity, 1 + I]", "HarmonicNumber[Infinity, 1 + I]");
    /* sanity: the convergent orders and symbolic r */
    audit_chk(id, "HarmonicNumber[Infinity, 2]", "1/6 Pi^2");
    audit_chk(id, "HarmonicNumber[Infinity, 3 + I]", "Zeta[3 + I]");
    audit_chk(id, "HarmonicNumber[Infinity, r]", "Zeta[r]");
}

/* LerchPhi#1 (high): exact reductions at Re a < 0 used the principal power
 * (k + a)^-s instead of the symmetric ((k + a)^2)^(-s/2) the numeric path
 * (and the convention) uses. Reference: the series itself, summed exactly or
 * numerically, and the machine path. */
static void test_a2_lerchphi_negative_a(void) {
    const char* id = "LerchPhi#1";
    /* k=0: 1/8, k=1: 1/2, k=2 dropped, k>=3: z^2 Li3(z) */
    audit_true(id, "PossibleZeroQ[LerchPhi[1/2, 3, -2] - (5/8 + PolyLog[3, 1/2]/4)]");
    audit_true(id, "Abs[N[LerchPhi[1/2, 3, -2]] - LerchPhi[0.5, 3., -2.]] < 1.*^-14");
    audit_chk(id, "LerchPhi[1/2, -1, -2]", "3");
    audit_chk(id, "LerchPhi[0, 1/2, -1/2]", "Sqrt[2]");
    audit_true(id, "Abs[N[LerchPhi[1/2, 1/2, -1]] - LerchPhi[0.5, 0.5, -1.]] < 1.*^-14");
    /* non-integer negative a: Sum |k - 1/2| 2^-k = 2;  Sum (-1)^k/|k - 1/2| */
    audit_chk(id, "LerchPhi[1/2, -1, -1/2]", "2");
    audit_true(id, "Abs[N[LerchPhi[-1, 1, -1/2]] - "
                   "(2 + NSum[(-1)^k/(k - 1/2), {k, 1, Infinity}])] < 1.*^-10");
    audit_true(id, "Abs[N[LerchPhi[1/3, 2, -3/2]] - "
                   "N[Sum[(1/3)^k/((k - 3/2)^2), {k, 0, 200}], 20]] < 1.*^-14");
    /* the k = -a term of a z = 0 call is dropped, not 0^0 */
    audit_chk(id, "LerchPhi[0, 2, -1]", "1");
    /* sanity: positive a and irreducible negative a unchanged */
    audit_true(id, "PossibleZeroQ[LerchPhi[1/2, 3, 2] - 4 (PolyLog[3, 1/2] - 1/2)]");
    audit_chk(id, "LerchPhi[z, s, -1/2]", "LerchPhi[z, s, -1/2]");
    audit_chk(id, "LerchPhi[1/2, s, 0]", "PolyLog[s, 1/2]");
}

/* InverseErf#1/#2, InverseErfc#1 (high): the arbitrary-precision Newton was
 * seeded from s rounded to a double, which is 1.0 (or 2.0) next to the
 * endpoints: Indeterminate, or digits wrong far below the tagged precision.
 * Checked by the residual of Erf / Erfc and against N of the exact argument. */
static void test_a2_inverse_erf_tails(void) {
    const char* id = "InverseErf#1";
    audit_true(id, "Abs[InverseErf[N[1 - 10^-20, 40]] - 6.6015806223551425615163916] < 10^-24");
    audit_true(id, "Abs[InverseErf[N[-1 + 10^-20, 40]] + 6.6015806223551425615163916] < 10^-24");
    audit_true(id, "Abs[N[InverseErf[1 - 10^-20], 30] - 6.60158062235514256151639163242] < 10^-27");
    audit_true(id, "Abs[Erf[InverseErf[N[1 - 10^-17, 40]]] - N[1 - 10^-17, 40]] < 10^-38");
    id = "InverseErf#2";
    audit_true(id, "Abs[InverseErf[N[1 - 10^-16, 40]] - "
                   "5.87237009045396314506648428187] < 10^-25");
    audit_true(id, "Abs[InverseErf[N[1 - 10^-15, 40]] - "
                   "5.67584634746764698500497080641] < 10^-25");
    audit_true(id, "Abs[Erf[InverseErf[N[1 - 10^-16, 40]]] - N[1 - 10^-16, 40]] < 10^-38");
    /* sanity: the interior and the machine path */
    audit_true(id, "Abs[Erf[InverseErf[N[1/3, 30]]] - 1/3] < 10^-28");
    audit_true(id, "Abs[Erf[InverseErf[N[7/10, 30]]] - 7/10] < 10^-28");
    audit_true(id, "Abs[InverseErf[0.7] - InverseErf[N[7/10, 30]]] < 1.*^-15");

    id = "InverseErfc#1";
    audit_true(id, "Abs[InverseErfc[N[2 - 10^-30, 40]] + 8.1486162231698646] < 10^-15");
    audit_true(id, "Abs[Erfc[InverseErfc[N[2 - 10^-30, 40]]] - N[2 - 10^-30, 40]] < 10^-38");
    audit_true(id, "Abs[InverseErfc[N[2 - 10^-16, 40]] + 5.8723700904539631450664] < 10^-21");
    /* below the double range (InverseErfc#2 shape, same root cause) */
    audit_true(id, "30 < InverseErfc[N[10^-400, 30]] < 304/10");
    audit_true(id, "Abs[Erfc[InverseErfc[N[10^-400, 30]]]/N[10^-400, 30] - 1] < 10^-27");
    /* sanity */
    audit_true(id, "Abs[InverseErfc[N[10^-30, 40]] - 8.1486162231698646073845666606481] < 10^-30");
    audit_true(id, "Abs[InverseErfc[1.5] + 0.4769362762044699] < 1.*^-15");
}

/* AiryAi#1 (high): the large-negative-x machine asymptotics took sin/cos of
 * zeta = zh + zl to first order in zl, and zh + zl itself carries the phase
 * only to ~2^-104 relative -- AiryAi[-1.*^15] was 69. The arbitrary-precision
 * path on the same (exact) double is the reference. */
static void test_a2_airyai_large_negative(void) {
    const char* id = "AiryAi#1";
    audit_true(id, "Abs[AiryAi[-1.*^12] - AiryAi[-10^12`30]] < 1.*^-18");
    audit_true(id, "Abs[AiryAi[-1.*^15] - AiryAi[-10^15`30]] < 1.*^-19");
    audit_true(id, "Abs[AiryAi[-1.*^20] - AiryAi[-10^20`30]] < 1.*^-20");
    audit_true(id, "Abs[AiryAi[-1.*^8]/AiryAi[-10^8`30] - 1] < 1.*^-14");
    audit_true(id, "Abs[AiryAi[-1.*^10]/AiryAi[-10^10`30] - 1] < 1.*^-13");
    audit_true(id, "Abs[AiryBi[-1.*^8]/AiryBi[-10^8`30] - 1] < 1.*^-14");
    /* scalar == packed == NDArray == Compile, across the decline threshold */
    audit_true(id, "With[{v = {-1.*^12, -1.*^8, -5., -1.*^15}}, "
                   "(AiryAi /@ v) === AiryAi[v] === Normal[AiryAi[NDArray[v]]] === "
                   "(Compile[{{x, _Real}}, AiryAi[x]] /@ v)]");
    /* sanity: moderate arguments */
    audit_true(id, "Abs[AiryAi[-50.] - AiryAi[-50`30]] < 1.*^-16");
    audit_true(id, "Abs[AiryAi[-5.] - AiryAi[-5`30]] < 1.*^-16");
}

/* DiracDelta#1/#2 (high): definite integrals with a nonlinear delta argument,
 * a product of deltas or a multivariate delta declined the sifting and a later
 * method answered 0. Reference: the composition rule Sum_r h(r)/|g'(r)| and
 * DiracDelta[x, y] = DiracDelta[x] DiracDelta[y]. */
static void test_a2_diracdelta_integrals(void) {
    const char* id = "DiracDelta#1";
    audit_chk(id, "Integrate[DiracDelta[x^2 - 1], {x, -2, 2}]", "1");
    audit_chk(id, "Integrate[DiracDelta[x^2 - 4], {x, -3, 3}]", "1/2");
    audit_chk(id, "Integrate[DiracDelta[Sin[x]], {x, -1, 1}]", "1");
    audit_chk(id, "Integrate[x^2 DiracDelta[x^2 - 1], {x, 0, 2}]", "1/2");
    audit_chk(id, "Integrate[DiracDelta[x^2 + 1], {x, -1, 1}]", "0");
    /* a multiple root: DiracDelta[x^2] is undefined, stays unevaluated */
    audit_chk(id, "Integrate[DiracDelta[x^2], {x, -1, 1}]",
                  "Integrate[DiracDelta[x^2], {x, -1, 1}]");
    audit_chk(id, "Integrate[DiracDelta[x^3], {x, -1, 1}]",
                  "Integrate[DiracDelta[x^3], {x, -1, 1}]");
    id = "DiracDelta#2";
    audit_chk(id, "Integrate[DiracDelta[x] DiracDelta[y], {x, -1, 1}]", "DiracDelta[y]");
    audit_chk(id, "Integrate[DiracDelta[x] DiracDelta[y], {x, -1, 1}, {y, -1, 1}]", "1");
    audit_chk(id, "Integrate[DiracDelta[x, y], {x, -1, 1}, {y, -1, 1}]", "1");
    audit_chk(id, "Integrate[DiracDelta[x - 1/2, y] Exp[x], {x, 0, 1}, {y, -1, 1}]", "Sqrt[E]");
    /* sanity: the linear sifting */
    audit_chk(id, "Integrate[DiracDelta[x - 1] Exp[x], {x, 0, 2}]", "E");
    audit_chk(id, "Integrate[1 + DiracDelta[t - 2], {t, 0, 3}]", "4");
    audit_chk(id, "Integrate[DiracDelta[x] Cos[x], {x, -Infinity, Infinity}]", "1");
}

/* HaarWavelet#1 (high): Integrate split a Piecewise at no breakpoint when its
 * condition was an Inequality (0 <= x < 1) or a chained Less[1/2, x, 1], and
 * integrated the whole interval on one branch. */
static void test_a2_haar_norms(void) {
    const char* id = "HaarWavelet#1";
    audit_chk(id, "Integrate[WaveletPhi[HaarWavelet[], x], {x, -1, 2}]", "1");
    audit_chk(id, "Integrate[WaveletPsi[HaarWavelet[], x]^2, {x, 0, 1}]", "1");
    audit_chk(id, "Integrate[WaveletPsi[HaarWavelet[], x], {x, 0, 3/4}]", "1/4");
    audit_chk(id, "Integrate[Piecewise[{{1, 0 <= x < 1}}, 0], {x, -1, 2}]", "1");
    audit_chk(id, "Integrate[Piecewise[{{1, Less[0, x, 1]}}, 0], {x, -1, 2}]", "1");
    /* sanity: two-argument relations */
    audit_chk(id, "Integrate[Piecewise[{{x, x < 1}}, 2], {x, 0, 3}]", "9/2");
}

/* JacobiP#1 (high): a degenerate three-term recurrence (a + b <= -2) fell back
 * to Binomial[n+a, n] 2F1(..., a+1, ...), a 0 * Infinity read as an exact 0.
 * Reference: the exact polynomial JacobiP[3, -1, -1, x] = (x^3 - x)/2. */
static void test_a2_jacobip_degenerate(void) {
    const char* id = "JacobiP#1";
    audit_chk(id, "JacobiP[3, -1, -1, x]", "-1/2 x + 1/2 x^3");
    audit_chk(id, "JacobiP[3, -1, -1, 0.5]", "-0.1875");
    audit_chk(id, "JacobiP[3, -2, -1, 0.5]", "0.046875");
    audit_true(id, "Abs[JacobiP[3, -1, -1, 0.5`30] + 3/16] < 10^-28");
    audit_true(id, "Abs[JacobiP[3, -1, -1, 0.5 + 0.1 I] - (x^3 - x)/2 /. x -> 0.5 + 0.1 I] < 1.*^-15");
    /* JacobiP#2's shape (real a + b = -2): no messages, the same value */
    audit_chk(id, "Check[JacobiP[3, -1., -1., 0.5], \"msg\"]", "-0.1875");
    /* sanity: a regular recurrence */
    audit_true(id, "Abs[JacobiP[3, 1/2, 1/3, 0.5] - N[JacobiP[3, 1/2, 1/3, 1/2]]] < 1.*^-15");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("BesselJ#4", "Quiet[Check[BesselJ[x], fired, BesselJ::argr]]", "fired");
    audit_chk("BesselJ#4", "Quiet[BesselJ[x]]", "BesselJ[x]");
    audit_chk("Gamma#3", "Quiet[Check[Gamma[], fired, Gamma::argb]]", "fired");
    audit_chk("Gamma#3", "Quiet[Gamma[]]", "Gamma[]");
    audit_chk("AiryBiPrimeZero#3", "Quiet[Check[AiryBiPrimeZero[1, 2, 3], fired, AiryBiPrimeZero::argt]]", "fired");
    audit_chk("AiryBiPrimeZero#3", "Quiet[AiryBiPrimeZero[1, 2, 3]]", "AiryBiPrimeZero[1, 2, 3]");
    audit_chk("AngerJ#1", "Quiet[Check[AngerJ[1, 2, 3, 4], fired, AngerJ::argt]]", "fired");
    audit_chk("AngerJ#1", "Quiet[AngerJ[1, 2, 3, 4]]", "AngerJ[1, 2, 3, 4]");
    audit_chk("BesselJZero#5", "Quiet[Check[BesselJZero[1], fired, BesselJZero::argtu]]", "fired");
    audit_chk("BesselJZero#5", "Quiet[BesselJZero[1]]", "BesselJZero[1]");
    audit_chk("MathieuCharacteristicA#1", "Quiet[Check[MathieuCharacteristicA[1, 2, 3], fired, MathieuCharacteristicA::argrx]]", "fired");
    audit_chk("MathieuCharacteristicA#1", "Quiet[MathieuCharacteristicA[1, 2, 3]]", "MathieuCharacteristicA[1, 2, 3]");
    audit_chk("QPochhammer#4", "Quiet[Check[QPochhammer[a, b, c, d], fired, QPochhammer::argb]]", "fired");
    audit_chk("QPochhammer#4", "Quiet[QPochhammer[a, b, c, d]]", "QPochhammer[a, b, c, d]");
    audit_chk("WaveletPhi#2", "Quiet[Check[WaveletPhi[DaubechiesWavelet[2], 1., Foo -> 2], fired, WaveletPhi::optx]]", "fired");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* BesselJ#1: BesselJ at the real infinities */
    audit_true("BesselJ#1", "BesselJ[0, Infinity] === 0");
    audit_true("BesselJ#1", "BesselJ[2, -Infinity] === 0");
    audit_true("BesselJ#1", "BesselJ[1/2, Infinity] === 0");
    /* BesselK#1: BesselK at Infinity */
    audit_true("BesselK#1", "BesselK[0, Infinity] === 0");
    audit_true("BesselK#1", "BesselK[2, Infinity] === 0");
    /* ArcSin#1: ArcSin at +-Infinity (branch of ArcSin[10.^8] = 1.57 - 19.1 I) */
    audit_true("ArcSin#1", "ArcSin[Infinity] === DirectedInfinity[-I]");
    audit_true("ArcSin#1", "ArcSin[-Infinity] === DirectedInfinity[I]");
    audit_true("ArcSin#1", "ArcSin[I Infinity] === DirectedInfinity[I]");
    audit_true("ArcSin#1", "Im[ArcSin[10.^8]] < -19");
    /* ChebyshevT#1: ChebyshevT at the infinities is its leading term */
    audit_true("ChebyshevT#1", "Check[ChebyshevT[3, Infinity], \"m\"] === Infinity");
    audit_true("ChebyshevT#1", "Check[ChebyshevT[4, -Infinity], \"m\"] === Infinity");
    audit_true("ChebyshevT#1", "ChebyshevT[3, -Infinity] === -Infinity");
    audit_true("ChebyshevT#1", "ChebyshevT[0, Infinity] === 1");
    /* ChebyshevU#1: ChebyshevU at the infinities */
    audit_true("ChebyshevU#1", "Check[ChebyshevU[3, Infinity], \"m\"] === Infinity");
    audit_true("ChebyshevU#1", "Check[ChebyshevU[3, -Infinity], \"m\"] === -Infinity");
    audit_true("ChebyshevU#1", "ChebyshevU[2, I Infinity] === -Infinity");
    /* EllipticPi#1: EllipticPi[n, 1] is a pole like EllipticK[1] */
    audit_true("EllipticPi#1", "EllipticPi[0.5, 1.] === ComplexInfinity");
    audit_true("EllipticPi#1", "EllipticPi[-1., 1.] === ComplexInfinity");
    audit_true("EllipticPi#1", "EllipticK[1.] === ComplexInfinity");
    audit_true("EllipticPi#1", "Abs[EllipticPi[0.5, 0.5] - 2.7012877620953506] < 10.^-12");
    /* ExpIntegralE#2: ExpIntegralE[n, 0] for n <= 1 */
    audit_true("ExpIntegralE#2", "ExpIntegralE[1, 0] === ComplexInfinity");
    audit_true("ExpIntegralE#2", "ExpIntegralE[-2, 0] === ComplexInfinity");
    audit_true("ExpIntegralE#2", "ExpIntegralE[3, 0] === 1/2");
    /* HermiteH#1: HermiteH at the infinities */
    audit_true("HermiteH#1", "Check[HermiteH[3, Infinity], \"m\"] === Infinity");
    audit_true("HermiteH#1", "Check[HermiteH[3, -Infinity], \"m\"] === -Infinity");
    /* Hypergeometric1F1#2: a lower parameter at a non-positive integer is a pole */
    audit_true("Hypergeometric1F1#2", "Hypergeometric1F1[1, 0, z] === ComplexInfinity");
    audit_true("Hypergeometric1F1#2", "Hypergeometric1F1[1, -2, z] === ComplexInfinity");
    audit_true("Hypergeometric1F1#2", "Hypergeometric1F1[1, 0, 1.] === ComplexInfinity");
    audit_true("Hypergeometric1F1#2", "Hypergeometric1F1[1, 0, 0] === 1");
    audit_true("Hypergeometric1F1#2", "Expand[Hypergeometric1F1[-1, -2, z]] === 1 + z/2");
    /* InverseJacobiSN#1: InverseJacobiSN[1.0, 1] is the logarithmic singularity, not Overflow[] */
    audit_true("InverseJacobiSN#1", "MatchQ[Check[InverseJacobiSN[1.0, 1], \"m\"], Infinity | ComplexInfinity]");
    audit_true("InverseJacobiSN#1", "InverseJacobiSN[1, 1] === ComplexInfinity");
    /* Pochhammer#2: Pochhammer at a Gamma pole */
    audit_true("Pochhammer#2", "Pochhammer[1/2, -3/2] === ComplexInfinity");
    audit_true("Pochhammer#2", "Pochhammer[-5/2, 1/2] === ComplexInfinity");
    /* PolyGamma#1: PolyGamma at an inexact non-positive integer */
    audit_true("PolyGamma#1", "Check[PolyGamma[0, -1.], \"m\"] === ComplexInfinity");
    audit_true("PolyGamma#1", "Check[PolyGamma[1, 0.], \"m\"] === ComplexInfinity");
    audit_true("PolyGamma#1", "Check[PolyGamma[3, -2.], \"m\"] === ComplexInfinity");
    audit_true("PolyGamma#1", "Check[PolyGamma[0, -1.`30], \"m\"] === ComplexInfinity");
    /* QPochhammer#2: QPochhammer[a, q, Infinity] is QPochhammer[a, q] */
    audit_true("QPochhammer#2", "QPochhammer[a, q, Infinity] === QPochhammer[a, q]");
    audit_true("QPochhammer#2", "QPochhammer[0.5, 0.3, Infinity] === QPochhammer[0.5, 0.3]");
    audit_true("QPochhammer#2", "Abs[QPochhammer[0.5, 0.3, Infinity] - 0.3980822043018776] < 10.^-14");
    /* Surd#1: Surd of an infinity with a negative index */
    audit_true("Surd#1", "Surd[Infinity, -3] === 0");
    audit_true("Surd#1", "Surd[-Infinity, -3] === 0");
    audit_true("Surd#1", "Surd[Infinity, -2] === 0");
    audit_true("Surd#1", "Surd[Infinity, 3] === Infinity");
    /* WaveletPhi#1: WaveletPhi at Indeterminate */
    audit_true("WaveletPhi#1", "WaveletPhi[DaubechiesWavelet[2], Indeterminate] === Indeterminate");
    audit_true("WaveletPhi#1", "WaveletPhi[HaarWavelet[], Indeterminate] === Indeterminate");
    /* BesselK#2: BesselK at an inexact zero */
    audit_true("BesselK#2", "BesselK[0, 0.] === Infinity");
    audit_true("BesselK#2", "BesselK[1, 0.] === ComplexInfinity");
    audit_true("BesselK#2", "BesselK[0, 0] === Infinity");
    /* EllipticE#2: EllipticE[-Infinity] */
    audit_true("EllipticE#2", "EllipticE[-Infinity] === Infinity");
    audit_true("EllipticE#2", "EllipticE[-10.^10] > 10^4");
    /* EllipticK#1: EllipticK[Indeterminate] */
    audit_true("EllipticK#1", "EllipticK[Indeterminate] === Indeterminate");
    /* EllipticK#2: EllipticK decays in every direction */
    audit_true("EllipticK#2", "EllipticK[-Infinity] === 0");
    audit_true("EllipticK#2", "EllipticK[I Infinity] === 0");
    audit_true("EllipticK#2", "EllipticK[Infinity] === 0");
    /* Erf#2: Erf[ComplexInfinity] has no limit */
    audit_true("Erf#2", "Erf[ComplexInfinity] === Indeterminate");
    audit_true("Erf#2", "Erf[Infinity] === 1");
    /* Zeta#2: Zeta with an Indeterminate argument */
    audit_true("Zeta#2", "Zeta[Indeterminate] === Indeterminate");
    audit_true("Zeta#2", "Zeta[2, Indeterminate] === Indeterminate");
    audit_true("Zeta#2", "Zeta[Indeterminate, 2] === Indeterminate");
    /* AppellF2#1: AppellF2 with an Indeterminate argument */
    audit_true("AppellF2#1", "AppellF2[1, 1, 1, 2, 2, 0.3, Indeterminate] === Indeterminate");
    /* ArcSin#2: ArcSin[Indeterminate] */
    audit_true("ArcSin#2", "ArcSin[Indeterminate] === Indeterminate");
    /* BernoulliB#4: BernoulliB[n, x] at the infinities is x^n */
    audit_true("BernoulliB#4", "Check[BernoulliB[3, Infinity], \"m\"] === Infinity");
    audit_true("BernoulliB#4", "BernoulliB[3, -Infinity] === -Infinity");
    audit_true("BernoulliB#4", "BernoulliB[2, -Infinity] === Infinity");
    /* BesselI#2: BesselI[-1/2, 0.] is a pole like the exact call */
    audit_true("BesselI#2", "BesselI[-0.5, 0.] === ComplexInfinity");
    audit_true("BesselI#2", "BesselI[-1/2, 0] === ComplexInfinity");
    /* BesselI#3: BesselI at the infinities */
    audit_true("BesselI#3", "BesselI[0, Infinity] === Infinity");
    audit_true("BesselI#3", "BesselI[1, -Infinity] === -Infinity");
    audit_true("BesselI#3", "BesselI[2, -Infinity] === Infinity");
    /* BesselY#2: BesselY[n, Infinity] */
    audit_true("BesselY#2", "BesselY[0, Infinity] === 0");
    /* CubeRoot#5: CubeRoot[Indeterminate] */
    audit_true("CubeRoot#5", "CubeRoot[Indeterminate] === Indeterminate");
    /* ExpIntegralE#5: ExpIntegralE[n, Infinity] */
    audit_true("ExpIntegralE#5", "ExpIntegralE[2, Infinity] === 0");
    /* ExpIntegralE#6: ExpIntegralE[1, Indeterminate] */
    audit_true("ExpIntegralE#6", "ExpIntegralE[1, Indeterminate] === Indeterminate");
    /* GegenbauerC#3: GegenbauerC with an Indeterminate argument */
    audit_true("GegenbauerC#3", "GegenbauerC[Indeterminate, 1, 0.5] === Indeterminate");
    audit_true("GegenbauerC#3", "GegenbauerC[3, 1, Indeterminate] === Indeterminate");
    /* HankelH1#1: HankelH1 at 0 for a rational order */
    audit_true("HankelH1#1", "HankelH1[1/2, 0] === ComplexInfinity");
    audit_true("HankelH1#1", "HankelH1[0.5, 0] === ComplexInfinity");
    /* HankelH1#2: HankelH1[0, Indeterminate] */
    audit_true("HankelH1#2", "HankelH1[0, Indeterminate] === Indeterminate");
    /* HankelH2#1: HankelH2 at 0 for a rational order */
    audit_true("HankelH2#1", "HankelH2[1/2, 0] === ComplexInfinity");
    audit_true("HankelH2#1", "BesselY[1/2, 0] === ComplexInfinity");
    /* HurwitzZeta#2: HurwitzZeta[Infinity, a] term by term */
    audit_true("HurwitzZeta#2", "Check[HurwitzZeta[Infinity, 2], \"m\"] === 0");
    audit_true("HurwitzZeta#2", "HurwitzZeta[Infinity, 1] === 1");
    audit_true("HurwitzZeta#2", "HurwitzZeta[Infinity, 1/2] === Infinity");
    /* Hypergeometric0F1#6: Hypergeometric0F1 with an Indeterminate argument */
    audit_true("Hypergeometric0F1#6", "Hypergeometric0F1[1, Indeterminate] === Indeterminate");
    /* Hypergeometric1F1#4: 1F1 at the infinities */
    audit_true("Hypergeometric1F1#4", "Hypergeometric1F1[1, 2, ComplexInfinity] === Indeterminate");
    audit_true("Hypergeometric1F1#4", "Check[Hypergeometric1F1[1, 2, Infinity], \"m\"] === Infinity");
    audit_true("Hypergeometric1F1#4", "Hypergeometric1F1[-2, 1, Infinity] === Infinity");
    audit_true("Hypergeometric1F1#4", "Hypergeometric1F1[-1, 1, Infinity] === -Infinity");
    /* HypergeometricPFQ#1: 2F1(1, 1; 2; z) at the divergent z = 1.0 agrees with z = 1 */
    audit_true("HypergeometricPFQ#1", "HypergeometricPFQ[{1, 1}, {2}, 1] === Infinity");
    audit_true("HypergeometricPFQ#1", "HypergeometricPFQ[{1, 1}, {2}, 1.] === Infinity");
    audit_true("HypergeometricPFQ#1", "Hypergeometric2F1[1, 1, 2, 1.] === Infinity");
    audit_true("HypergeometricPFQ#1", "Hypergeometric2F1[1, 1, -2, 0.5] === ComplexInfinity");
    audit_true("HypergeometricPFQ#1", "Abs[Hypergeometric2F1[1, 1, 2, 0.5] - 2 Log[2.]] < 10.^-14");
    /* InverseJacobiCN#2: InverseJacobiCN with an Indeterminate argument */
    audit_true("InverseJacobiCN#2", "InverseJacobiCN[Indeterminate, 0.5] === Indeterminate");
    /* InverseJacobiSN#3: InverseJacobiSN with an Indeterminate argument */
    audit_true("InverseJacobiSN#3", "InverseJacobiSN[Indeterminate, 0.3] === Indeterminate");
    /* JacobiCS#2: JacobiCS with an Indeterminate argument */
    audit_true("JacobiCS#2", "JacobiCS[Indeterminate, 0.3] === Indeterminate");
    audit_true("JacobiCS#2", "JacobiCS[0.5, Indeterminate] === Indeterminate");
    /* LogGamma#2: LogGamma at a bignum non-positive integer */
    audit_true("LogGamma#2", "LogGamma[-10^30] === Infinity");
    audit_true("LogGamma#2", "LogGamma[-3] === Infinity");
    /* RiemannR#1: RiemannR[Indeterminate] */
    audit_true("RiemannR#1", "RiemannR[Indeterminate] === Indeterminate");
    /* RiemannSiegelTheta#1: RiemannSiegelTheta[Indeterminate] */
    audit_true("RiemannSiegelTheta#1", "RiemannSiegelTheta[Indeterminate] === Indeterminate");
    /* RiemannSiegelZ#2: no half-evaluated E^(I ComplexInfinity) */
    audit_true("RiemannSiegelZ#2", "FreeQ[RiemannSiegelZ[0.5 I], ComplexInfinity]");
    audit_true("RiemannSiegelZ#2", "FreeQ[RiemannSiegelZ[-0.5 I], ComplexInfinity]");
    /* SphericalBesselJ#5: SphericalBesselJ[0, Indeterminate] */
    audit_true("SphericalBesselJ#5", "SphericalBesselJ[0, Indeterminate] === Indeterminate");
    /* StieltjesGamma#2: StieltjesGamma at an arbitrary-precision pole */
    audit_true("StieltjesGamma#2", "StieltjesGamma[2, -2.`30] === ComplexInfinity");
    audit_true("StieltjesGamma#2", "StieltjesGamma[2, 0``10] === ComplexInfinity");
    audit_true("StieltjesGamma#2", "StieltjesGamma[2, -2.] === ComplexInfinity");
    /* Surd#5: Surd[Indeterminate, 3] */
    audit_true("Surd#5", "Surd[Indeterminate, 3] === Indeterminate");
    /* WaveletPsi#2: WaveletPsi at Indeterminate */
    audit_true("WaveletPsi#2", "WaveletPsi[DaubechiesWavelet[2], Indeterminate] === Indeterminate");
}

/* ================================================================== */
/* Wave 4 G1-specfun: special-functions medium findings.               */
/* ================================================================== */

static void test_wave4_exact_forms(void) {
    audit_true("Erf#1", "Erf[-2] === -Erf[2]");
    audit_true("Erf#1", "Erf[-2] + Erf[2] === 0");
    audit_true("LegendreP#1", "LegendreP[2, -1, x] === x Sqrt[1 - x^2]/2");
    audit_true("LegendreP#1", "LegendreP[2, -1, 1/2] === Sqrt[3]/8");
    audit_true("Hyperfactorial#1", "{Hyperfactorial[-1], Hyperfactorial[-2], Hyperfactorial[-3]} === {1, -1, -4}");
    audit_true("Hyperfactorial#2", "Quiet[Check[{Hyperfactorial[-1.], Hyperfactorial[-2.]}, \"m\"]] === {1., -1.}");
    audit_true("BarnesG#4", "BarnesG[-2^70] === 0");
    audit_true("EulerE#1", "{EulerE[5001], EulerE[10^30 + 1], EulerE[1001, 1/2]} === {0, 0, 0}");
    audit_true("QPochhammer#1", "Abs[QPochhammer[0.5] - 0.288788095086602] < 10^-12");
    audit_true("Sinc#1", "{Sinc[Pi], Sinc[Pi/2], Sinc[Pi/4], Sinc[3 Pi], Sinc[I]} === "
                         "{0, 2/Pi, 2 Sqrt[2]/Pi, 0, Sinh[1]}");
    audit_true("CubeRoot#1", "CubeRoot[1/10^30] === 1/10^10");
    audit_true("CubeRoot#1", "CubeRoot[1/2^300] === 1/2^100");
    audit_true("CubeRoot#2", "D[CubeRoot[x], x] === 1/(3 Surd[x, 3]^2)");
    audit_true("Surd#2", "D[Surd[x, 3], x] === 1/(3 Surd[x, 3]^2)");
    audit_true("Pochhammer#3", "D[Pochhammer[a, n], a] === Pochhammer[a, n] (PolyGamma[0, a + n] - PolyGamma[0, a])");
    audit_true("HeavisideTheta#3", "D[HeavisideTheta[x, y], x] === DiracDelta[x] HeavisideTheta[y]");
    audit_true("Hypergeometric0F1#3", "Quiet[Head[Hypergeometric0F1[{1, 2}, {1., 2., 3.}]]] === Hypergeometric0F1");
    audit_true("Hypergeometric0F1#3", "Quiet[Length[Hypergeometric0F1[{1, 2}, {1., 2., 3.}]]] === 2");
    audit_true("Hypergeometric0F1#3", "Quiet[Check[Hypergeometric0F1[{1, 2}, {1., 2., 3.}], \"m\", Thread::tdlen]] === \"m\"");
    audit_true("StruveL#1", "Quiet[Check[StruveL[-3/2, 0], \"m\"]] === 0");
    audit_true("ChebyshevU#2", "Quiet[Check[{ChebyshevU[2.5, -1.], ChebyshevU[2.5, -1]}, \"m\"]] === "
                               "{ComplexInfinity, ComplexInfinity}");
    audit_true("HurwitzZeta#1", "HurwitzZeta[-1000, -1] === 1");
    audit_true("Beta#2", "Beta[10^20, 2] === 1/10000000000000000000100000000000000000000");
    audit_true("Beta#3", "Quiet[Check[Beta[-2., 1.], \"m\"]] === -0.5");
    audit_true("Beta#4", "Beta[-2, 1.] === -0.5 && Beta[2, -5.] === 0.05");
    audit_true("FunctionExpand#1", "FunctionExpand[Log[1 + I]] === Log[2]/2 + I Pi/4");
    audit_true("FunctionExpand#1", "FunctionExpand[Log[1 - I]] === Log[2]/2 - I Pi/4");
    audit_true("FunctionExpand#2", "FunctionExpand[Beta[a, b]] === Gamma[a] Gamma[b]/Gamma[a + b]");
    audit_true("FunctionExpand#2", "FunctionExpand[Binomial[n, k]] === Gamma[1 + n]/(Gamma[1 + k] Gamma[1 - k + n])");
    audit_true("FunctionExpand#2", "FunctionExpand[Pochhammer[a, n]] === Gamma[a + n]/Gamma[a]");
    audit_true("FunctionExpand#2", "FunctionExpand[n!] === Gamma[1 + n]");
    audit_true("SphericalBesselJ#3", "Normal[Series[SphericalBesselJ[0, x], {x, 0, 4}]] === 1 - x^2/6 + x^4/120");
    audit_true("HeavisideTheta#1", "Simplify[HeavisideTheta[x], x > 0] === 1");
    audit_true("HeavisideTheta#1", "Refine[HeavisideTheta[x], x < 0] === 0");
    audit_true("HeavisideTheta#1", "Refine[HeavisideTheta[x - 1], x > 2] === 1");
    audit_true("HeavisideTheta#2", "Limit[HeavisideTheta[x], x -> Infinity] === 1");
    audit_true("HeavisideTheta#2", "Limit[HeavisideTheta[x - 1], x -> 2] === 1");
    audit_true("HeavisideTheta#2", "Limit[HeavisideTheta[x], x -> 0, Direction -> -1] === 1");
    audit_true("MeijerG#1", "Quiet[Check[{MeijerG[{{}, {}}, {{0}, {}}, Infinity], MeijerG[{{}, {}}, {{1}, {}}, Infinity]}, \"m\"]] === {0, 0}");
}

static void test_wave4_numerics(void) {
    audit_true("Zeta#1", "Abs[Zeta[2. + 1000. I] - (0.9532621843464252 - 0.1107231074605998 I)] < 2*10^-15");
    audit_true("Zeta#1", "Abs[Zeta[0.5 + 100. I] - (2.692619885681324 - 0.02038602960259816 I)] < 4*10^-15");
    audit_true("LerchPhi#3", "LerchPhi[0.5, 0., -2.] === 1.75");
    audit_true("LerchPhi#4", "Quiet[Check[{LerchPhi[-1, 1., 0.5], LerchPhi[-1., 1., 0.5]}, \"m\"]] === "
                             "{1.5707963267948966, 1.5707963267948966}");
    audit_true("SphericalBesselJ#1", "SphericalBesselJ[0, 0.`20] == 1 && Precision[SphericalBesselJ[0, 0.`20]] > 15");
    audit_true("SphericalBesselJ#1", "SphericalBesselJ[2, 0.`20] == 0 && !MachineNumberQ[SphericalBesselJ[2, 0.`20]]");
    audit_true("SinIntegral#1", "Re[SinIntegral[100. I]] == 0 && Abs[Im[SinIntegral[100. I]]/1.3577763724269401*^41 - 1] < 10^-14");
    audit_true("SinIntegral#1", "Re[SinIntegral[100`20 I]] == 0 && SinIntegral[-100`20 I] == -SinIntegral[100`20 I]");
    audit_true("InverseJacobiCN#1", "Abs[InverseJacobiCN[-0.3, 1.] - (1.8738202425274146 + 3.141592653589793 I)] < 10^-14");
    audit_true("InverseJacobiCN#1", "InverseJacobiCN[0., 1.] === ComplexInfinity && FreeQ[InverseJacobiCN[-1., 1.], EllipticF]");
    audit_true("InverseJacobiSN#2", "{InverseJacobiSN[1.0, 1.0], InverseJacobiSN[1`30, 1`30], InverseJacobiSN[1, 1.0]} === "
                                    "{ComplexInfinity, ComplexInfinity, ComplexInfinity}");
    audit_true("MathieuCharacteristicB#1", "MachineNumberQ[MathieuCharacteristicB[Pi, 1.]] && "
                                           "MathieuCharacteristicB[Pi, 1.] == MathieuCharacteristicA[Pi, 1.]");
    audit_true("MathieuCharacteristicB#1", "MachineNumberQ[MathieuCharacteristicB[Sqrt[2], 1.]]");
    audit_true("WeierstrassHalfPeriods#1", "Abs[WeierstrassHalfPeriods[{3.`20, 1}][[1]] - N[Pi/Sqrt[6], 20]] < 10^-18 && "
                                           "WeierstrassHalfPeriods[{3.`20, 1}][[2]] === DirectedInfinity[I]");
    audit_true("WeierstrassHalfPeriods#1", "Precision[WeierstrassHalfPeriods[{12.`30, -8}][[1]]] > 25");
    audit_true("WeierstrassSigma#1", "Abs[WeierstrassSigma[0.3`20, {3, 1}] - N[WeierstrassSigma[3/10, {3, 1}], 20]] < 10^-18");
    audit_true("WeierstrassSigma#1", "NumberQ[WeierstrassSigma[0.3`20, {3, -1}]]");
    audit_true("WeierstrassZeta#1", "Abs[WeierstrassZeta[0.3`20, {3, 1}] - N[WeierstrassZeta[3/10, {3, 1}], 20]] < 10^-18");
    audit_true("WeierstrassZeta#1", "NumberQ[WeierstrassZeta[0.3`20, {3, -1}]]");
    audit_true("FresnelC#1", "Precision[FresnelC[1.5`10]] == 10. && Precision[FresnelC[1.5`5]] == 5.");
    audit_true("FresnelC#1", "Precision[FresnelC[1.5`10 + I]] == 10.");
    audit_true("FresnelS#1", "Precision[FresnelS[1.5`10]] == 10. && Abs[FresnelS[1.5`10] - 0.6975049601] < 10^-9");
    audit_true("FresnelS#1", "Precision[FresnelS[1.5`10 + 0.5`10 I]] == 10.");
    audit_true("StieltjesGamma#1", "Precision[StieltjesGamma[1, 1.5`10]] == 10. && Precision[StieltjesGamma[1, 1.5`2]] == 2.");
    audit_true("HaarWavelet#2", "{WaveletPhi[HaarWavelet[], 1 - 10^-30], WaveletPsi[HaarWavelet[], 1/2 - 10^-30], "
                                "WaveletPsi[HaarWavelet[], 1/2 + 10^-30]} == {1, 1, -1}");
    audit_true("WaveletPsi#1", "{WaveletPsi[HaarWavelet[], 1/2 + 10^-30], WaveletPsi[HaarWavelet[], 1 - 10^-30]} == {-1, -1}");
    audit_true("ZetaZero#1", "Abs[N[ZetaZero[1, 20]] - (0.5 + 21.022039638771556 I)] < 10^-12");
    audit_true("AiryBiPrimeZero#1", "Abs[N[AiryBiPrimeZero[1, -3]] + 4.07316] < 10^-5");
    audit_true("BesselJZero#2", "Abs[N[BesselJZero[0, 1, 10]] - 11.791534439014281] < 10^-12");
    audit_true("BesselJZero#1", "Abs[BesselJZero[1000., 1.] - 1018.66] < 0.01 && MachineNumberQ[BesselJZero[700., 1.]]");
    audit_true("BesselYZero#1", "708 < N[BesselYZero[700, 1]] < 709");
    audit_true("AiryAiPrimeZero#1", "MachineNumberQ[N[AiryAiPrimeZero[2^62]]] && N[AiryAiPrimeZero[2^62]] < -10^12");
    audit_true("AiryBiZero#1", "MachineNumberQ[N[AiryBiZero[2^62]]] && N[AiryBiZero[2^62]] < -10^12");
    audit_true("BarnesG#3", "Abs[BarnesG[50.]/N[BarnesG[50], 30] - 1] < 10^-15");
    audit_true("BarnesG#3", "Abs[BarnesG[200.]/N[BarnesG[200], 30] - 1] < 10^-15");
    audit_true("RiemannSiegelZ#1", "Head[RiemannSiegelZ[1.*^50]] === RiemannSiegelZ && "
                                   "Head[RiemannSiegelZ[1.*^300]] === RiemannSiegelZ");
    audit_true("RiemannSiegelZ#1", "Abs[RiemannSiegelZ[100.] - Abs[Zeta[0.5 + 100. I]]] < 10^-12");
    audit_true("DiracDelta#3", "Integrate[x^2 Derivative[1][DiracDelta][x - 1], {x, -Infinity, Infinity}] === -2");
    audit_true("DiracDelta#3", "Integrate[x Derivative[1][DiracDelta][x], {x, -1, 1}] === -1");
    audit_true("DaubechiesWavelet#1", "With[{h = WaveletFilterCoefficients[DaubechiesWavelet[24], \"PrimalLowpass\"]}, "
                                      "Length[h] == 48 && Abs[Total[h[[All, 2]]] - 1] < 10^-12]");
    audit_true("DaubechiesWavelet#1", "MachineNumberQ[WaveletPhi[DaubechiesWavelet[30], 3.3]]");
    audit_true("KelvinKer#1", "Abs[KelvinKer[2, 1.*^-8] - 0.5] < 10^-15");
    audit_true("Hypergeometric2F1#2", "Abs[Hypergeometric2F1[1 + I, 2, 3, 2.] - "
                                      "Hypergeometric2F1[1 + I, 2, 3, 2. - 1.*^-13 I]] < 10^-10");
    audit_true("Hypergeometric2F1#2", "Abs[Hypergeometric2F1[1/2 + I/10, 1/2, 1, 2.] - "
                                      "Hypergeometric2F1[0.5 + 0.1 I, 0.5, 1, 2. - 1.*^-13 I]] < 10^-10");
    audit_true("JacobiAmplitude#1", "Abs[JacobiAmplitude[1.*^20, 0.5]/N[Pi 10^20/(2 EllipticK[1/2]), 20] - 1] < 10^-15");
    audit_true("AiryAiPrime#1", "Abs[AiryAiPrime[-10.^30]/N[AiryAiPrime[SetPrecision[-10.^30, 80]], 20] - 1] < 10^-12");
    audit_true("AiryAiPrime#1", "Abs[AiryAiPrime[-1.*^12]/N[AiryAiPrime[SetPrecision[-1.*^12, 60]], 20] - 1] < 10^-12");
    audit_true("AiryAiPrime#1", "Abs[AiryBi[-10.^30]/N[AiryBi[SetPrecision[-10.^30, 80]], 20] - 1] < 10^-12");
    audit_true("HypergeometricU#3", "Abs[HypergeometricU[-200, 1, 1.]/N[HypergeometricU[-200, 1, 1], 30] - 1] < 10^-15");
    audit_true("JacobiCS#1", "Abs[JacobiCS[1.*^6, 3.*^-17]/N[JacobiCS[10^6, 3*10^-17], 30] - 1] < 10^-15");
    audit_true("JacobiCS#1", "Abs[JacobiCS[5.*^9, 10.^-20]/N[JacobiCS[5*10^9, 10^-20], 30] - 1] < 10^-15");
}

/* ================================================================== */
/* Wave 5 G1-specfun: special-functions low findings.                  */
/* ================================================================== */

static void test_wave5_messages(void) {
    audit_chk("FunctionExpand#4", "Quiet[Check[FunctionExpand[], fired, FunctionExpand::argt]]", "fired");
    /* FunctionExpand takes options, so a third positional argument is
     * reported as a non-option (General::nonopt), as for every such head */
    audit_chk("FunctionExpand#4", "Quiet[Check[FunctionExpand[1, 2, 3], fired, FunctionExpand::nonopt]]", "fired");
    audit_chk("FunctionExpand#4", "Quiet[FunctionExpand[1, 2, 3]]", "FunctionExpand[1, 2, 3]");
    audit_chk("KelvinBer#3", "Quiet[Check[KelvinBer[], fired, KelvinBer::argt]]", "fired");
    audit_chk("KelvinBer#3", "Quiet[Check[KelvinBer[1, 2, 3], fired, KelvinBer::argt]]", "fired");
    audit_chk("StieltjesGamma#3", "Quiet[Check[StieltjesGamma[], fired, StieltjesGamma::argt]]", "fired");
    audit_chk("StieltjesGamma#3", "Quiet[Check[StieltjesGamma[1, 2, 3], fired, StieltjesGamma::argt]]", "fired");
    audit_chk("WeberE#2", "Quiet[Check[WeberE[], fired, WeberE::argt]]", "fired");
    audit_chk("WeberE#2", "Quiet[Check[WeberE[1], fired, WeberE::argtu]]", "fired");
    audit_true("HeavisideTheta#4", "MemberQ[Attributes[HeavisideTheta], NumericFunction]");
    audit_true("HeavisideTheta#5", "StringQ[HeavisideTheta::usage] && "
                                   "StringContainsQ[HeavisideTheta::usage, \"Heaviside theta function\"]");
    audit_true("SphericalBesselJ#6", "StringContainsQ[Information[SphericalBesselJ], \"Sqrt[Pi]/Sqrt[2 z]\"]");
}

static void test_wave5_exact(void) {
    audit_true("Erf#3", "Erf[InverseErf[x]] === x");
    audit_true("Erfc#2", "Erfc[InverseErfc[x]] === x");
    audit_true("Zeta#3", "Zeta[-10^20] === 0 && Zeta[-1000] === 0");
    audit_true("BernoulliB#2", "{BernoulliB[5001], BernoulliB[10^30 + 1]} === {0, 0} && BernoulliB[5001.] === 0.");
    audit_true("BernoulliB#3", "D[BernoulliB[n, x], x] === n BernoulliB[-1 + n, x]");
    audit_true("EulerE#3", "D[EulerE[n, x], x] === n EulerE[-1 + n, x]");
    audit_true("FresnelS#2", "FresnelS[-2] === -FresnelS[2] && FresnelS[-1/2] === -FresnelS[1/2] && FresnelS[-Pi] === -FresnelS[Pi]");
    audit_true("InverseErf#3", "InverseErf[-1/2] === -InverseErf[1/2]");
    audit_true("SinIntegral#2", "SinIntegral[-2] === -SinIntegral[2] && SinIntegral[-1/2] === -SinIntegral[1/2]");
    audit_true("SinIntegral#2", "SinIntegral[-2] + SinIntegral[2] === 0");
    audit_true("SinhIntegral#2", "SinhIntegral[I x] === I SinIntegral[x] && SinhIntegral[-I] === -I SinIntegral[1]");
    audit_true("SinhIntegral#2", "SinhIntegral[-I x] === -I SinIntegral[x] && SinIntegral[I x] === I SinhIntegral[x]");
    audit_true("LerchPhi#5", "LerchPhi[-1, 1, 1/2] === Pi/2");
    audit_true("QPochhammer#3", "QPochhammer[0, q] === 1");
    audit_true("Surd#4", "Surd[0, 10^20] === 0 && 0^(1/10^20) === 0");
    audit_true("ExpIntegralE#4", "ExpIntegralE[5/2, 0] === 2/3 && ExpIntegralE[1/2, 0] === ComplexInfinity");
    audit_true("ExpIntegralE#7", "ExpIntegralE[-1, z] === E^-z (1 + z)/z^2");
    audit_true("ExpIntegralE#7", "Abs[N[ExpIntegralE[-2, 3]] - ExpIntegralE[-2, 3.]] < 10^-15");
    audit_true("WhittakerM#2", "WhittakerM[1., 1/2, 0] === 0. && Precision[WhittakerM[1., 1/2, 0]] === MachinePrecision");
    audit_true("WhittakerM#3", "Quiet[Check[WhittakerM[1, -3/4, 0.], \"m\"]] === ComplexInfinity");
    audit_true("WhittakerM#3", "WhittakerM[1, -3/4, 0] === ComplexInfinity");
}

static void test_wave5_numerics(void) {
    audit_true("PolyLog#1", "{PolyLog[-3, 2.], PolyLog[-2, 3.], PolyLog[-1, 2.]} === {26., -1.5, 2.}");
    audit_true("PolyLog#1", "Abs[PolyLog[-3, -2.] - 2/27] < 10^-15 && Abs[PolyLog[-2, 3. + I] - N[PolyLog[-2, 3 + I]]] < 10^-14");
    /* sub-machine precision is kept */
    audit_true("CosIntegral#1", "Precision[CosIntegral[1.5`3]] == 3. && Abs[CosIntegral[1.5`3] - 0.4703563] < 10^-2");
    audit_true("CoshIntegral#1", "Precision[CoshIntegral[2.`5]] == 5. && Precision[CoshIntegral[2.`16]] == 16.");
    audit_true("Erfi#1", "Precision[Erfi[1.5`3]] == 3. && Precision[Erf[1.5`3]] == 3.");
    audit_true("HurwitzZeta#3", "Precision[HurwitzZeta[3`5, 2]] == 5.");
    audit_true("LogIntegral#1", "Precision[LogIntegral[2.`10]] == 10. && Abs[LogIntegral[2.`10] - 1.045163780] < 10^-8");
    audit_true("LogIntegral#1", "Precision[LogIntegral[2.`2]] == 2.");
    audit_true("MathieuCharacteristicB#3", "Precision[MathieuCharacteristicB[2, 1.`5]] == 5.");
    audit_true("SinhIntegral#1", "Precision[SinhIntegral[1.5`3]] == 3. && MachineNumberQ[SinhIntegral[1.5]]");
    audit_true("Sinc#3", "Precision[Sinc[1.`30 + I]] == 30.");
    audit_true("HarmonicNumber#2", "HarmonicNumber[0., 2] === 0.");
    audit_true("KelvinBei#1", "KelvinBei[-2, 0.] === 0. && Precision[KelvinBei[-2, 0.`20]] === Precision[KelvinBei[2, 0.`20]]");
    audit_true("KelvinKei#1", "Abs[KelvinKei[0, 0.`20] + N[Pi/4, 20]] < 10^-19 && Precision[KelvinKei[0, 0.`20]] == 20.");
    audit_true("KelvinKei#1", "KelvinKei[1, 0.] === ComplexInfinity");
    audit_true("StruveL#3", "StruveL[-1.5, 0] === 0.");
    audit_true("SphericalBesselJ#4", "SphericalBesselJ[0.5, 0.] === 0.");
    audit_true("Hypergeometric1F1#3", "Hypergeometric1F1[1, b, 0.] === 1.");
    audit_true("QPochhammer#5", "With[{v = QPochhammer[0.5, 0.99999]}, NumberQ[v] && -25288 < Log[10, v] < -25286]");
    audit_true("QPochhammer#5", "With[{v = QPochhammer[0.5, -0.99999]}, NumberQ[v] && -2908 < Log[10, v] < -2905]");
    audit_true("QPochhammer#5", "Abs[QPochhammer[0.5`20, 0.999`20]/"
                                "Exp[Total[Log[1 - 0.5`20 0.999`20^Range[0, 50000]]]] - 1] < 10^-17");
    audit_true("QPochhammer#5", "Abs[QPochhammer[0.5, 0.9995]/"
                                "QPochhammer[SetPrecision[0.5, 30], SetPrecision[0.9995, 30]] - 1] < 10^-14");
    audit_true("LogGamma#3", "Abs[LogGamma[1.*^307]/N[LogGamma[10^307], 20] - 1] < 10^-15");
    audit_true("ProductLog#1", "With[{w = ProductLog[10^20, 1.]}, Abs[Re[w] + 47.8896] < 10^-3 && Abs[Im[w]/(2 Pi 10^20) - 1] < 10^-15]");
    audit_true("ProductLog#1", "Im[ProductLog[2^62, 1.]] > 0 && Im[ProductLog[-10^20, 1.]] < 0");
    audit_true("ProductLog#1", "With[{w = ProductLog[10^20, 1.`30]}, Abs[w Exp[w] - 1] < 10^-10]");
}

static void test_wave5_more(void) {
    audit_true("ArcSin#3", "ArcSin[Sin[Pi/12]] === Pi/12 && ArcSin[Sin[Pi/10]] === Pi/10");
    audit_true("ArcSin#3", "ArcSin[(1 + Sqrt[5])/4] === 3 Pi/10 && ArcCos[Cos[Pi/12]] === Pi/12");
    audit_true("LegendreQ#6", "Quiet[Check[LegendreQ[10^20, 0.3], \"m\"]] =!= \"m\"");
    audit_true("Hypergeometric0F1#4", "Quiet[Hypergeometric0F1[2, 1.*^300]] === Overflow[]");
    audit_true("Hypergeometric0F1#4", "Quiet[Check[Hypergeometric0F1[2, 1.*^300], \"m\", General::ovfl]] === \"m\"");
    audit_true("AppellF1#1", "Quiet[Check[AppellF1[-2, 1, 1, -1, x, y], \"m\"]] === AppellF1[-2, 1, 1, -1, x, y]");
    audit_true("AppellF1#2", "Head[AppellF1[-1, 1, 1, -2, 3/10, 1/5]] === Head[AppellF1[-1, 1, 1, -2, 0.3, 0.2]] === AppellF1");
    audit_true("AppellF1#2", "AppellF1[-2, 1, 1, 3, x, y] === Expand[AppellF1[-2, 1, 1, 3, x, y]]");
    audit_true("WeierstrassInvariants#2", "Quiet[Check[WeierstrassInvariants[{0, I}], \"m\", Power::infy]] === WeierstrassInvariants[{0, I}]");
    audit_true("WeierstrassInvariants#2", "Quiet[Check[WeierstrassInvariants[{0, I}], \"m\", EllipticReducedHalfPeriods::nind]] === \"m\"");
    audit_true("AiryBi#1", "Abs[AiryBi[-10.^8]/N[AiryBi[-10^8], 30] - 1] < 2*10^-15");
    audit_true("AiryBiPrime#1", "Abs[AiryBiPrime[-1.*^7]/N[AiryBiPrime[-10^7], 30] - 1] < 2*10^-15 && "
                                "Abs[AiryBiPrime[-1.*^8]/N[AiryBiPrime[-10^8], 30] - 1] < 2*10^-15");
    audit_true("BarnesG#5", "BarnesG[5.] === 12. && BarnesG[1.] === 1.");
    audit_true("BarnesG#5", "Abs[BarnesG[-0.5]/N[BarnesG[-1/2], 30] - 1] < 10^-15");
    audit_true("ChebyshevT#2", "Abs[ChebyshevT[2000001, 0.3] - N[Cos[2000001 ArcCos[SetPrecision[0.3, 60]]], 20]] < 10^-14");
    audit_true("ChebyshevT#2", "Precision[ChebyshevT[2000001, 0.3`30]] == 30.");
    audit_true("ChebyshevT#3", "ChebyshevT[2001, 1/2] === -1 && ChebyshevT[2001, 2] == ChebyshevT[2001, 2.`3000]");
    audit_true("ChebyshevT#3", "ChebyshevT[2003, 3/7] === 2 (3/7) ChebyshevT[2002, 3/7] - ChebyshevT[2001, 3/7]");
    audit_true("ChebyshevT#4", "Abs[ChebyshevT[10^30, 0.5`50] + 1/2] < 10^-45");
    audit_true("ChebyshevU#3", "ChebyshevU[2001, 1/2] === -1 && IntegerQ[ChebyshevU[2001, 2]]");
    audit_true("ChebyshevU#3", "ChebyshevU[2003, 3/7] === 2 (3/7) ChebyshevU[2002, 3/7] - ChebyshevU[2001, 3/7]");
    audit_true("ChebyshevU#4", "Abs[ChebyshevU[10^7, 0.5] + 1] < 10^-8 && MachineNumberQ[ChebyshevU[10^7, 0.5]]");
    audit_true("GegenbauerC#2", "Abs[GegenbauerC[10^6, 1.5, 0.3]/N[GegenbauerC[10^6, 3/2, 3/10], 20] - 1] < 10^-12");
    audit_true("EulerE#2", "IntegerQ[EulerE[6000]] && EulerE[6000] != 0 && Abs[EulerE[6000]/N[EulerE[6000.]] - 1] < 10^-14");
    audit_true("BesselJZero#3", "Abs[BesselJZero[1/2, 1.] - Pi] < 10^-14 && BesselJZero[1/2, 1.] == BesselJZero[0.5, 1.]");
    audit_true("BesselJZero#4", "Abs[BesselJZero[0, 1000001.] - 3141595.009784] < 10^-5");
    audit_true("BesselJZero#4", "Abs[BesselJ[0, BesselJZero[0, 1000001.]]] < 10^-12");
    audit_true("BesselYZero#2", "Abs[N[BesselYZero[0, 10^6 + 1]] - 3141593.43898896] < 10^-6");
    audit_true("FunctionExpand#3", "FunctionExpand[HarmonicNumber[n]] === EulerGamma + PolyGamma[0, 1 + n]");
    audit_true("FunctionExpand#3", "FunctionExpand[ChebyshevT[n, x]] === Cos[n ArcCos[x]]");
    audit_true("FunctionExpand#3", "FunctionExpand[GoldenRatio] === (1 + Sqrt[5])/2");
    audit_true("CubeRoot#3", "Integrate[CubeRoot[x], x] === 3/4 x Surd[x, 3]");
    audit_true("CubeRoot#3", "Integrate[x^2 Surd[x, 5], x] === 5/16 x^3 Surd[x, 5]");
    audit_true("CubeRoot#4", "CubeRoot[x]^3 === x && Simplify[CubeRoot[x]^3] === x && Surd[x, 3]^6 === x^2");
    audit_true("HeavisideTheta#6", "With[{s = Series[HeavisideTheta[x], {x, 1, 2}]}, "
                                   "Normal[s] === 1 && s[[5]] === 3 && FreeQ[s, DiracDelta]]");
    audit_true("WeierstrassPPrime#1", "Quiet[Check[Normal[Series[WeierstrassPPrime[u, {g2, g3}], {u, 0, 3}]], \"m\"]] === "
                                      "-2/u^3 + g2 u/10 + g3 u^3/7");
    audit_true("WeierstrassPPrime#1", "Abs[(Normal[Series[WeierstrassP[u, {2., 1.}], {u, 0, 10}]] /. u -> 0.3) - "
                                      "WeierstrassP[0.3, {2., 1.}]] < 10^-8");
    audit_true("HankelH2#2", "Series[HankelH2[0, x], {x, 0, 1}] === Series[BesselJ[0, x] - I BesselY[0, x], {x, 0, 1}] && "
                             "FreeQ[Series[HankelH2[0, x], {x, 0, 1}], HankelH2]");
    audit_true("HankelH2#2", "Abs[(Normal[Series[HankelH2[0, x], {x, 0, 8}]] /. x -> 0.3) - HankelH2[0, 0.3]] < 10^-9");
    audit_true("KelvinKei#2", "Simplify[Normal[Series[KelvinKei[x], {x, 0, 2}]] - "
                              "(-Pi/4 + (1/4 - EulerGamma/4 + Log[2]/4 - Log[x]/4) x^2)] === 0");
    audit_true("KelvinKei#2", "Abs[(Normal[Series[KelvinKei[x], {x, 0, 12}]] /. x -> 0.5) - KelvinKei[0.5]] < 10^-12 && "
                              "Abs[(Normal[Series[KelvinKer[x], {x, 0, 12}]] /. x -> 0.5) - KelvinKer[0.5]] < 10^-12");
    audit_true("ExpIntegralE#8", "Normal[Series[ExpIntegralE[2, x], {x, 0, 2}]] === "
                                 "1 + (-1 + EulerGamma + Log[x]) x - x^2/2");
    audit_true("ExpIntegralE#8", "Abs[(Normal[Series[ExpIntegralE[2, x], {x, 0, 4}]] /. x -> 0.01) - "
                                 "ExpIntegralE[2, 0.01]] < 10^-12");
    audit_true("SphericalBesselY#2", "Normal[Series[SphericalBesselY[0, x], {x, 0, 2}]] === -1/x + x/2");
    audit_true("SphericalBesselY#2", "Normal[Series[SphericalBesselY[1, x], {x, 0, 1}]] === -1/x^2 - 1/2");
    audit_true("HermiteH#2", "HermiteH[2^70 + 1, 0] === 0 && HermiteH[7, 0] === 0");
    audit_true("HermiteH#3", "Precision[HermiteH[200, 100.]] === MachinePrecision && "
                             "Abs[HermiteH[200, 100.]/N[HermiteH[200, 100], 30] - 1] < 10^-14");
    audit_true("ParabolicCylinderD#2", "NumberQ[TimeConstrained[ParabolicCylinderD[10.^6, 1.], 20]]");
    audit_true("ParabolicCylinderD#2", "Abs[Log10[Abs[N[ParabolicCylinderD[10^6, 1], 20]]] - 2782852.6594] < 10^-3");
    audit_true("ParabolicCylinderD#2", "Abs[N[ParabolicCylinderD[2000, 3/2], 20]/ParabolicCylinderD[2000., 1.5] - 1] < 10^-13");
    audit_true("BarnesG#6", "Abs[N[BarnesG[13/2], 30] - 2548.74576956849898973590610636`30] < 10^-25");
    audit_true("Hypergeometric2F1#3", "Hypergeometric2F1[1, 2, 3, 1] === Infinity && Hypergeometric2F1[1, 2, 3, 1.] === Infinity");
    audit_true("Hypergeometric2F1#3", "Hypergeometric2F1[1, 1, 2, 1] === Infinity && Hypergeometric2F1[1, 1, -2, 0.5] === ComplexInfinity");
    audit_true("InverseJacobiDN#1", "Quiet[Check[InverseJacobiDN[1.0, 0], \"m\"]] === 0. && InverseJacobiDN[1, 0] === 0");
    audit_true("InverseJacobiDN#1", "InverseJacobiDN[1., 0.5] === 0.");
    audit_true("MathieuCharacteristicB#2", "Abs[MathieuCharacteristicB[2.00000000000000000001`30, 1.`30] - "
                                           "MathieuCharacteristicA[2, 1.`30]] < 10^-18");
    audit_true("MathieuCharacteristicB#2", "Abs[MathieuCharacteristicB[1.99999999999999999999`30, 1.`30] - "
                                           "MathieuCharacteristicB[2, 1.`30]] < 10^-18");
    audit_true("WeierstrassPPrime#2", "Integrate[WeierstrassPPrime[u, {g2, g3}], u] === WeierstrassP[u, {g2, g3}]");
    audit_true("WeierstrassPPrime#2", "Integrate[3 WeierstrassPPrime[2 u + 1, {g2, g3}], u] === 3/2 WeierstrassP[1 + 2 u, {g2, g3}]");
    audit_true("HeavisideTheta#6", "Derivative[1][DiracDelta][1] === 0 && Head[Derivative[1][DiracDelta][0]] === Derivative[1][DiracDelta]");
}

/* Wave 6 G2: real type-3 Legendre values, named hypergeometric heads kept,
 * accuracy near singular points and for large parameters. */
static void test_wave6_g2(void) {
    audit_true("LegendreQ#1", "Head[LegendreQ[1/2, 0, 3, 2.]] === Real");
    audit_true("LegendreQ#1", "Abs[LegendreQ[1/2, 0, 3, 2.] - 0.22401429283641563] < 10^-15");
    audit_true("LegendreQ#1", "Head[LegendreQ[1/2, 0, 3, 2.`30]] === Real");
    audit_true("LegendreQ#1", "Head[LegendreQ[0, 0, 3, 2.]] === Real && Abs[LegendreQ[0, 0, 3, 2.] - ArcCoth[2.]] < 10^-15");
    audit_true("LegendreQ#1", "LegendreQ[0, 0, 3, x] === Log[1 + x]/2 - Log[-1 + x]/2");
    audit_chk("Hypergeometric0F1#2", "Hypergeometric0F1[1, 1]", "Hypergeometric0F1[1, 1]");
    audit_chk("Hypergeometric0F1#2", "Hypergeometric0F1[b, z]", "Hypergeometric0F1[b, z]");
    audit_true("Hypergeometric0F1#2", "D[Hypergeometric0F1[b, z], z] === Hypergeometric0F1[1 + b, z]/b");
    audit_chk("Hypergeometric1F1#1", "Hypergeometric1F1[a, b, z]", "Hypergeometric1F1[a, b, z]");
    audit_chk("Hypergeometric1F1#1", "Head[Hypergeometric1F1[1/3, 1/4, z]]", "Hypergeometric1F1");
    audit_true("Hypergeometric1F1#1", "Hypergeometric1F1[1, 2, z] === (-1 + E^z)/z");
    audit_true("Hypergeometric1F1#1", "Abs[Hypergeometric1F1[1/3, 1/4, 2.] - N[Hypergeometric1F1[1/3, 1/4, 2], 20]] < 10^-13");
    audit_chk("Hypergeometric2F1#1", "Hypergeometric2F1[1/2, 1, 2, 1/2]", "Hypergeometric2F1[1/2, 1, 2, 1/2]");
    audit_true("Hypergeometric2F1#1", "D[Hypergeometric2F1[a, b, c, z], z] === a b Hypergeometric2F1[1 + a, 1 + b, 1 + c, z]/c");
    audit_true("Hypergeometric2F1#1", "Simplify[Hypergeometric2F1[a, b, c, z]] === Hypergeometric2F1[a, b, c, z]");
    audit_true("Hypergeometric2F1#1", "Hypergeometric2F1[1, 1, 2, z] === -Log[1 - z]/z");
    audit_true("GegenbauerC#1", "Abs[GegenbauerC[2.5, 1.5, -0.999999]/(-318339.06387664295) - 1] < 4*10^-16");
    audit_true("GegenbauerC#1", "Abs[GegenbauerC[2.5, 1.5, -0.99999]/(-31854.393615169893) - 1] < 4*10^-16");
    audit_true("BetaRegularized#3", "Abs[BetaRegularized[0.5, 1.*^6, 1.*^6] - 0.5] < 10^-15");
    audit_true("BetaRegularized#3", "Abs[BetaRegularized[0.5, 100000., 100000.] - 0.5] < 10^-15");
    audit_true("BarnesG#6", "Abs[N[BarnesG[13/2], 30] - 254874576956849898973590610636/10^26] < 10^-25");
}

int main(void) {
    audit_init();
    TEST(test_legendreq_type3);
    TEST(test_legendreq_large_degree);
    TEST(test_hypergeometric0f1_large_negative);
    TEST(test_hypergeometricu_zero);
    TEST(test_whittakerw_zero);
    TEST(test_webere_zero);
    TEST(test_mathieu_characteristic);
    TEST(test_betaregularized_large_b);
    TEST(test_expintegrale_large_order);
    TEST(test_struveh_large_z);
    TEST(test_mathieu_far_imaginary);
    TEST(test_bernoullib_large_index);
    TEST(test_qpochhammer_large_n);
    /* group A2 */
    TEST(test_a2_gamma_incomplete_at_zero);
    TEST(test_a2_loggamma_small_complex);
    TEST(test_a2_pochhammer_mpfr_pole);
    TEST(test_a2_harmonicnumber_infinity);
    TEST(test_a2_lerchphi_negative_a);
    TEST(test_a2_inverse_erf_tails);
    TEST(test_a2_airyai_large_negative);
    TEST(test_a2_diracdelta_integrals);
    TEST(test_a2_haar_norms);
    TEST(test_a2_jacobip_degenerate);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_exact_forms);
    TEST(test_wave4_numerics);
    TEST(test_wave5_messages);
    TEST(test_wave5_exact);
    TEST(test_wave5_numerics);
    TEST(test_wave5_more);
    TEST(test_wave6_g2);
    return audit_finish("audit_special_functions");
}
