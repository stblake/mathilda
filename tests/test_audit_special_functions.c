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
    return audit_finish("audit_special_functions");
}
