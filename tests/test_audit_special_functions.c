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

int main(void) {
    audit_init();
    TEST(test_betaregularized_large_b);
    TEST(test_expintegrale_large_order);
    TEST(test_struveh_large_z);
    TEST(test_mathieu_far_imaginary);
    TEST(test_bernoullib_large_index);
    TEST(test_qpochhammer_large_n);
    return audit_finish("audit_special_functions");
}
