/* test_adversarial_regressions.c -- regressions from the adversarial reviews of
 * the Bronstein, Trager and ParallelMixedTower integrators (2026-11-14).  Every
 * case was a false non-elementarity certificate, a wrong decline or a leaked
 * message; the open findings are in tests/bronstein_paper/open_gaps.m.
 *
 *  - ParallelMixedTower, cause 1: a residue at Root-valued roots read as
 *    non-constant by syntax (pmtAlgebraicZeroQ, pmtDerivNumeratorZeroQ,
 *    pmtDependsModP: Together does not combine fractions with Root objects);
 *  - ParallelMixedTower, cause 2: dependent tower generators; certificates now
 *    need independent confirmation (pmtConfirmCertificate);
 *  - Risch`ElementaryIntegralQ: a transcendental constant (Log[2]) beside the
 *    tower made a "no" over Q(i) a false certificate (rt_has_transcendental_constant);
 *    its ::undec message now goes through the Quiet/Check funnel;
 *  - Bronstein: one placeholder for every radical power collapsed a denominator
 *    (a wrong "not algebraic" decline and a leaked Power::infy);
 *  - Trager: the trace certificate over Q(i) used a Q(i)-valued "functional"
 *    (traceCertificate gives no verdict over Q(i)): Chebyshev-elementary
 *    binomials with a Gaussian leading coefficient were certified.
 */
#include "bronstein_test_utils.h"

static void test_pmt_root_residues(void) {
    br_run("Quiet[Integrate`ParallelMixedTower[Log[x], x]]; True");   /* load the package (lazy) */
    BR_FALSE("ParallelMixed`Private`pmtAlgebraicZeroQ[Root[1 + #1 + #1^4 &, 1] + 1, {}]");
    BR_TRUE("ParallelMixed`Private`pmtAlgebraicZeroQ[With[{r = Root[1 + #1 + #1^4 &, 1]}, t (r + r^4 + 1)], {t}]");
    BR_TRUE("With[{r = Root[1 + #1 + #1^4 &, 1]}, ParallelMixed`Private`pmtDerivNumeratorZeroQ["
            "(1 + t + 2 t r + 2 t r^2 + 4 r^3 + 4 t r^3 + 2 t r^5)/(1 + t + 4 r^3 + 4 t r^3), t, {t}]]");
    BR_FALSE("ParallelMixed`CertifyNonconstant[With[{r = Root[1 + #1 + #1^4 &, 1]},"
             " (1 + t + 2 t r + 2 t r^2 + 4 r^3 + 4 t r^3 + 2 t r^5)/(1 + t + 4 r^3 + 4 t r^3)], {t}]");
    BR_FALSE("ParallelMixed`Private`pmtDependsModP[t g (1 + g + g^4), t, 1 + g + g^4, g]");
    BR_TRUE("ParallelMixed`Private`pmtDependsModP[t g, t, 1 + g + g^4, g]");
    const char* elem[] = { "(3 x^2 - 1)/(x^3 - x + 1) + Log[x]", "(3 x^2 - 1)/(x^3 - x + 1) + Sqrt[Log[x]]/x",
                           "Together[D[Log[x^2 + 1 + Sqrt[x]] + Log[E^x + 1], x]]", NULL };
    for (int i = 0; elem[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "With[{F = Quiet[Integrate`ParallelMixedTower[%s, x]]}, !ListQ[F] &&"
                 " PossibleZeroQ[D[F, x] - (%s), Assumptions -> Element[x, Complexes]]]", elem[i], elem[i]);
        BR_TRUE(buf);
    }
}

static void test_pmt_confirmation(void) {
    /* dependent generators: no certificate any more */
    const char* dep[] = { "(I E^(I x) + 2 E^(2 x) + I E^((2 + I) x))/(1 + E^(2 x))",
                          "Together[D[Log[E^(2 x) + x] + E^(I x), x]]",
                          "Together[D[ArcTan[Sqrt[E^(2 x) + x]] + E^(I x), x]]", NULL };
    for (int i = 0; dep[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "MatchQ[Quiet[Integrate`ParallelMixedTower[%s, x]], {\"not elementary\", ___}]", dep[i]);
        BR_FALSE(buf);
    }
    /* genuine certificates are confirmed and stand */
    BR_EQ("Quiet[Integrate`ParallelMixedTower[E^(x^2), x]][[1]]", "\"not elementary\"");
    BR_EQ("Quiet[Integrate`ParallelMixedTower[1/Log[x], x]][[1]]", "\"not elementary\"");
    BR_EQ("Quiet[Integrate`ParallelMixedTower[Sqrt[x^3 + 1], x]][[1]]", "\"not elementary\"");
    BR_EQ("ParallelMixed`Private`pmtConfirmCertificate[x, x, {\"failed\", \"r\"}]", "{\"failed\", \"r\"}");
}

static void test_eiq_transcendental_constant(void) {
    /* D[ArcTan[u] Log[2 u]] in logarithmic form: not certified any more */
    BR_FALSE("Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + (I/2) (Log[1 - I u] - Log[1 + I u])/u + Log[u]/(1 + u^2), u]] === False");
    BR_FALSE("Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + ArcTan[u]/u + Log[u]/(1 + u^2), u]] === False");
    BR_EQ("Risch`ElementaryIntegralQ[(I/2) (Log[1 - I u] - Log[1 + I u])/u + Log[u]/(1 + u^2), u]", "True");
    BR_EQ("Risch`ElementaryIntegralQ[E^(x^2), x]", "False");
    /* the ::undec message honours Quiet and Check */
    BR_EQ("Check[Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + Log[u]/(1 + u^2), u]]; \"ok\", \"msg\"]", "\"ok\"");
    BR_EQ("Check[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + Log[u]/(1 + u^2), u]; \"ok\", \"msg\"]", "\"msg\"");
}

static void test_bronstein_placeholder(void) {
    BR_TRUE("With[{f = (2 + Log[x])/(-2 x Sqrt[Log[x]] + 2 x Log[x]^(3/2))}, With[{F = BronsteinIntegrate[f, x]},"
            " !ListQ[F] && PossibleZeroQ[D[F, x] - f]]]");
    BR_EQ("Check[BronsteinIntegrate[(2 + Log[x])/(-2 x Sqrt[Log[x]] + 2 x Log[x]^(3/2)), x]; \"ok\", \"msg\"]", "\"ok\"");
}

static void test_trager_gaussian_trace(void) {
    const char* elem[] = { "x^2 (1 + I x^4)^(1/4)", "(1 + I x^4)^(-1/4)", "x^2 (2 - I x^4)^(1/4)", NULL };
    for (int i = 0; elem[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "With[{F = Trager`TragerIntegrate[%s, x]}, !ListQ[F] &&"
                 " PossibleZeroQ[D[F, x] - (%s), Assumptions -> Element[x, Complexes]]]", elem[i], elem[i]);
        BR_TRUE(buf);
    }
    /* the rational analogues are unchanged: Chebyshev decides both ways */
    BR_EQ("Trager`TragerIntegrate[(1 + x^4)^(1/4), x][[1]]", "\"not elementary\"");
    BR_TRUE("!ListQ[Trager`TragerIntegrate[x^2 (1 + x^4)^(1/4), x]]");
    BR_EQ("Block[{tragerConstantI = True}, traceCertificate[<|\"Z\" -> z, \"R\" -> z^2 - I|>]]", "Null");
}

/* G-nf-7: an algebraic constant that is zero but unsimplified must not reach
 * the integrator as a nonzero coefficient (E^(0 x^2)/(2*0) was returned and
 * PossibleZeroQ accepted it) */
static void test_zero_algebraic_constant(void) {
    BR_TRUE("Integrate[Exp[(GoldenRatio^2 - GoldenRatio - 1) x^2] x + 1, x] === x + x^2/2");
    BR_TRUE("Integrate[Exp[x^2 + (GoldenRatio^2 - GoldenRatio - 1) x] x, x] === E^x^2/2");
    BR_TRUE("Integrate[Exp[(Sqrt[2]^2 - 2 + 1) x], x] === E^x");
    BR_EQ("Integrate[((Sqrt[3] - 1) (Sqrt[3] + 1) - 2) Log[x], x]", "0");
    /* the caller's integrand is not mutated (expr_copy is a refcount bump) */
    BR_EQ("f0 = Exp[(GoldenRatio^2 - GoldenRatio - 1) x^2] x; Integrate[f0, x]; Head[f0]", "Times");
}

/* G-core-2: rt_try_erf took the constant factor as f|_{x=0}, a singular
 * point of Log[2 x] - Log[x]: the "antiderivative" was Indeterminate */
static void test_erf_constant_at_singular_point(void) {
    BR_FALSE("FreeQ[Integrate[(Log[x^2] - 2 Log[x]) Exp[x^2], x], Erfi | Erf]");
    BR_FALSE("Integrate[(Log[x^2] - 2 Log[x]) Exp[x^2], x] === Indeterminate");
    BR_TRUE("With[{F = Integrate`RischTranscendental[(Log[2 x] - Log[x]) Exp[x^2], x]}, "
            "!FreeQ[F, Erf] && PossibleZeroQ[D[F, x] - (Log[2 x] - Log[x]) Exp[x^2]]]");
    BR_TRUE("!FreeQ[Integrate`RischTranscendental[c Exp[x^2], x], Erf]");
}

/* G-core-3: a Throw[..., "PIM"] from the ParallelMixed valuation routines
 * escaped ParallelMixedSpecial to the top level */
static void test_pim_throw_caught(void) {
    BR_EQ("Check[Head[Integrate[Exp[x + Log[x^2]/4], x]], \"msg\"]", "Integrate");
}

/* G-core-5: machine ArcSech of a negative real took the lower side of the cut */
static void test_arcsech_negative_real(void) {
    BR_TRUE("Sign[Im[ArcSech[-3.7]]] == Sign[Im[N[ArcSech[-37/10], 20]]] == 1");
    BR_TRUE("Abs[ArcSech[-0.5] - ArcCosh[-2.]] < 10^-14");
    BR_TRUE("Abs[ArcSech[-0.5] - N[ArcSech[-1/2], 20]] < 10^-14");
    BR_TRUE("Max[Abs[Normal[ArcSech[NDArray[{-3.7, -0.5, 0.5, 2.}]]] - (ArcSech /@ {-3.7, -0.5, 0.5, 2.})]] == 0");
}

/* G-adv-b3: the transcendental integrator expands Log arguments (Log[x^2] ->
 * 2 Log[x]) and drops the locally constant kappa; the kappa split now makes
 * that expansion explicit (brLogExpand) */
static void test_log_expansion_kappa(void) {
    BR_EQ("brLogExpand[Log[x^2 (1 + x)]]", "2 Log[x] + Log[1 + x]");
    BR_EQ("brLogExpand[Log[-x]]", "I Pi + Log[x]");
    BR_EQ("brLogExpand[Log[1 + x]]", "Log[1 + x]");
    BR_EQ("brLogExpand[Sin[x]]", "Sin[x]");
    BR_TRUE("brLogExpandExactQ[3 x] && brLogExpandExactQ[x/2]");
    BR_FALSE("brLogExpandExactQ[x^2] || brLogExpandExactQ[-x] || brLogExpandExactQ[2 x^2]");
    BR_ANTIDERIV("Log[x^2]/(1 + x) + 2 Log[1 + x]/x", "Bronstein");
    BR_ANTIDERIV("Log[-x]/(1 + x) + Log[1 + x]/x", "Bronstein");
    BR_ANTIDERIV("Log[1/x]/(1 + x) - Log[1 + x]/x", "Bronstein");
    BR_ANTIDERIV("Log[x^3]", "Bronstein");
    /* one kappa, not a second (zero) one for Log[x] met again on the surface */
    BR_EQ("Length[brReducibleLogSplit[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]], x][\"kappas\"]]", "1");
}

/* G-sign-2: PossibleZeroQ read a residual as non-zero at a real point ON a
 * radical's cut reached through complex arithmetic (Tan[x] as a TrigToExp
 * quotient), where rounding picks the side at every precision; such a point
 * is now re-drawn */
static void test_ambiguous_cut_point(void) {
    BR_ANTIDERIV("Sqrt[Tan[x]] + Sqrt[Cot[x]]", "Bronstein");
    /* a cut reached through REAL arithmetic is still seen (x < -1) */
    BR_FALSE("PossibleZeroQ[Sqrt[(1 - x)/(1 + x)] Sqrt[(1 + x)/(1 - x)] - 1, Assumptions -> Element[x, Complexes]]");
}

/* G-adv-p4: the PMT surface entry answers in its documented protocol only */
static void test_pmt_status_protocol(void) {
    BR_EQ("Quiet[Integrate`ParallelMixedTower[Together[D[ArcTan[Log[x]] + Sqrt[x^3 - x] (2 x - 3), x]], x]][[1]]",
          "\"failed\"");
}

/* G-core-6: Together over Q(Root) (poly.c: Root content of a univariate gcd,
 * exact division by a Root-carrying divisor) */
static void test_together_root(void) {
    /* with the extension, the squarefree decomposition over Q(2^(1/3)) */
    BR_EQ("FactorSquareFreeList[z^2 - 2^(1/3) z^3, Extension -> Automatic]", "{{1, 1}, {z, 2}, {1 - 2^(1/3) z, 1}}");
    BR_TRUE("With[{r = Root[1 + #1 + #1^4 &, 1]}, With[{e = 1/(1 + 2 t + t^2) + 1/(1 + t + 4 r^3 + 4 t r^3)},"
            " Head[Together[e]] =!= Plus && PossibleZeroQ[Together[e] - e]]]");
    BR_TRUE("With[{r = Root[1 + #1 + #1^4 &, 1]}, With[{e = D[(1 + t + 2 t r + 2 t r^2 + 4 r^3 + 4 t r^3 + 2 t r^5)/(1 + t + 4 r^3 + 4 t r^3), t]},"
            " Head[Together[e]] =!= Plus && PossibleZeroQ[Together[e] - e]]]");
}

/* G-adv-b2: a literal I in the integrand, residues in Q(i), over a Log tower:
 * the log-part remainder met I as a literal on one side and as the residue
 * field's generator on the other */
static void test_logpart_remainder_qi(void) {
    BR_ANTIDERIV("(I/2)/(x Sqrt[Log[x]]) + 1/(2 x Sqrt[Log[x]] (1 + Log[x]))", "Bronstein");
    BR_ANTIDERIV("I/(x Sqrt[Log[x]]) + 1/(x Sqrt[Log[x]] (1 + Log[x]))", "Bronstein");
}

int main(void) {
    br_init();
    test_zero_algebraic_constant();
    test_erf_constant_at_singular_point();
    test_pim_throw_caught();
    test_arcsech_negative_real();
    test_log_expansion_kappa();
    test_ambiguous_cut_point();
    test_pmt_status_protocol();
    test_together_root();
    test_logpart_remainder_qi();
    test_pmt_root_residues();
    test_pmt_confirmation();
    test_eiq_transcendental_constant();
    test_bronstein_placeholder();
    test_trager_gaussian_trace();
    BR_DONE("adversarial_regressions_tests");
    return 0;
}
