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

int main(void) {
    br_init();
    test_pmt_root_residues();
    test_pmt_confirmation();
    test_eiq_transcendental_constant();
    test_bronstein_placeholder();
    test_trager_gaussian_trace();
    BR_DONE("adversarial_regressions_tests");
    return 0;
}
