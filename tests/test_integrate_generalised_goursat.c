/* test_integrate_generalised_goursat.c
 *
 * Unit tests for Integrate's Method -> "GeneralisedGoursat": the .m package
 * src/internal/goursat/ (a port of the sympy reference goursat.py of the
 * cube-root Goursat paper, 2026-11-26).
 *
 * Every reduction is checked two ways: its verdict (elementary /
 * non-elementary / undecided) against independent ground truth -- Chebyshev's
 * theorem, the paper's worked examples and its "Sample results" table, the
 * review's counterexamples -- and its answer by the package's exact verifier
 * (formal differentiation on y^n = R), which is itself checked to REJECT
 * tampered answers.  The end-to-end tests go through Integrate[..., Method ->
 * "GeneralisedGoursat"] and re-check by PossibleZeroQ.  The large constructed
 * corpora (tests/data/generalised_goursat/) run in the nightly
 * generalised_goursat_stress_tests.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/* Evaluate `input`, assert its printed form equals `expected`. */
static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

static void run(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* GoursatReduce[F, R, t, p] has status `status` and its answer verifies. */
static void reduce(const char* F, const char* R, const char* p, const char* status) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "With[{gr = GeneralisedGoursat`GoursatReduce[%s, %s, t, %s]},"
             " {gr[\"Status\"], GeneralisedGoursat`GoursatVerify[gr, %s, t]}]",
             F, R, p, F);
    char want[128];
    snprintf(want, sizeof(want), "{\"%s\", True}", status);
    check_eq(buf, want);
}

/* ... and its residual is `residual` (printed form). */
static void reduce_residual(const char* F, const char* R, const char* p,
                            const char* status, const char* residual) {
    reduce(F, R, p, status);
    char buf[2048];
    snprintf(buf, sizeof(buf), "GeneralisedGoursat`GoursatReduce[%s, %s, t, %s][\"Residual\"]",
             F, R, p);
    check_eq(buf, residual);
}

/* GoursatIntegrateAlgebraic[f, t] has status `status` and verifies. */
static void algebraic(const char* f, const char* status) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "With[{gr = GeneralisedGoursat`GoursatIntegrateAlgebraic[%s, t]},"
             " {gr[\"Status\"], GeneralisedGoursat`GoursatVerifyExpr[gr, %s, t]}]", f, f);
    char want[128];
    snprintf(want, sizeof(want), "{\"%s\", True}", status);
    check_eq(buf, want);
}

/* Integrate[f, t, Method -> "GeneralisedGoursat"] closes and differentiates back. */
static void ok(const char* f) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "FreeQ[Integrate[%s, t, Method -> \"GeneralisedGoursat\"], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[D[Integrate[%s, t, Method -> \"GeneralisedGoursat\"], t] - (%s)]", f, f);
    check_eq(buf, "True");
}

/* The method declines: Integrate stays unevaluated. */
static void declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Head[Integrate[%s, t, Method -> \"GeneralisedGoursat\"]]", f);
    check_eq(buf, "Integrate");
}

/* ------------------------------------------------------------------ */
static void test_building_blocks(void) {
    /* descent along a rational map: (t + 1/t)^2 */
    check_eq("GeneralisedGoursat`GoursatDescend[(t^2 + 1)^2/t^2, t + 1/t, t, x]", "x^2");
    check_eq("GeneralisedGoursat`GoursatDescend[(t^4 + 1)/t^2, t + 1/t, t, x]", "-2 + x^2");
    /* coordinates over Q(z^3), z = t: t^4 = t * t^3 */
    check_eq("GeneralisedGoursat`GoursatCoords[1 + t^4 + t^5, {1, t, t^2}, t^3, t, x]",
             "{1, x, x}");
    /* normal form in Q(t)[y]/(y^3 - R) */
    check_eq("Together[GeneralisedGoursat`GoursatNormalize[y^4 + 1/y, t^3 - 1, t, y, 3] - "
             "((-1 + t^3) y + y^2/(-1 + t^3))]", "0");
    /* covariants: H and G of a cubic satisfy the syzygy G^2 = a3 H^3 + a0 R^2 with
     * x0^2 = 4 w^3 - 27 disc  (H = 3 (b^2-3ac)...): check via the second-level identity */
    check_eq("With[{c = GeneralisedGoursat`GoursatCubicCovariants[t^3 - 1, t]}, c]",
             "{9 t, -27 - 27 t^3}");
    check_eq("With[{c = GeneralisedGoursat`GoursatQuarticCovariants[t^4 - 1, t]}, Expand[c[[1]]]]",
             "-t^2");
    /* the involutions of y^2 = (t^2-1)(t^2-4): t -> -t, 2/t, -2/t */
    check_eq("Sort[GeneralisedGoursat`GoursatInvolutions[(t^2 - 1)(t^2 - 4), t]]",
             "{{0, 1, 0}, {1, 0, -2}, {1, 0, 2}}");
    /* Goursat's F0 is a projection: F0[F0] == F0, F0[F - F0] == 0 (irreducible resolvent) */
    run("gg$F0 = GeneralisedGoursat`GoursatF0[1/(t - 2), t^4 + t + 1, t]");
    check_eq("Together[GeneralisedGoursat`GoursatF0[gg$F0, t^4 + t + 1, t] - gg$F0]", "0");
    check_eq("GeneralisedGoursat`GoursatF0[Together[1/(t - 2) - gg$F0], t^4 + t + 1, t]", "0");
    /* verified rational integration over an algebraic-number field */
    check_eq("With[{P = GeneralisedGoursat`GoursatRationalIntegrate[1/(m^2 - Sqrt[2]), m]},"
             " PossibleZeroQ[D[P, m] - 1/(m^2 - Sqrt[2])]]", "True");
}

static void test_hermite(void) {
    /* an exact differential reduces to zero: F = d(A y^(n-k)) */
    check_eq("Rest[GeneralisedGoursat`GoursatHermite[3 t^4 - 10 t^2 + 4, (t^2-1)(t^2-4), t, 2, 1]]",
             "{0, 0, 1}");
    check_eq("GeneralisedGoursat`GoursatHermite[3 t^4 - 10 t^2 + 4, (t^2-1)(t^2-4), t, 2, 1][[1]]",
             "t");
    check_eq("GeneralisedGoursat`GoursatHermite[t (2 t^3 - 1), t^3 - 1, t, 3, 1][[1]]", "1/2 t^2");
    /* poles at branch points are removed entirely; a pole elsewhere stays simple */
    check_eq("GeneralisedGoursat`GoursatHermite[1/(t - 1)^2, t^3 - 1, t, 3, 1][[3]]", "0");
    check_eq("GeneralisedGoursat`GoursatHermite[1/(t + 3)^2, t^3 - 1, t, 3, 1][[4]]", "3 + t");
}

static void test_classify(void) {
    /* y^2 = cubic: dt/y first kind (not exact, no residues) */
    check_eq("GeneralisedGoursat`GoursatClassify[1, 0, t^3 + 1, t, 2, 1]", "\"second-kind\"");
    /* y^3 = t^3 - 1: dt/y has poles at the three points over infinity */
    check_eq("GeneralisedGoursat`GoursatClassify[1, 0, t^3 - 1, t, 3, 1]", "\"third-kind\"");
    /* t dt/y on y^3 = t^3 - 1: second kind */
    check_eq("GeneralisedGoursat`GoursatClassify[t, 0, t^3 - 1, t, 3, 1]", "\"second-kind\"");
    /* simple poles at non-branch points */
    check_eq("GeneralisedGoursat`GoursatClassify[0, 1, t^3 - 1, t, 3, 1]", "\"third-kind\"");
    check_eq("GeneralisedGoursat`GoursatClassify[0, 0, t^3 - 1, t, 3, 1]", "\"zero\"");
}

/* The paper's "Sample results" table (section Implementation). */
static void test_sample_results(void) {
    reduce_residual("t", "(t^2-1)(t^2-4)", "1/2", "elementary", "0");
    reduce_residual("3 t^4 - 10 t^2 + 4", "(t^2-1)(t^2-4)", "1/2", "elementary", "0");
    reduce("t/(t^3 + 8)", "t^3 - 1", "1/2", "elementary");
    reduce("1/t", "t^3 - 1", "1/2", "elementary");
    reduce_residual("t", "t^3 + 1", "1/2", "non-elementary", "t");
    reduce_residual("1", "t^4 + 1", "1/2", "non-elementary", "1");
    reduce("1", "t^3 - 1", "1/3", "elementary");
    reduce_residual("t", "t^3 - 1", "1/3", "non-elementary", "t");
    reduce_residual("1", "t^2 - 1", "1/3", "non-elementary", "1");
    reduce("t (2 t^3 - 1)", "t^3 - 1", "1/3", "elementary");
    check_eq("GeneralisedGoursat`GoursatReduce[t (2 t^3 - 1), t^3 - 1, t, 1/3][\"Antiderivative\"]",
             "1/2 t^2 y^2");
    reduce("1", "t^3 + t", "1/3", "elementary");
    reduce("1 + 5 t^2 - 7 t/(t^3 + 1)", "t^3 - 1", "1/3", "elementary");
    reduce("1", "t^3 + t + 1", "1/3", "undecided");
    /* second level: the residual is a square-root integral on 4 w^3 - 27 disc */
    reduce("(t^6 - 2 t^3 + 4)/((t^3 - 1) (t^3 - 4))", "t^3 - 2", "2/3", "elementary");
    reduce("t", "t^3 - 1", "2/3", "elementary");
    reduce_residual("1", "t^3 - 1", "2/3", "non-elementary", "1");
}

/* Chebyshev: Int t^m (t^n - 1)^(-p) dt is elementary iff (m+1)/n or
 * (m+1)/n - p is an integer -- independent ground truth for every verdict. */
static void test_chebyshev(void) {
    static const int ns[] = { 3, 3, 2, 2, 3, 4 };
    static const char* ps[] = { "1/3", "2/3", "1/3", "2/3", "1/2", "1/2" };
    static const int pnum[] = { 1, 2, 1, 2, 1, 1 };
    static const int pden[] = { 3, 3, 3, 3, 2, 2 };
    for (int c = 0; c < 6; c++) {
        for (int m = -5; m <= 7; m++) {
            /* (m+1)/n integer, or (m+1)/n - p integer:  (m+1) pden - n pnum  = 0 mod n pden */
            int n = ns[c];
            bool a = ((m + 1) % n) == 0;
            int num = (m + 1) * pden[c] - n * pnum[c];
            bool b = (num % (n * pden[c])) == 0;
            char F[32], R[32];
            snprintf(F, sizeof(F), "t^(%d)", m);
            snprintf(R, sizeof(R), "t^%d - 1", n);
            reduce(F, R, ps[c], (a || b) ? "elementary" : "non-elementary");
        }
    }
}

/* test_goursat.py: the papers' examples and the review's counterexamples. */
static void test_reference_cases(void) {
    reduce("(t^3 - 4)/(t^3 + 8)", "t^3 - 1", "1/2", "elementary");
    reduce("t^2/(t^3 + 8)", "t^3 - 1", "1/2", "elementary");
    reduce("(t + 1)^2/((t + 1)^3 + 2)", "(t + 1)^3 - 2", "1/2", "elementary");
    reduce("t^2/(t^3 + 5)", "t^3 - 2", "1/2", "elementary");
    reduce("t/(t^3 + 2)", "t^4 - t", "1/2", "elementary");
    reduce("1", "t^4 - t", "1/2", "non-elementary");
    reduce("(t - 1)/(t + 2)", "t^3 - 1", "1/2", "elementary");
    reduce("(t - 1)/(t + 2) + t^2/(t^3 + 8) + 1/t", "t^3 - 1", "1/2", "elementary");
    reduce("t^2", "t^4 - 1", "1/2", "non-elementary");
    reduce("1/(t^4 + 3)", "t^4 - 1", "1/2", "undecided");
    reduce("(t^2 - 1)/(t^2 + 1)", "t^4 + 1", "1/2", "elementary");
    reduce("2 t^3 - 1", "t^3 - 1", "2/3", "elementary");
    reduce("(2 t + 5)/(3 t^2 - 12 t + 13)", "(t - 1)(t - 2)(t - 3)", "1/3", "elementary");
    reduce("1", "t^3 + t + 1", "1/2", "non-elementary");
    reduce("(t^3 + 1)/(t^3 - 1)", "t^3 - 1", "2/3", "elementary");
    reduce("1", "(t - 1)(t - 2)(t - 3)", "1/3", "elementary");
    reduce("t", "(t - 1)(t - 2)(t - 3)", "2/3", "non-elementary");
    reduce("1", "t^2 + t + 1", "2/3", "non-elementary");
}

/* The irreducible-resolvent case, which goursat.py leaves unevaluated:
 * carried out over the cubic field (CubicField.m). */
static void test_cubic_field(void) {
    run("gg$Fc = Together[1/(t - 2) - GeneralisedGoursat`GoursatF0[1/(t - 2), t^4 + t + 1, t]]");
    run("gg$rc = GeneralisedGoursat`GoursatReduce[gg$Fc, t^4 + t + 1, t, 1/2]");
    /* Goursat's criterion is decided over Q (the trace): elementary */
    check_eq("gg$rc[\"Status\"]", "\"elementary\"");
    /* the explicit reduction over the cubic field is returned only verified;
     * otherwise (budget, or no verifiable form) it stays goursat.py's
     * unevaluated integral -- never an unverified closed form */
    check_eq("If[FreeQ[gg$rc[\"Antiderivative\"], GeneralisedGoursat`Private`ggIntegral],"
             " GeneralisedGoursat`GoursatVerify[gg$rc, gg$Fc, t], True]", "True");
    /* the V4-invariant part of 1/(t-2) is of the third kind: undecided */
    check_eq("GeneralisedGoursat`GoursatReduce[1/(t - 2), t^4 + t + 1, t, 1/2][\"Status\"]",
             "\"undecided\"");
}

static void test_verifier_rejects(void) {
    /* the verifier must say False for a wrong antiderivative */
    run("gg$r = GeneralisedGoursat`GoursatReduce[t, (t^2-1)(t^2-4), t, 1/2]");
    check_eq("GeneralisedGoursat`GoursatVerify[gg$r, t, t]", "True");
    check_eq("GeneralisedGoursat`GoursatVerify[Append[gg$r, \"Antiderivative\" -> 2 gg$r[\"Antiderivative\"]], t, t]",
             "False");
    check_eq("GeneralisedGoursat`GoursatVerify[Append[gg$r, \"Antiderivative\" -> gg$r[\"Antiderivative\"] + Log[t]], t, t]",
             "False");
    run("gg$r = GeneralisedGoursat`GoursatReduce[t/(t^3 + 8), t^3 - 1, t, 1/2]");
    check_eq("GeneralisedGoursat`GoursatVerify[gg$r, t/(t^3 + 8), t]", "True");
    check_eq("GeneralisedGoursat`GoursatVerify[Append[gg$r, \"Antiderivative\" -> 3 gg$r[\"Antiderivative\"]], t/(t^3 + 8), t]",
             "False");
    check_eq("GeneralisedGoursat`GoursatVerify[Append[gg$r, \"Residual\" -> 1], t/(t^3 + 8), t]", "False");
}

/* The front end: any element of Q(t, R^(1/n)), split into eigencomponents. */
static void test_front_end(void) {
    algebraic("(t^2 - 1)/((t^2 + 1) Sqrt[t^4 + 1])", "elementary");
    algebraic("Sqrt[t^4 + 1]/t", "elementary");
    algebraic("1/Sqrt[1 - t^4]", "non-elementary");
    algebraic("1/(t^3 - 1)^(1/3) + t^2", "elementary");
    algebraic("(t^3 + 1)^(2/3)", "elementary");
    algebraic("(t^3 + 1)^(1/3)", "non-elementary");
    algebraic("1/(t (t^2 - 1)^(1/3))", "elementary");
    /* fourth roots (welz_cases.py) and y^2 = Sqrt[R] inside a fourth-root integrand */
    algebraic("((-1 + t) (-3 + 2 t^4)^(1/4))/t^2", "elementary");
    algebraic("(t^4 (-1 + 2 t^4)^(3/4))/(2 + t^4)", "elementary");
    algebraic("1/(t^3 (1 + 2 t^3)^(1/3) (2 + t^6))", "elementary");
    /* Chebyshev's torsion example and the second-level literature case */
    algebraic("t/Sqrt[t^4 + 10 t^2 - 96 t - 71]", "elementary");
    algebraic("(t^6 - 2 t^3 + 4)/((t^3 - 4) (t^3 - 2)^(2/3) (t^3 - 1))", "elementary");
    /* pseudo-elliptic integrands with algebraic coefficients (goursat_examples*.wl) */
    algebraic("(t/(t^2 + 3 - 2 Sqrt[2]))/Sqrt[t^3 + t]", "elementary");
    algebraic("((t - I)/(t + I))/Sqrt[t^4 + 1]", "elementary");
    algebraic("((t + 3 I)^(-1))/(t^2 - 3)^(1/3)", "elementary");
}

static void test_integrate_method(void) {
    ok("t/Sqrt[(t^2-1)(t^2-4)]");
    ok("1/(t^3-1)^(1/3)");
    ok("t/((t^3 + 8) Sqrt[t^3 - 1])");
    ok("(t^2 - 1)/((t^2 + 1) Sqrt[t^4 + 1])");
    ok("t (2 t^3 - 1)/(t^3 - 1)^(1/3)");
    ok("1/(t^4 - 1)^(1/4)");
}

static void test_integrate_partial(void) {
    /* certified non-elementary with an elementary part: the method symbol returns
     * {"not elementary", partial, residual} ... */
    check_eq("Integrate`GeneralisedGoursat[(3 t^4 - 10 t^2 + 5)/Sqrt[(t^2-1)(t^2-4)], t]",
             "{\"not elementary\", t Sqrt[4 - 5 t^2 + t^4], 1/Sqrt[4 - 5 t^2 + t^4]}");
    /* ... and Integrate returns partial + Integrate[residual], the residual left to
     * the rest of the cascade (here an elliptic integral) */
    check_eq("With[{r = Integrate[(t^3 + 3)/Sqrt[t^4 + 1], t, Method -> \"GeneralisedGoursat\"]},"
             " {FreeQ[r, Integrate], PossibleZeroQ[D[r, t] - (t^3 + 3)/Sqrt[t^4 + 1]]}]",
             "{True, True}");
    /* nothing elementary to return: unevaluated */
    declines("1/Sqrt[t^3 + t + 1]");
    declines("1/(t^3 + t + 1)^(1/3)");       /* undecided */
    declines("Sin[t]");                        /* out of shape */
    declines("1/(t^5 + 1)^(1/3)");             /* genus > 1 */
    /* the method symbol answers in the package's own protocol */
    check_eq("Integrate`GeneralisedGoursat[t/(t^3 - 1)^(1/3), t]",
             "{\"not elementary\", 0, t/(-1 + t^3)^(1/3)}");
    check_eq("Integrate`GeneralisedGoursat[Sin[t], t]", "$Failed");
}

int main(void) {
    symtab_init();
    core_init();
    /* load the package through its method symbol (the lazy loader under test) */
    run("Integrate`GeneralisedGoursat[1/Sqrt[t^3 - 1], t]");

    TEST(test_building_blocks);
    TEST(test_hermite);
    TEST(test_classify);
    TEST(test_sample_results);
    TEST(test_chebyshev);
    TEST(test_reference_cases);
    TEST(test_cubic_field);
    TEST(test_verifier_rejects);
    TEST(test_front_end);
    TEST(test_integrate_method);
    TEST(test_integrate_partial);

    printf("All Integrate GeneralisedGoursat tests passed!\n");
    return 0;
}
