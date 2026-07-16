/* test_risch_canonical.c — Bronstein §3.4–3.5 splitting factorization.
 *
 * Verifies the Risch` differential-field builtins against:
 *   - the worked examples 3.5.1 (SplitFactor) and 3.5.2 (SplitSquarefreeFactor)
 *     from Bronstein, Symbolic Integration I, 2nd ed.;
 *   - the normal/special classification for the three Liouvillian monomial
 *     kinds (exponential t' = t, hypertangent t' = 1 + t^2, primitive/log
 *     t' = 1/x).
 *
 * Polynomial equality is checked with Expand[Together[a - b]] == 0 (Together
 * alone combines denominators but does not expand the numerator), so the tests
 * are robust to cosmetic ordering/expansion differences.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert that the rational expression `diff` is identically zero. */
static void run_zero(const char* diff) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* The monomial of Bronstein Examples 3.5.1/3.5.2: Dt = -t^2 - 3/(2x) t + 1/(2x). */
#define DERIV "{x -> 1, t -> -t^2 - (3/(2 x)) t + 1/(2 x)}"

/* ---- Monomial derivation ---------------------------------------------- */
static void test_derivation(void) {
    /* Exponential monomial t' = t:  D[t^2] = 2 t t = 2 t^2. */
    run_zero("Risch`Derivation[t^2, {x -> 1, t -> t}] - 2 t^2");
    /* Hypertangent t' = 1 + t^2:  D[t] = 1 + t^2. */
    run_zero("Risch`Derivation[t, {x -> 1, t -> 1 + t^2}] - (1 + t^2)");
    /* Log/primitive t' = 1/x:  D[x t] = t + x (1/x) = t + 1. */
    run_zero("Risch`Derivation[x t, {x -> 1, t -> 1/x}] - (t + 1)");
}

/* ---- Example 3.5.1: SplitFactor --------------------------------------- */
#define PBIG "(4 x^4 t^5 - 4 x^3 (x + 1) t^4 + x^2 (2 x - 3) t^3 " \
             "+ x (2 x^2 + 7 x + 2) t^2 - (4 x^2 + 4 x - 1) t + 2 x - 1)"
#define PS_351 "(t^2 + t/x - (2 x - 1)/(4 x^2))"
#define PN_351 "(4 x^4 t^3 - 4 x^3 (x + 2) t^2 + 4 x^2 (2 x + 1) t - 4 x^2)"
#define SF351 "Risch`SplitFactor[" PBIG ", t, " DERIV "]"

static void test_split_factor_bronstein_351(void) {
    run_zero("Last[" SF351 "] - " PS_351);            /* special part */
    run_zero("First[" SF351 "] - " PN_351);           /* normal part */
    run_zero("First[" SF351 "] Last[" SF351 "] - " PBIG);  /* product = p */
    run_test("Risch`SpecialQ[" PS_351 ", t, " DERIV "]", "True");
    run_test("Risch`SpecialQ[" PN_351 ", t, " DERIV "]", "False");
}

/* ---- Example 3.5.2: SplitSquarefreeFactor ----------------------------- */
#define P352 "((4 x^2 t^3 - 4 x (x - 1) t^2 - (6 x - 1) t + 2 x - 1)(x t - 1)^2)"
#define SSF "Risch`SplitSquarefreeFactor[" P352 ", t, " DERIV "]"

static void test_split_squarefree_bronstein_352(void) {
    /* Two squarefree factors of multiplicity 1 and 2. */
    run_test("Length[" SSF "[[1]]]", "2");
    /* Reconstruction: p = content * (N1 N2^2)(S1 S2^2), content = lc_t(p) = 4x^4
     * (the factors are monic over k, so their product is the monic p). */
    run_zero("4 x^4 " SSF "[[1, 1]] " SSF "[[1, 2]]^2 " SSF "[[2, 1]] " SSF "[[2, 2]]^2 - " P352);
    /* The special part S1 is canonical (monic) and matches the book. */
    run_zero(SSF "[[2, 1]] - " PS_351);
    run_test(SSF "[[2, 2]]", "1");
    /* Each normal factor is normal; each special factor is special. */
    run_test("Risch`NormalQ[" SSF "[[1, 1]], t, " DERIV "]", "True");
    run_test("Risch`NormalQ[" SSF "[[1, 2]], t, " DERIV "]", "True");
    run_test("Risch`SpecialQ[" SSF "[[2, 1]], t, " DERIV "]", "True");
    /* N1, N2 are the monic (over k) representatives: book's 4x^2(t-1), x(t-1/x). */
    run_zero(SSF "[[1, 1]] - (t - 1)");
    run_zero(SSF "[[1, 2]] - (t - 1/x)");
}

/* ---- Exponential monomial: special polynomials are powers of t --------- */
static void test_exponential_monomial(void) {
    run_test("Risch`SpecialQ[t, t, {x -> 1, t -> t}]", "True");
    run_test("Risch`NormalQ[t, t, {x -> 1, t -> t}]", "False");
    /* t^2 + 1 is normal in an exponential extension. */
    run_test("Risch`SpecialQ[t^2 + 1, t, {x -> 1, t -> t}]", "False");
    run_test("Risch`NormalQ[t^2 + 1, t, {x -> 1, t -> t}]", "True");
    /* SplitFactor[t^3] = {1, t^3}. */
    run_test("First[Risch`SplitFactor[t^3, t, {x -> 1, t -> t}]]", "1");
    run_zero("Last[Risch`SplitFactor[t^3, t, {x -> 1, t -> t}]] - t^3");
    /* SplitFactor[t^2 + 1] = {t^2 + 1, 1}. */
    run_test("Last[Risch`SplitFactor[t^2 + 1, t, {x -> 1, t -> t}]]", "1");
}

/* ---- Hypertangent monomial: the only special irreducible is t^2 + 1 ---- */
static void test_hypertangent_monomial(void) {
    run_test("Risch`SpecialQ[t^2 + 1, t, {x -> 1, t -> 1 + t^2}]", "True");
    run_test("Risch`NormalQ[t, t, {x -> 1, t -> 1 + t^2}]", "True");
    run_test("Risch`SpecialQ[t, t, {x -> 1, t -> 1 + t^2}]", "False");
    /* SplitFactor[(t^2+1)^2] = {1, (t^2+1)^2}. */
    run_test("First[Risch`SplitFactor[(t^2 + 1)^2, t, {x -> 1, t -> 1 + t^2}]]", "1");
    run_zero("Last[Risch`SplitFactor[(t^2 + 1)^2, t, {x -> 1, t -> 1 + t^2}]] - (t^2 + 1)^2");
    /* A mixed denominator: t (t^2+1) splits into normal t, special t^2+1. */
    run_zero("First[Risch`SplitFactor[t (t^2 + 1), t, {x -> 1, t -> 1 + t^2}]] - t");
    run_zero("Last[Risch`SplitFactor[t (t^2 + 1), t, {x -> 1, t -> 1 + t^2}]] - (t^2 + 1)");
}

/* ---- Primitive/log monomial: no non-constant special polynomials ------- */
static void test_log_monomial(void) {
    run_test("Risch`NormalQ[t, t, {x -> 1, t -> 1/x}]", "True");
    run_test("Risch`SpecialQ[t, t, {x -> 1, t -> 1/x}]", "False");
    /* Everything is normal: SplitFactor[t^2 + 1] = {t^2 + 1, 1}. */
    run_test("Last[Risch`SplitFactor[t^2 + 1, t, {x -> 1, t -> 1/x}]]", "1");
    run_zero("First[Risch`SplitFactor[t^2 + 1, t, {x -> 1, t -> 1/x}]] - (t^2 + 1)");
    run_test("Last[Risch`SplitFactor[t^2, t, {x -> 1, t -> 1/x}]]", "1");
}

/* ---- Field polynomial division in k[t] over C(x) ----------------------- */
static void test_polydivide(void) {
    run_test("Risch`PolyDivide[t^3 + t + 1, t^2 + 1, t]", "List[t, 1]");
    /* Over the field C(x): t^2 = (t - 1/x)(t + 1/x) + 1/x^2. */
    run_zero("Risch`PolyDivide[t^2, t - 1/x, t][[1]] - (t + 1/x)");
    run_zero("Risch`PolyDivide[t^2, t - 1/x, t][[2]] - 1/x^2");
    /* Non-monic divisor with an x-coefficient leading term: reconstruction
     * (the leading coefficient x is inverted in the field). */
    run_zero("Risch`PolyDivide[x t^3 + t, x t + 1, t][[1]] (x t + 1) "
             "+ Risch`PolyDivide[x t^3 + t, x t + 1, t][[2]] - (x t^3 + t)");
}

/* ---- CanonicalRepresentation f = f_p + f_s + f_n ------------------------ */
#define EXPD "{x -> 1, t -> t}"
#define TAND "{x -> 1, t -> 1 + t^2}"
#define FEXP "(t^2 + (2 t + 3)/t^2 + (t + 1)/(t^2 + 1))"
#define CREXP "Risch`CanonicalRepresentation[" FEXP ", t, " EXPD "]"
#define FTAN "((3 t)/((t^2 + 1)(t - 1)))"
#define CRTAN "Risch`CanonicalRepresentation[" FTAN ", t, " TAND "]"

static void test_canonical_representation(void) {
    /* Exponential monomial: special denominators are powers of t. */
    run_zero(CREXP "[[1]] + " CREXP "[[2]] + " CREXP "[[3]] - " FEXP);  /* reconstruction */
    /* The decomposition is unique: poly t^2, special (2t+3)/t^2, normal (t+1)/(t^2+1). */
    run_zero(CREXP "[[1]] - t^2");
    run_zero(CREXP "[[2]] - (2 t + 3)/t^2");
    run_zero(CREXP "[[3]] - (t + 1)/(t^2 + 1)");
    run_test("PolynomialQ[" CREXP "[[1]], t]", "True");
    run_test("Risch`SpecialQ[Denominator[Together[" CREXP "[[2]]]], t, " EXPD "]", "True");
    run_test("Risch`NormalQ[Denominator[Together[" CREXP "[[3]]]], t, " EXPD "]", "True");

    /* Pure polynomial: f_s = f_n = 0, f_p = f. */
    run_zero("Risch`CanonicalRepresentation[t^2 + t + 1, t, " EXPD "][[1]] - (t^2 + t + 1)");
    run_zero("Risch`CanonicalRepresentation[t^2 + t + 1, t, " EXPD "][[2]]");
    run_zero("Risch`CanonicalRepresentation[t^2 + t + 1, t, " EXPD "][[3]]");

    /* Hypertangent monomial: special t^2+1, normal t-1; proper fraction (f_p=0). */
    run_zero(CRTAN "[[1]] + " CRTAN "[[2]] + " CRTAN "[[3]] - " FTAN);  /* reconstruction */
    run_zero(CRTAN "[[1]]");                                           /* f_p = 0 */
    run_test("Risch`SpecialQ[Denominator[Together[" CRTAN "[[2]]]], t, " TAND "]", "True");
    run_test("Risch`NormalQ[Denominator[Together[" CRTAN "[[3]]]], t, " TAND "]", "True");
}

/* ---- SplitFactor: edge cases + the defining property ------------------ */
#define LOGD "{x -> 1, t -> 1/x}"
static void test_split_factor_edge(void) {
    /* Pure special (exp): SplitFactor[t^3] = {1, t^3}. */
    run_test("First[Risch`SplitFactor[t^3, t, " EXPD "]]", "1");
    run_zero("Last[Risch`SplitFactor[t^3, t, " EXPD "]] - t^3");
    /* Pure normal (exp): SplitFactor[t^2 + 1] = {t^2 + 1, 1}. */
    run_test("Last[Risch`SplitFactor[t^2 + 1, t, " EXPD "]]", "1");
    /* Constant: SplitFactor[5] = {5, 1}. */
    run_test("First[Risch`SplitFactor[5, t, " EXPD "]]", "5");
    run_test("Last[Risch`SplitFactor[5, t, " EXPD "]]", "1");
    /* Mixed exp: t^2 (t - 1) -> normal t - 1, special t^2. */
    run_zero("First[Risch`SplitFactor[t^2 (t - 1), t, " EXPD "]] - (t - 1)");
    run_zero("Last[Risch`SplitFactor[t^2 (t - 1), t, " EXPD "]] - t^2");

    /* DEFINING PROPERTY over a deeper exp case p = t^3 (t^2 + 1)^2:
     * p_n p_s = p, p_s special, p_n has no non-constant special factor. */
#define PDEEP "(t^3 (t^2 + 1)^2)"
    run_zero("First[Risch`SplitFactor[" PDEEP ", t, " EXPD "]] "
             "Last[Risch`SplitFactor[" PDEEP ", t, " EXPD "]] - " PDEEP);
    run_test("Risch`SpecialQ[Last[Risch`SplitFactor[" PDEEP ", t, " EXPD "]], t, " EXPD "]", "True");
    run_test("Risch`SpecialQ[First[Risch`SplitFactor[" PDEEP ", t, " EXPD "]], t, " EXPD "]", "False");
    run_zero("Last[Risch`SplitFactor[" PDEEP ", t, " EXPD "]] - t^3");
#undef PDEEP

    /* Hypertangent with multiplicity: t (t^2 + 1)^2 -> normal t, special (t^2+1)^2. */
    run_zero("First[Risch`SplitFactor[t (t^2 + 1)^2, t, " TAND "]] - t");
    run_zero("Last[Risch`SplitFactor[t (t^2 + 1)^2, t, " TAND "]] - (t^2 + 1)^2");
    run_test("Risch`SpecialQ[Last[Risch`SplitFactor[t (t^2 + 1)^2, t, " TAND "]], t, " TAND "]", "True");
}

/* ---- SplitSquarefreeFactor: three distinct multiplicities ------------- */
static void test_split_squarefree_edge(void) {
    /* Exp, all-normal: p = (t-1)(t-2)^2 (t-3)^3, multiplicities 1,2,3. */
#define P3 "((t - 1)(t - 2)^2 (t - 3)^3)"
#define SSF3 "Risch`SplitSquarefreeFactor[" P3 ", t, " EXPD "]"
    run_test("Length[" SSF3 "[[1]]]", "3");
    /* Reconstruction N1 N2^2 N3^3 = p (monic; content = 1). */
    run_zero(SSF3 "[[1, 1]] " SSF3 "[[1, 2]]^2 " SSF3 "[[1, 3]]^3 - " P3);
    /* Monic squarefree normal factors. */
    run_zero(SSF3 "[[1, 1]] - (t - 1)");
    run_zero(SSF3 "[[1, 2]] - (t - 2)");
    run_zero(SSF3 "[[1, 3]] - (t - 3)");
    /* No special factors (all-normal): each S_i = 1. */
    run_test(SSF3 "[[2, 1]]", "1");
    run_test(SSF3 "[[2, 2]]", "1");
    run_test(SSF3 "[[2, 3]]", "1");
    /* Each normal factor is normal. */
    run_test("Risch`NormalQ[" SSF3 "[[1, 1]], t, " EXPD "]", "True");
    run_test("Risch`NormalQ[" SSF3 "[[1, 3]], t, " EXPD "]", "True");
#undef SSF3
#undef P3
}

/* ---- CanonicalRepresentation across all monomial kinds ---------------- */
static void test_canonical_all_kinds(void) {
    /* Hypertangent, all three parts non-zero:
     * f = t + (2t+1)/(t^2+1) + 1/(t-1)  ->  {t, (2t+1)/(t^2+1), 1/(t-1)}. */
#define FT "(t + (2 t + 1)/(t^2 + 1) + 1/(t - 1))"
#define CT "Risch`CanonicalRepresentation[" FT ", t, " TAND "]"
    run_zero(CT "[[1]] + " CT "[[2]] + " CT "[[3]] - " FT);   /* reconstruction */
    run_zero(CT "[[1]] - t");                                  /* f_p */
    run_zero(CT "[[2]] - (2 t + 1)/(t^2 + 1)");                /* f_s (special denom) */
    run_zero(CT "[[3]] - 1/(t - 1)");                          /* f_n (normal denom) */
    run_test("Risch`SpecialQ[Denominator[Together[" CT "[[2]]]], t, " TAND "]", "True");
    run_test("Risch`NormalQ[Denominator[Together[" CT "[[3]]]], t, " TAND "]", "True");
#undef CT
#undef FT

    /* Log/primitive: NOTHING is special, so f_s = 0 always.
     * f = t^2 + (t+1)/(t^2+1) -> {t^2, 0, (t+1)/(t^2+1)}. */
#define FL "(t^2 + (t + 1)/(t^2 + 1))"
#define CL "Risch`CanonicalRepresentation[" FL ", t, " LOGD "]"
    run_zero(CL "[[1]] + " CL "[[2]] + " CL "[[3]] - " FL);
    run_zero(CL "[[1]] - t^2");
    run_zero(CL "[[2]]");                                      /* f_s = 0 */
    run_zero(CL "[[3]] - (t + 1)/(t^2 + 1)");
#undef CL
#undef FL

    /* Special-only (exp): f = (2t+3)/t^2 -> f_p = 0, f_n = 0. */
    run_zero("Risch`CanonicalRepresentation[(2 t + 3)/t^2, t, " EXPD "][[1]]");
    run_zero("Risch`CanonicalRepresentation[(2 t + 3)/t^2, t, " EXPD "][[3]]");
    run_zero("Risch`CanonicalRepresentation[(2 t + 3)/t^2, t, " EXPD "][[2]] - (2 t + 3)/t^2");

    /* Normal-only (exp): f = 1/(t^2+1) -> f_p = 0, f_s = 0. */
    run_zero("Risch`CanonicalRepresentation[1/(t^2 + 1), t, " EXPD "][[1]]");
    run_zero("Risch`CanonicalRepresentation[1/(t^2 + 1), t, " EXPD "][[2]]");
    run_zero("Risch`CanonicalRepresentation[1/(t^2 + 1), t, " EXPD "][[3]] - 1/(t^2 + 1)");
}

/* ---- Robustness: malformed input stays unevaluated -------------------- */
static void test_canonical_robustness(void) {
    run_test("Head[Risch`SplitFactor[t, t, {x, t}]] === Risch`SplitFactor", "True");
    run_test("Head[Risch`CanonicalRepresentation[t, t, {x, t}]] === Risch`CanonicalRepresentation", "True");
    run_test("Head[Risch`SplitFactor[t, 5, " EXPD "]] === Risch`SplitFactor", "True");  /* t not symbol */
    run_test("Head[Risch`NormalQ[t, t]] === Risch`NormalQ", "True");                     /* wrong arity */
}

int main(void) {
    core_init();

    TEST(test_derivation);
    TEST(test_split_factor_bronstein_351);
    TEST(test_split_squarefree_bronstein_352);
    TEST(test_exponential_monomial);
    TEST(test_hypertangent_monomial);
    TEST(test_log_monomial);
    TEST(test_polydivide);
    TEST(test_canonical_representation);
    TEST(test_split_factor_edge);
    TEST(test_split_squarefree_edge);
    TEST(test_canonical_all_kinds);
    TEST(test_canonical_robustness);

    printf("All risch_canonical tests passed.\n");
    return 0;
}
