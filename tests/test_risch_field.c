/* test_risch_field.c — the differential-field primitives (Bronstein §3.4).
 *
 * These are the arithmetic foundation the whole transcendental Risch stack
 * stands on: the monomial derivation D, field gcd / exact division / extended
 * Euclidean / Diophantine solve in k[t] over k = C(x, lower monomials), the
 * numerator/denominator split, polynomial division, the normal/special
 * classification (Def. 3.4.2), and PolynomialReduce (§5.4).  Each is exercised
 * in isolation with:
 *   - concrete values (including the field-unit behaviour of pure-k factors),
 *   - property round-trips (a = q b + r; u a + v b = g; b dn + c ds = r; a/d = f),
 *   - degree bounds (deg r < deg b, deg b < deg ds, deg r < delta), and
 *   - robustness (malformed derivation, wrong arity, division by zero) — every
 *     such call must stay UNEVALUATED, never crash.
 *
 * Polynomial equality is checked with Expand[Together[a - b]] == 0.
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

/* Assert the rational expression `diff` is identically zero. */
static void run_zero(const char* diff) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* Assert `expr` stays unevaluated with the given Risch` head symbol. */
static void run_held(const char* expr, const char* head) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Head[%s] === %s", expr, head);
    run_test(buf, "True");
}

/* The three Liouvillian monomial derivations plus a genuine nonlinear one. */
#define EXPD "{x -> 1, t -> t}"          /* exponential   t' = t         */
#define LOGD "{x -> 1, t -> 1/x}"        /* primitive/log t' = 1/x       */
#define TAND "{x -> 1, t -> 1 + t^2}"    /* hypertangent  t' = 1 + t^2   */
#define NLD  "{x -> 1, t -> t^2}"        /* nonlinear     t' = t^2       */

/* =====================================================================
 * Monomial derivation  D[p] = Sum_i Dvar_i (dp/dvar_i)
 * ===================================================================*/
static void test_derivation(void) {
    /* Power/product rule, exponential t' = t:  D[t^2] = 2 t t = 2 t^2. */
    run_zero("Risch`Derivation[t^2, " EXPD "] - 2 t^2");
    run_zero("Risch`Derivation[t^3, " EXPD "] - 3 t^3");
    /* Quotient rule: D[1/t] = -t^-2 t = -1/t. */
    run_zero("Risch`Derivation[1/t, " EXPD "] + 1/t");
    /* A symbol absent from the derivation is a constant: D[a t] = a t. */
    run_zero("Risch`Derivation[a t, " EXPD "] - a t");

    /* Hypertangent t' = 1 + t^2:  D[t] = 1 + t^2,  D[t^2] = 2 t (1+t^2). */
    run_zero("Risch`Derivation[t, " TAND "] - (1 + t^2)");
    run_zero("Risch`Derivation[t^2, " TAND "] - 2 t (1 + t^2)");

    /* Log/primitive t' = 1/x:  D[x t] = t + x/x = t + 1;  D[t] = 1/x. */
    run_zero("Risch`Derivation[x t, " LOGD "] - (t + 1)");
    run_zero("Risch`Derivation[t, " LOGD "] - 1/x");

    /* Base-field element (Leibniz w.r.t. x):  D[x^2] = 2 x. */
    run_zero("Risch`Derivation[x^2, " EXPD "] - 2 x");

    /* Two-monomial tower  {t1 = log x, Dt1 = 1/x ; t2 = exp x, Dt2 = t2}:
     * D[t1 t2] = (1/x) t2 + t1 t2. */
#define TOWER "{x -> 1, t1 -> 1/x, t2 -> t2}"
    run_zero("Risch`Derivation[t1 t2, " TOWER "] - (t2/x + t1 t2)");
    run_zero("Risch`Derivation[t1^2, " TOWER "] - 2 t1/x");
    run_zero("Risch`Derivation[t2^2, " TOWER "] - 2 t2^2");
#undef TOWER

    /* Constant maps to zero. */
    run_test("Risch`Derivation[5, " EXPD "]", "0");
}

/* =====================================================================
 * Field gcd in k[t]  (monic in t; pure-k factors are units)
 * ===================================================================*/
static void test_field_gcd(void) {
    /* Coprime -> 1. */
    run_test("Risch`FieldGCD[t, t + 1, t]", "1");
    run_test("Risch`FieldGCD[t^2 + 1, t - 1, t]", "1");
    /* Shared factor, made monic. */
    run_zero("Risch`FieldGCD[t^2 - 1, t - 1, t] - (t - 1)");
    run_zero("Risch`FieldGCD[(t - 1)(t - 2), (t - 1)(t - 3), t] - (t - 1)");
    /* gcd(a, 0) = a (made monic); gcd(0, b) = b. */
    run_zero("Risch`FieldGCD[t - 1, 0, t] - (t - 1)");
    run_zero("Risch`FieldGCD[0, t^2 + 1, t] - (t^2 + 1)");
    /* Pure-k factor is a UNIT: gcd((x+1) t, x, t) = 1 (x shares nothing in t). */
    run_test("Risch`FieldGCD[(x + 1) t, x, t]", "1");
    /* The x-content is divided out: gcd(x(t^2-1), x(t-1)) = t - 1, not x(t-1). */
    run_zero("Risch`FieldGCD[x (t^2 - 1), x (t - 1), t] - (t - 1)");
    /* Leading-coefficient normalisation: gcd(2 t^2 - 2, 3 t - 3) = t - 1 monic. */
    run_zero("Risch`FieldGCD[2 t^2 - 2, 3 t - 3, t] - (t - 1)");
    /* Result is monic in t. */
    run_test("Coefficient[Risch`FieldGCD[x (t^2 - 1), x (t - 1), t], t, 1]", "1");
}

/* =====================================================================
 * Field exact division  a/b in k[t] (NULL/unevaluated when inexact)
 * ===================================================================*/
static void test_div_exact(void) {
    run_zero("Risch`DivExact[t^2 - 1, t - 1, t] - (t + 1)");
    run_zero("Risch`DivExact[t^3 - 1, t - 1, t] - (t^2 + t + 1)");
    /* Pure-k divisor is a unit: (x t^2 - x)/x = t^2 - 1. */
    run_zero("Risch`DivExact[x t^2 - x, x, t] - (t^2 - 1)");
    /* Non-monic (leading x-coeff) divisor, exact: (x t - 1) t = x t^2 - t. */
    run_zero("Risch`DivExact[x t^2 - t, x t - 1, t] - t");
    /* Inexact -> stays unevaluated. */
    run_held("Risch`DivExact[t, t + 1, t]", "Risch`DivExact");
    run_held("Risch`DivExact[t^2, t^2 + 1, t]", "Risch`DivExact");
}

/* =====================================================================
 * Numerator/denominator over k(t)  (denominator monic in t)
 * ===================================================================*/
static void test_num_den(void) {
    /* Proper fraction. */
    run_zero("Risch`NumDen[(t + 1)/(t^2 + 1), t][[1]] - (t + 1)");
    run_zero("Risch`NumDen[(t + 1)/(t^2 + 1), t][[2]] - (t^2 + 1)");
    /* Reconstruction a/d = f. */
    run_zero("Risch`NumDen[(t + 1)/(t^2 + 1), t][[1]] / Risch`NumDen[(t + 1)/(t^2 + 1), t][[2]] - (t + 1)/(t^2 + 1)");
    /* Improper fraction: numerator keeps its high degree. */
    run_zero("Risch`NumDen[t^3/(t^2 + 1), t][[1]] - t^3");
    run_zero("Risch`NumDen[t^3/(t^2 + 1), t][[2]] - (t^2 + 1)");
    /* A pure-k denominator carries no t: it is absorbed, d = 1. */
    run_test("Risch`NumDen[(t + 1)/x, t][[2]]", "1");
    run_zero("Risch`NumDen[(t + 1)/x, t][[1]] - (t + 1)/x");
    /* Denominator made monic in t: den 2 t^2 + 2 -> t^2 + 1, unit pushed to num. */
    run_zero("Risch`NumDen[t/(2 t^2 + 2), t][[2]] - (t^2 + 1)");
    run_test("Coefficient[Risch`NumDen[t/(2 t^2 + 2), t][[2]], t, 2]", "1");
    run_zero("Risch`NumDen[t/(2 t^2 + 2), t][[1]] / Risch`NumDen[t/(2 t^2 + 2), t][[2]] - t/(2 t^2 + 2)");
    /* A polynomial input: denominator 1. */
    run_test("Risch`NumDen[t^2 + t + 1, t][[2]]", "1");
}

/* =====================================================================
 * Polynomial division in k[t]  (a = q b + r, deg_t r < deg_t b)
 * ===================================================================*/
static void test_poly_divide(void) {
    run_test("Risch`PolyDivide[t^3 + t + 1, t^2 + 1, t]", "List[t, 1]");
    /* Exact division: remainder 0. */
    run_test("Risch`PolyDivide[t^2 - 1, t - 1, t][[2]]", "0");
    run_zero("Risch`PolyDivide[t^2 - 1, t - 1, t][[1]] - (t + 1)");
    /* Degree-0 divisor (a field unit): q = a/b, r = 0. */
    run_zero("Risch`PolyDivide[t^2 + t, 3, t][[1]] - (t^2 + t)/3");
    run_test("Risch`PolyDivide[t^2 + t, 3, t][[2]]", "0");
    /* Over C(x): t^2 = (t - 1/x)(t + 1/x) + 1/x^2. */
    run_zero("Risch`PolyDivide[t^2, t - 1/x, t][[1]] - (t + 1/x)");
    run_zero("Risch`PolyDivide[t^2, t - 1/x, t][[2]] - 1/x^2");
    /* PROPERTY a = q b + r for a non-monic divisor over C(x). */
    run_zero("Risch`PolyDivide[x t^3 + t, x t + 1, t][[1]] (x t + 1) "
             "+ Risch`PolyDivide[x t^3 + t, x t + 1, t][[2]] - (x t^3 + t)");
    run_zero("Risch`PolyDivide[t^4 + 3 t^2 + 1, t^2 - t + 2, t][[1]] (t^2 - t + 2) "
             "+ Risch`PolyDivide[t^4 + 3 t^2 + 1, t^2 - t + 2, t][[2]] - (t^4 + 3 t^2 + 1)");
    /* Degree bound: deg_t r < deg_t b = 2. */
    run_zero("Coefficient[Risch`PolyDivide[t^4 + 3 t^2 + 1, t^2 - t + 2, t][[2]], t, 2]");
    run_zero("Coefficient[Risch`PolyDivide[t^4 + 3 t^2 + 1, t^2 - t + 2, t][[2]], t, 3]");
    /* Division by zero -> stays unevaluated. */
    run_held("Risch`PolyDivide[t, 0, t]", "Risch`PolyDivide");
}

/* =====================================================================
 * Extended Euclidean  u a + v b = g in k[t]
 * ===================================================================*/
static void test_extended_euclidean(void) {
    /* PROPERTY u a + v b = g. */
    run_zero("Risch`ExtendedEuclidean[t^2 - 1, t - 1, t][[2]] (t^2 - 1) "
             "+ Risch`ExtendedEuclidean[t^2 - 1, t - 1, t][[3]] (t - 1) "
             "- Risch`ExtendedEuclidean[t^2 - 1, t - 1, t][[1]]");
    run_zero("Risch`ExtendedEuclidean[t, t + 1, t][[2]] t "
             "+ Risch`ExtendedEuclidean[t, t + 1, t][[3]] (t + 1) "
             "- Risch`ExtendedEuclidean[t, t + 1, t][[1]]");
    run_zero("Risch`ExtendedEuclidean[t^3 - t, t^2 - 1, t][[2]] (t^3 - t) "
             "+ Risch`ExtendedEuclidean[t^3 - t, t^2 - 1, t][[3]] (t^2 - 1) "
             "- Risch`ExtendedEuclidean[t^3 - t, t^2 - 1, t][[1]]");
    /* Over C(x): u a + v b = g still holds. */
    run_zero("Risch`ExtendedEuclidean[x t^2 - x, t - 1/x, t][[2]] (x t^2 - x) "
             "+ Risch`ExtendedEuclidean[x t^2 - x, t - 1/x, t][[3]] (t - 1/x) "
             "- Risch`ExtendedEuclidean[x t^2 - x, t - 1/x, t][[1]]");
    /* Coprime a, b: gcd g is a nonzero constant (degree 0 in t). */
    run_test("Coefficient[Risch`ExtendedEuclidean[t, t + 1, t][[1]], t, 1]", "0");
    run_test("Coefficient[Risch`ExtendedEuclidean[t^2 + 1, t - 1, t][[1]], t, 1]", "0");
}

/* =====================================================================
 * Diophantine  b dn + c ds = r,  deg_t b < deg_t ds
 * ===================================================================*/
static void test_diophantine(void) {
    /* PROPERTY b dn + c ds = r. */
    run_zero("Risch`Diophantine[t - 1, t + 1, 1, t][[1]] (t - 1) "
             "+ Risch`Diophantine[t - 1, t + 1, 1, t][[2]] (t + 1) - 1");
    run_zero("Risch`Diophantine[t - 1, t + 1, t^2 + 1, t][[1]] (t - 1) "
             "+ Risch`Diophantine[t - 1, t + 1, t^2 + 1, t][[2]] (t + 1) - (t^2 + 1)");
    run_zero("Risch`Diophantine[t^2 + 1, t - 1, 3 t, t][[1]] (t^2 + 1) "
             "+ Risch`Diophantine[t^2 + 1, t - 1, 3 t, t][[2]] (t - 1) - 3 t");
    /* Over C(x). */
    run_zero("Risch`Diophantine[t - 1/x, t + 1/x, 1, t][[1]] (t - 1/x) "
             "+ Risch`Diophantine[t - 1/x, t + 1/x, 1, t][[2]] (t + 1/x) - 1");
    /* Degree bound deg_t b < deg_t ds: with ds = t + 1 (deg 1), b is constant. */
    run_zero("Coefficient[Risch`Diophantine[t - 1, t + 1, t^2 + 1, t][[1]], t, 1]");
    /* With ds = t^2 + 1 (deg 2), deg_t b <= 1. */
    run_zero("Coefficient[Risch`Diophantine[t - 1, t^2 + 1, t^3, t][[1]], t, 2]");
    run_zero("Risch`Diophantine[t - 1, t^2 + 1, t^3, t][[1]] (t - 1) "
             "+ Risch`Diophantine[t - 1, t^2 + 1, t^3, t][[2]] (t^2 + 1) - t^3");
}

/* =====================================================================
 * Normal / special classification (Def. 3.4.2) across all monomial kinds
 * ===================================================================*/
static void test_normal_special(void) {
    /* Exponential t' = t: special polynomials are exactly the powers of t. */
    run_test("Risch`SpecialQ[t, t, " EXPD "]", "True");
    run_test("Risch`SpecialQ[t^2, t, " EXPD "]", "True");
    run_test("Risch`SpecialQ[t^3, t, " EXPD "]", "True");
    run_test("Risch`NormalQ[t, t, " EXPD "]", "False");
    run_test("Risch`NormalQ[t^2 + 1, t, " EXPD "]", "True");
    run_test("Risch`NormalQ[t - 1, t, " EXPD "]", "True");
    run_test("Risch`SpecialQ[t^2 + 1, t, " EXPD "]", "False");
    run_test("Risch`SpecialQ[t^2 - 1, t, " EXPD "]", "False");

    /* Hypertangent t' = 1 + t^2: the only special irreducible is t^2 + 1. */
    run_test("Risch`SpecialQ[t^2 + 1, t, " TAND "]", "True");
    run_test("Risch`SpecialQ[(t^2 + 1)^2, t, " TAND "]", "True");
    run_test("Risch`NormalQ[t, t, " TAND "]", "True");
    run_test("Risch`NormalQ[t - 1, t, " TAND "]", "True");
    run_test("Risch`NormalQ[t^2 - 1, t, " TAND "]", "True");
    run_test("Risch`SpecialQ[t, t, " TAND "]", "False");
    run_test("Risch`NormalQ[t^2 + 1, t, " TAND "]", "False");

    /* Nonlinear t' = t^2:  t | Dt = t^2, so t is special;  t^2+1 is normal. */
    run_test("Risch`SpecialQ[t, t, " NLD "]", "True");
    run_test("Risch`NormalQ[t^2 + 1, t, " NLD "]", "True");
    run_test("Risch`NormalQ[t - 1, t, " NLD "]", "True");

    /* Primitive/log t' = 1/x:  nothing non-constant is special; all normal. */
    run_test("Risch`NormalQ[t, t, " LOGD "]", "True");
    run_test("Risch`NormalQ[t^2 + 1, t, " LOGD "]", "True");
    run_test("Risch`NormalQ[t - 1, t, " LOGD "]", "True");
    run_test("Risch`SpecialQ[t, t, " LOGD "]", "False");
    run_test("Risch`SpecialQ[t^2 + 1, t, " LOGD "]", "False");
}

/* =====================================================================
 * PolynomialReduce (§5.4): p = D[q] + r, deg_t r < deg_t(Dt), nonlinear t
 * ===================================================================*/
static void test_polynomial_reduce(void) {
    /* Nonlinear t' = t^2 (delta = 2):  reduce t^3 -> q = t^2/2, r = 0. */
    run_zero("Risch`PolynomialReduce[t^3, t, " NLD "][[1]] - t^2/2");
    run_test("Risch`PolynomialReduce[t^3, t, " NLD "][[2]]", "0");
    /* Reconstruction p = D[q] + r for higher degree. */
    run_zero("Risch`Derivation[Risch`PolynomialReduce[t^5 + t + x, t, " NLD "][[1]], " NLD "] "
             "+ Risch`PolynomialReduce[t^5 + t + x, t, " NLD "][[2]] - (t^5 + t + x)");
    /* Remainder degree < delta = 2. */
    run_zero("Coefficient[Risch`PolynomialReduce[t^5 + t + x, t, " NLD "][[2]], t, 2]");
    run_zero("Coefficient[Risch`PolynomialReduce[t^5 + t + x, t, " NLD "][[2]], t, 3]");

    /* Hypertangent t' = 1 + t^2 (delta = 2): reconstruction for t^4. */
    run_zero("Risch`Derivation[Risch`PolynomialReduce[t^4 + t + x, t, " TAND "][[1]], " TAND "] "
             "+ Risch`PolynomialReduce[t^4 + t + x, t, " TAND "][[2]] - (t^4 + t + x)");
    run_zero("Coefficient[Risch`PolynomialReduce[t^4 + t + x, t, " TAND "][[2]], t, 2]");

    /* Non-nonlinear monomials (delta < 2) are rejected -> unevaluated. */
    run_held("Risch`PolynomialReduce[t^2, t, " EXPD "]", "Risch`PolynomialReduce");  /* exp: delta 1 */
    run_held("Risch`PolynomialReduce[t^2, t, " LOGD "]", "Risch`PolynomialReduce");  /* log: delta 0 */
}

/* =====================================================================
 * Robustness: malformed input stays unevaluated (never crashes)
 * ===================================================================*/
static void test_robustness(void) {
    /* Wrong arity. */
    run_held("Risch`FieldGCD[t, t, t, t]", "Risch`FieldGCD");
    run_held("Risch`Derivation[t]", "Risch`Derivation");
    run_held("Risch`Diophantine[t - 1, t + 1, 1]", "Risch`Diophantine");
    /* t argument not a symbol. */
    run_held("Risch`FieldGCD[t^2, t, 5]", "Risch`FieldGCD");
    run_held("Risch`PolyDivide[t^2, t, 5]", "Risch`PolyDivide");
    /* Malformed derivation (a bare list, not a list of rules). */
    run_held("Risch`Derivation[t, {x, t}]", "Risch`Derivation");
    run_held("Risch`SplitFactor[t, t, {x, t}]", "Risch`SplitFactor");
    /* Empty derivation. */
    run_held("Risch`Derivation[t, {}]", "Risch`Derivation");
}

int main(void) {
    core_init();

    TEST(test_derivation);
    TEST(test_field_gcd);
    TEST(test_div_exact);
    TEST(test_num_den);
    TEST(test_poly_divide);
    TEST(test_extended_euclidean);
    TEST(test_diophantine);
    TEST(test_normal_special);
    TEST(test_polynomial_reduce);
    TEST(test_robustness);

    printf("All risch_field tests passed.\n");
    return 0;
}
