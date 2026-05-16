/* Tests for Rationalize.
 *
 * Cover: numeric leaf conversion, the c/q^2 default-mode bound, the
 * tolerance-mode "simplest rational in interval" semantics, dx = 0,
 * threading over Plus / Times / Power / Complex, exact-input handling,
 * symbolic pass-through, error / edge-case behaviour, and direct calls
 * to the internal_rationalize_double API.
 *
 * Comparisons against Mathematica reference outputs are pinned to the
 * ones in the project's spec discussion of Rationalize.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "rationalize.h"
#include "symtab.h"
#include "test_utils.h"

#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  internal_rationalize_double — the C-level core
 * ---------------------------------------------------------------------- */

static void test_internal_simple_terminating_decimal(void) {
    mpz_t n, d;
    bool ok = internal_rationalize_double(0.5, 0.0, RATIONALIZE_ZERO, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 1) == 0);
    ASSERT(mpz_cmp_si(d, 2) == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_default_finds_675(void) {
    mpz_t n, d;
    bool ok = internal_rationalize_double(6.75, 0.0, RATIONALIZE_DEFAULT, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 27) == 0);
    ASSERT(mpz_cmp_si(d, 4) == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_default_rejects_pi(void) {
    /* No convergent of N[Pi] satisfies |p/q - Pi| < 10^-4 / q^2; the
     * core must report failure (caller leaves x alone). */
    mpz_t n, d;
    bool ok = internal_rationalize_double(M_PI, 0.0, RATIONALIZE_DEFAULT, n, d);
    ASSERT(!ok);
    /* On failure, n/d are not initialised; nothing to clear. */
}

static void test_internal_tolerance_pi(void) {
    mpz_t n, d;
    bool ok = internal_rationalize_double(M_PI, 0.01, RATIONALIZE_TOLERANCE, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 22) == 0);
    ASSERT(mpz_cmp_si(d, 7)  == 0);
    mpz_clear(n); mpz_clear(d);

    ok = internal_rationalize_double(M_PI, 0.001, RATIONALIZE_TOLERANCE, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 201) == 0);
    ASSERT(mpz_cmp_si(d, 64)  == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_negative_tolerance_input(void) {
    /* Negative dx is treated as |dx| internally. */
    mpz_t n, d;
    bool ok = internal_rationalize_double(M_PI, -0.01, RATIONALIZE_TOLERANCE, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 22) == 0);
    ASSERT(mpz_cmp_si(d, 7)  == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_zero_input(void) {
    mpz_t n, d;
    bool ok = internal_rationalize_double(0.0, 0.0, RATIONALIZE_ZERO, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, 0) == 0);
    ASSERT(mpz_cmp_si(d, 1) == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_negative_value(void) {
    mpz_t n, d;
    bool ok = internal_rationalize_double(-6.75, 0.0, RATIONALIZE_DEFAULT, n, d);
    ASSERT(ok);
    ASSERT(mpz_cmp_si(n, -27) == 0);
    ASSERT(mpz_cmp_si(d, 4)   == 0);
    mpz_clear(n); mpz_clear(d);
}

static void test_internal_infinite_input(void) {
    mpz_t n, d;
    ASSERT(!internal_rationalize_double(INFINITY,  0.0, RATIONALIZE_DEFAULT, n, d));
    ASSERT(!internal_rationalize_double(-INFINITY, 0.0, RATIONALIZE_DEFAULT, n, d));
    ASSERT(!internal_rationalize_double(NAN,       0.0, RATIONALIZE_DEFAULT, n, d));
}

/* ------------------------------------------------------------------------
 *  Rationalize[x] — default mode, threading
 * ---------------------------------------------------------------------- */

static void test_default_simple_rational(void) {
    /* 6.75 has the smallest-denominator approximation 27/4 well within
     * the c/q^2 = 10^-4 / 16 ≈ 6.25e-6 bound. */
    assert_eval_eq("Rationalize[6.75]", "27/4", 0);
}

static void test_default_pi_unchanged(void) {
    /* No convergent of N[Pi] satisfies the bound, so the original Real
     * is returned unchanged. */
    assert_eval_eq("Rationalize[N[Pi]]", "3.14159", 0);
}

static void test_default_thread_over_plus(void) {
    /* Mathematica: Rationalize[1.2 + 6.7 x] == 6/5 + (67 x)/10
     * Mathilda prints rational coefficients as `n/d x` rather than
     * factoring x into the numerator — same value, different
     * presentation. */
    assert_eval_eq("Rationalize[1.2 + 6.7 x]", "6/5 + 67/10 x", 0);
}

static void test_default_passes_through_exact(void) {
    /* Exact inputs are returned unchanged in default mode. */
    assert_eval_eq("Rationalize[5]",     "5",   0);
    assert_eval_eq("Rationalize[3/4]",   "3/4", 0);
    assert_eval_eq("Rationalize[Pi]",    "Pi",  0);
    assert_eval_eq("Rationalize[Sqrt[2]]", "Sqrt[2]", 0);
    assert_eval_eq("Rationalize[x]",     "x",   0);
    assert_eval_eq("Rationalize[x + y]", "x + y", 0);
}

static void test_default_negative_real(void) {
    assert_eval_eq("Rationalize[-6.75]", "-27/4", 0);
}

static void test_default_threads_over_complex(void) {
    /* Default mode: only inexact components rationalise. With both
     * components inexact and within the c/q^2 bound. */
    assert_eval_eq("Rationalize[Complex[6.75, 1.5]]",
                   "Complex[Rational[27, 4], Rational[3, 2]]", 1);
}

/* ------------------------------------------------------------------------
 *  Rationalize[x, dx] — explicit positive tolerance
 * ---------------------------------------------------------------------- */

static void test_tolerance_pi(void) {
    /* 22/7 is the smallest-denominator approximation within 0.01 of Pi. */
    assert_eval_eq("Rationalize[Pi, 0.01]", "22/7", 0);
    /* And 201/64 within 0.001. */
    assert_eval_eq("Rationalize[Pi, 0.001]", "201/64", 0);
}

static void test_tolerance_works_on_exact(void) {
    /* Rationalize[Exp[Sqrt[2]], 2^-12] -> 218/53 in Mathematica. */
    assert_eval_eq("Rationalize[Exp[Sqrt[2]], 2^-12]", "218/53", 0);
}

static void test_tolerance_thread_over_polynomial(void) {
    /* Mathematica:
     *   Rationalize[2.718281828459045 x + 7.38905609893065 x^2
     *               + 20.085536923187668 x^3 + 54.598150033144236 x^4
     *               + 148.4131591025766 x^5 + 403.4287934927351 x^6, 0.1]
     * = (8 x)/3 + (22 x^2)/3 + 20 x^3 + (109 x^4)/2 + (297 x^5)/2 + (807 x^6)/2
     */
    /* Mathilda prints rational coefficients of x as `(n/d) x`, so the
     * displayed form differs from Mathematica's preference for
     * `(n x)/d` — the values are identical. */
    assert_eval_eq(
        "Rationalize[2.718281828459045 x + 7.38905609893065 x^2 + "
        "20.085536923187668 x^3 + 54.598150033144236 x^4 + "
        "148.4131591025766 x^5 + 403.4287934927351 x^6, 0.1]",
        "8/3 x + 22/3 x^2 + 20 x^3 + 109/2 x^4 + 297/2 x^5 + 807/2 x^6",
        0);
}

static void test_tolerance_zero_in_interval(void) {
    /* If 0 lies in [x - dx, x + dx] the simplest answer is 0. */
    assert_eval_eq("Rationalize[0.05, 0.1]", "0", 0);
    assert_eval_eq("Rationalize[-0.05, 0.1]", "0", 0);
}

static void test_tolerance_negative_value(void) {
    assert_eval_eq("Rationalize[-3.14, 0.01]", "-22/7", 0);
}

/* ------------------------------------------------------------------------
 *  Rationalize[x, 0] — exact-precision conversion
 * ---------------------------------------------------------------------- */

static void test_zero_pi_short_form(void) {
    /* Rationalize[N[Pi], 0] -> 245850922/78256779 (Mathematica). */
    assert_eval_eq("Rationalize[N[Pi], 0]", "245850922/78256779", 0);
}

static void test_zero_simple_decimal(void) {
    /* For an exact decimal that survives float conversion, dx=0 returns
     * the obvious small-denominator value rather than the bit-exact
     * (huge-denominator) representation of the double. */
    assert_eval_eq("Rationalize[0.5, 0]",  "1/2",  0);
    assert_eval_eq("Rationalize[0.25, 0]", "1/4",  0);
    assert_eval_eq("Rationalize[6.75, 0]", "27/4", 0);
}

static void test_zero_thread_over_complex(void) {
    /* Mathematica:
     *   Rationalize[N[Pi + 33/211 E I], 0]
     *   = 245850922/78256779 + (42051801 I)/98914198
     * Mathilda prints the imaginary coefficient inline as (n/d)*I rather
     * than absorbing I into the numerator — same value, different
     * presentation. */
    assert_eval_eq("Rationalize[N[Pi + 33/211 E I], 0]",
                   "245850922/78256779 + 42051801/98914198*I", 0);
}

static void test_zero_mode_on_zero(void) {
    assert_eval_eq("Rationalize[0.0, 0]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Edge cases / malformed input
 * ---------------------------------------------------------------------- */

static void test_zero_arg_count(void) {
    /* No args -> remain unevaluated. */
    assert_eval_eq("Rationalize[]", "Rationalize[]", 0);
    /* Three args -> remain unevaluated. */
    assert_eval_eq("Rationalize[1.5, 0.1, 0.01]",
                   "Rationalize[1.5, 0.1, 0.01]", 0);
}

static void test_symbolic_tolerance_holds(void) {
    /* Non-numeric dx -> remain unevaluated. */
    assert_eval_eq("Rationalize[1.5, x]", "Rationalize[1.5, x]", 0);
}

static void test_negative_tolerance_holds(void) {
    /* Negative tolerance is invalid -> remain unevaluated. */
    assert_eval_eq("Rationalize[1.5, -0.1]",
                   "Rationalize[1.5, -0.1]", 0);
}

static void test_rationalize_listable_via_threading(void) {
    /* Rationalize is *not* Listable, but it threads through compound
     * expressions including List, so each element is still rationalised. */
    assert_eval_eq("Rationalize[{1.2, 6.7, 0.5}]",
                   "{6/5, 67/10, 1/2}", 0);
}

static void test_rationalize_protected(void) {
    /* Built-in symbol cannot be redefined. */
    struct Expr* p = parse_expression("Rationalize[x_] := bogus");
    struct Expr* r = evaluate(p);
    expr_free(p); expr_free(r);
    /* Should still behave like the built-in. */
    assert_eval_eq("Rationalize[6.75]", "27/4", 0);
}

static void test_rationalize_rational_exact_in_zero_mode(void) {
    /* dx=0 on an exact integer / rational just returns it canonicalised
     * (it's already a rational). */
    assert_eval_eq("Rationalize[3, 0]",   "3",   0);
    assert_eval_eq("Rationalize[3/4, 0]", "3/4", 0);
}

static void test_rationalize_inside_list_with_constants(void) {
    /* Mix of inexact, exact, and constants. Default mode only converts
     * the inexact entries. */
    assert_eval_eq("Rationalize[{1.2, Pi, 6.7, x}]",
                   "{6/5, Pi, 67/10, x}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Internal API */
    TEST(test_internal_simple_terminating_decimal);
    TEST(test_internal_default_finds_675);
    TEST(test_internal_default_rejects_pi);
    TEST(test_internal_tolerance_pi);
    TEST(test_internal_negative_tolerance_input);
    TEST(test_internal_zero_input);
    TEST(test_internal_negative_value);
    TEST(test_internal_infinite_input);

    /* Default mode */
    TEST(test_default_simple_rational);
    TEST(test_default_pi_unchanged);
    TEST(test_default_thread_over_plus);
    TEST(test_default_passes_through_exact);
    TEST(test_default_negative_real);
    TEST(test_default_threads_over_complex);

    /* Tolerance mode */
    TEST(test_tolerance_pi);
    TEST(test_tolerance_works_on_exact);
    TEST(test_tolerance_thread_over_polynomial);
    TEST(test_tolerance_zero_in_interval);
    TEST(test_tolerance_negative_value);

    /* Zero mode */
    TEST(test_zero_pi_short_form);
    TEST(test_zero_simple_decimal);
    TEST(test_zero_thread_over_complex);
    TEST(test_zero_mode_on_zero);

    /* Edge cases */
    TEST(test_zero_arg_count);
    TEST(test_symbolic_tolerance_holds);
    TEST(test_negative_tolerance_holds);
    TEST(test_rationalize_listable_via_threading);
    TEST(test_rationalize_protected);
    TEST(test_rationalize_rational_exact_in_zero_mode);
    TEST(test_rationalize_inside_list_with_constants);

    printf("All rationalize_tests passed.\n");
    return 0;
}
