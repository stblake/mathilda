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
#include "numeric.h"
#include "symtab.h"
#include "test_utils.h"

#include <gmp.h>
#include <float.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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

/* ------------------------------------------------------------------------
 *  v0.6178: the exact core.  Every case below was wrong with the earlier
 *  double-precision continued fractions.
 * ---------------------------------------------------------------------- */

/* Evaluate `src` and require it to give True. */
static void assert_true(const char* src) {
    assert_eval_eq(src, "True", 0);
}

static void test_exact_tiny_and_huge_machine(void) {
    /* 1.*^-300 was 0 (the underflow guard took the zeroth convergent), so the
     * x term of 1.*^300 + 1.*^-300 x vanished */
    assert_true("With[{r = Rationalize[1.*^-300, 0]}, r > 0 && N[r] == 1.*^-300]");
    assert_true("With[{r = Rationalize[1.*^300 + 1.*^-300 x, 0]}, "
                "FreeQ[r, _Real] && !FreeQ[r, x] && N[Coefficient[r, x]] == 1.*^-300]");
    /* 1.*^20 and 1.*^300 stayed inexact: integers, held exactly */
    assert_eval_eq("Rationalize[1.*^20, 0]", "100000000000000000000", 0);
    assert_true("With[{r = Rationalize[1.*^300, 0]}, IntegerQ[r] && r == SetPrecision[1.*^300, Infinity]]");
    assert_eval_eq("Rationalize[1.*^20]", "100000000000000000000", 0);
    assert_eval_eq("Rationalize[1.*^12]", "1000000000000", 0);
    /* subnormals: nonzero, and strictly within half the (fixed) spacing
     * 2^-1074 of the value -- checked exactly, since N of a rational that
     * close to a rounding midpoint is not itself correctly rounded there */
    assert_true("With[{r = Rationalize[5.*^-324, 0]}, 1/2 < r 2^1074 < 3/2]");
    assert_true("With[{r = Rationalize[1.5*^-310, 0], k = Round[1.5*^-310/2.^-1074]}, "
                "r > 0 && Abs[r 2^1074 - k] < 1/2]");
    assert_true("With[{r = Rationalize[-1.*^-300, 0]}, r < 0 && N[r] == -1.*^-300]");
}

static void test_exact_round_trip_machine(void) {
    /* Rationalize[x, 0] rounds back to x, across magnitudes */
    assert_true("And @@ Table[With[{v = RandomReal[{-1, 1}] 10.^RandomInteger[{-300, 300}]}, "
                "N[Rationalize[v, 0]] == v], {200}]");
    /* the short forms are unchanged */
    assert_eval_eq("Rationalize[N[Pi], 0]", "245850922/78256779", 0);
    assert_eval_eq("Rationalize[0.1, 0]", "1/10", 0);
    assert_eval_eq("Rationalize[3.14159, 0]", "314159/100000", 0);
}

static void test_exact_arbitrary_precision(void) {
    /* N[Pi, 30] went through a double: accurate to 1e-17 only */
    assert_true("Abs[N[Pi - Rationalize[N[Pi, 30], 0], 50]] < 10^-29");
    assert_true("Abs[N[Pi - Rationalize[N[Pi, 60], 0], 80]] < 10^-59");
    assert_eval_eq("Rationalize[N[1/3, 30], 0]", "1/3", 0);
    assert_eval_eq("Rationalize[N[22/7, 40], 0]", "22/7", 0);
    /* default mode at the leaf's own precision: a 7-digit denominator past
     * the machine noise threshold (k_max ~ 1.7e6), so the machine value is
     * left alone and the 40-digit one is recovered */
    assert_eval_eq("Rationalize[N[1234567/7654321, 40]]", "1234567/7654321", 0);
    assert_true("Head[Rationalize[N[1234567/7654321]]] === Real");
    assert_true("Precision[Rationalize[N[Pi, 30]]] < Infinity");   /* left alone */
}

static void test_exact_tolerances(void) {
    /* a tolerance with a bignum denominator did not evaluate */
    assert_true("With[{r = Rationalize[N[Pi, 30], 10^-20]}, "
                "Head[r] === Rational && Abs[N[Pi - r, 40]] <= 10^-20]");
    assert_true("With[{r = Rationalize[Pi, 10^-30]}, "
                "Head[r] === Rational && Abs[N[Pi - r, 50]] <= 10^-30]");
    assert_true("With[{r = Rationalize[Sqrt[2], 10^-25]}, "
                "Head[r] === Rational && Abs[N[Sqrt[2] - r, 50]] <= 10^-25]");
    /* an exact bignum went through a double */
    assert_eval_eq("Rationalize[10^30 + 1, 1/2]", "1000000000000000000000000000001", 0);
    assert_eval_eq("Rationalize[(10^30 + 1)/3, 0]", "1000000000000000000000000000001/3", 0);
    assert_eval_eq("Rationalize[7/22, 1/10]", "1/3", 0);
    assert_eval_eq("Rationalize[Pi, 1/100]", "22/7", 0);
}

/* ------------------------------------------------------------------------
 *  v0.6179: N of an exact rational rounds once, to nearest, ties to even
 *  (numeric_mpq_to_double), and SetPrecision[x, Infinity] is exact.
 * ---------------------------------------------------------------------- */
#ifdef USE_MPFR
#include <mpfr.h>
/* Reference: an MPFR emulating an IEEE double exactly (53 bits, the double's
 * exponent range, subnormalised), so its rounding is IEEE round-to-nearest-
 * even, gradual underflow included. */
static double ref_ieee_double(const mpq_t q) {
    mpfr_exp_t emin = mpfr_get_emin(), emax = mpfr_get_emax();
    mpfr_set_emin(DBL_MIN_EXP - DBL_MANT_DIG + 1);   /* -1073 */
    mpfr_set_emax(DBL_MAX_EXP);                      /* 1024 */
    mpfr_t t;
    mpfr_init2(t, DBL_MANT_DIG);
    int tern = mpfr_set_q(t, q, MPFR_RNDN);
    tern = mpfr_check_range(t, tern, MPFR_RNDN);
    mpfr_subnormalize(t, tern, MPFR_RNDN);
    double d = mpfr_get_d(t, MPFR_RNDN);
    mpfr_clear(t);
    mpfr_set_emin(emin);
    mpfr_set_emax(emax);
    return d;
}

/* q = (num / den) 2^shift, built exactly. */
static void make_q(mpq_t q, const mpz_t num, const mpz_t den, long shift) {
    mpz_set(mpq_numref(q), num);
    mpz_set(mpq_denref(q), den);
    mpq_canonicalize(q);
    if (shift >= 0) mpq_mul_2exp(q, q, (mp_bitcnt_t)shift);
    else            mpq_div_2exp(q, q, (mp_bitcnt_t)(-shift));
}

static bool same_double(double a, double b) {
    return (a == b && signbit(a) == signbit(b)) || (isnan(a) && isnan(b));
}

static void test_mpq_to_double_matches_ieee(void) {
    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 20261105UL);
    mpz_t a, b;
    mpz_inits(a, b, NULL);
    mpq_t q;
    mpq_init(q);
    int bad = 0;
    /* random rationals across every exponent, overflow and underflow included */
    for (int i = 0; i < 20000 && bad < 5; i++) {
        mpz_urandomb(a, rs, 1 + (mp_bitcnt_t)(i % 200));
        mpz_add_ui(a, a, 1);
        mpz_urandomb(b, rs, 1 + (mp_bitcnt_t)((i * 7) % 150));
        mpz_add_ui(b, b, 1);
        if (i & 1) mpz_neg(a, a);
        long shift = (long)(i % 2300) - 1150;       /* 2^-1150 .. 2^1150 */
        make_q(q, a, b, shift);
        double got = numeric_mpq_to_double(q), want = ref_ieee_double(q);
        if (!same_double(got, want)) {
            bad++;
            gmp_fprintf(stderr, "  mpq->double %Qd: got %a, want %a\n", q, got, want);
        }
    }
    /* exact midpoints and their neighbours: (2k+1)/2 ulp, in the normal and
     * subnormal ranges and at the overflow edge -- ties must go to even */
    long exps[] = { -1074, -1073, -1060, -1023, -1022, -1000, -52, 0, 52, 900, 970, 971 };
    for (size_t j = 0; j < sizeof exps / sizeof exps[0]; j++) {
        for (unsigned long k = 0; k < 64; k++) {
            mpz_set_ui(a, 2 * (k * 2654435761UL % 9007199254740991UL) + 1);
            mpz_set_ui(b, 2);
            for (int nudge = -1; nudge <= 1; nudge++) {
                make_q(q, a, b, exps[j]);           /* (2m+1)/2 2^e: a midpoint */
                if (nudge) {
                    mpq_t eps;
                    mpq_init(eps);
                    mpq_set_ui(eps, 1, 1);
                    mpq_div_2exp(eps, eps, 1300);
                    if (nudge < 0) mpq_sub(q, q, eps); else mpq_add(q, q, eps);
                    mpq_clear(eps);
                }
                double got = numeric_mpq_to_double(q), want = ref_ieee_double(q);
                if (!same_double(got, want)) {
                    bad++;
                    if (bad < 8) fprintf(stderr, "  midpoint e=%ld k=%lu nudge=%d: got %a, want %a\n",
                                         exps[j], k, nudge, got, want);
                }
            }
        }
    }
    ASSERT(bad == 0);
    mpq_clear(q);
    mpz_clears(a, b, NULL);
    gmp_randclear(rs);
}
#endif

static void test_n_rounds_once(void) {
    /* the case that exposed it: just under 3/2 2^-1074 rounds DOWN to 2^-1074
     * (the 53-bit MPFR intermediate rounded it up to the midpoint, and then
     * ties-to-even took 2^-1073) */
    assert_true("N[Rationalize[5.*^-324, 0]] == 5.*^-324");
    assert_true("N[3/2 2^-1074] == 2^-1073 // N");          /* tie: to even */
    assert_true("N[5/2 2^-1074] == N[2 2^-1074]");          /* tie: to even */
    /* Rationalize[x, 0] round-trips at every magnitude, subnormals included */
    assert_true("And @@ Table[With[{v = RandomReal[] 2.^RandomInteger[{-1074, 1023}]}, "
                "N[Rationalize[v, 0]] == v], {500}]");
    /* int64 parts past 2^53: one rounding, not two */
    assert_true("N[(2^62 + 1)/3] == N[(2^62 + 1)/3, 40]");
    assert_true("N[(2^62 + 1)/3] == 1537228672809129301/1 // N");
    assert_true("N[1/(2^60 + 1)] == N[1/(2^60 + 1), 40]");
    /* other exact -> machine coercions round the same way as N */
    assert_true("And @@ Table[BernoulliB[1. k] === N[BernoulliB[k]], {k, 2, 60, 2}]");
}

static void test_setprecision_infinity_exact(void) {
    assert_eval_eq("SetPrecision[0.1, Infinity]", "3602879701896397/36028797018963968", 0);
    assert_eval_eq("SetPrecision[0.5, Infinity]", "1/2", 0);
    assert_true("SetPrecision[{0.5, Pi, x + 0.25, 1/3, 2.0 + 0.5 I}, Infinity] === "
                "{1/2, Pi, 1/4 + x, 1/3, 2 + I/2}");
    assert_eval_eq("SetAccuracy[0.25, Infinity]", "1/4", 0);
    assert_true("With[{r = SetPrecision[1.*^300, Infinity]}, IntegerQ[r] && N[r] === 1.*^300]");
    assert_true("SetPrecision[5.*^-324, Infinity] === 2^-1074");
    /* MPFR: the exact value of the binary number, which N gives back */
    assert_true("With[{r = SetPrecision[N[1/3, 30], Infinity]}, "
                "Head[r] === Rational && r != 1/3 && N[r, 30] == N[1/3, 30]]");
    /* exact round trip of machine reals */
    assert_true("And @@ Table[With[{v = RandomReal[{-1, 1}] 10.^RandomInteger[{-300, 300}]}, "
                "N[SetPrecision[v, Infinity]] === v], {300}]");
    assert_eval_eq("Hold[0.1] // SetPrecision[#, Infinity] &", "Hold[0.1]", 0);
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
    TEST(test_exact_tiny_and_huge_machine);
    TEST(test_exact_round_trip_machine);
    TEST(test_exact_arbitrary_precision);
    TEST(test_exact_tolerances);
#ifdef USE_MPFR
    TEST(test_mpq_to_double_matches_ieee);
#endif
    TEST(test_n_rounds_once);
    TEST(test_setprecision_infinity_exact);

    printf("All rationalize_tests passed.\n");
    return 0;
}
