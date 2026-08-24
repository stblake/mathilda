/* Tests for N[] on exact arguments that a double cannot hold.
 *
 * `N` numericalizes leaves and re-evaluates, so the leaf conversion is the
 * one step that must never be what loses information. Rounding an exact
 * argument to a double costs a relative 2^-53, and a function amplifies
 * that by its condition number — |x| for the trig family. The regression
 * this file guards is
 *
 *     N[Sin[3141592653589793238]]  ->  -0.641653   (wrong: sin of the
 *                                                   nearest double,
 *                                                   3141592653589793280)
 *                                  ->  -0.446315   (right)
 *
 * with Sin[3141592653589793238.] always having been right, because the
 * parser builds that literal as a 62-bit MPFR straight from its text.
 *
 * Everything here is checked against an MPFR oracle computed in-test from
 * the exact rational, at a precision far above the answer's — the same
 * shape as test_besselj.c's mpfr_jn oracle. Printed machine reals carry
 * only ~6 significant digits, so values are compared numerically rather
 * than as strings wherever the digits matter.
 *
 * Companion: test_numeric_stress.c sweeps the same failure mode across
 * every elementary and special function.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include <gmp.h>

/* Oracle precision: far enough above machine that a correctly-rounded
 * double falls out unambiguously even after the argument reduction that
 * Sin[10^40] needs. */
#define ORACLE_BITS 512

/* ------------------------------------------------------------------------
 *  Harness
 * ---------------------------------------------------------------------- */

static Expr* eval_str(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT_MSG(parsed != NULL, "parse failed: %s", input);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    ASSERT_MSG(r != NULL, "evaluate returned NULL: %s", input);
    return r;
}

static char* eval_to_string(const char* input) {
    Expr* r = eval_str(input);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

static void assert_prints(const char* input, const char* expected) {
    char* s = eval_to_string(input);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "%s\n  expected: %s\n  actual:   %s", input, expected, s);
    free(s);
}

static void assert_prints_prefix(const char* input, const char* prefix) {
    char* s = eval_to_string(input);
    ASSERT_MSG(strncmp(s, prefix, strlen(prefix)) == 0,
               "%s\n  expected prefix: %s\n  actual:          %s",
               input, prefix, s);
    free(s);
}

#ifdef USE_MPFR

/* Read a fully-evaluated numeric result into `out`. Deliberately accepts
 * EXPR_MPFR as well as EXPR_REAL: a machine number carries a 53-bit
 * mantissa but an arbitrary exponent, so N[Exp[1000]] is represented as a
 * DBL_MANT_DIG-bit MPFR rather than an IEEE double, exactly as N[1001!]
 * already was. */
static bool result_to_mpfr(const Expr* e, mpfr_t out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_REAL:    mpfr_set_d(out, e->data.real, MPFR_RNDN);   return true;
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z(out, e->data.bigint, MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set(out, e->data.mpfr, MPFR_RNDN);     return true;
        default:           return false;
    }
}

static void eval_to_mpfr(const char* input, mpfr_t out) {
    Expr* r = eval_str(input);
    char* s = expr_to_string(r);
    bool ok = result_to_mpfr(r, out);
    ASSERT_MSG(ok, "%s: expected a number, got %s", input, s);
    free(s);
    expr_free(r);
}

/* Exact value of a decimal integer or "num/den" rational literal. */
static void exact_of(const char* literal, mpfr_t out) {
    mpq_t q;
    mpq_init(q);
    ASSERT_MSG(mpq_set_str(q, literal, 10) == 0, "bad literal: %s", literal);
    mpq_canonicalize(q);
    mpfr_set_q(out, q, MPFR_RNDN);
    mpq_clear(q);
}

/* |got - want| <= reltol * |want|, evaluated in MPFR so the comparison
 * survives magnitudes outside IEEE range (Exp[1000] and friends). */
static void assert_close_mpfr(const char* what, const mpfr_t got,
                              const mpfr_t want, double reltol) {
    mpfr_t diff, scale;
    mpfr_inits2(ORACLE_BITS, diff, scale, (mpfr_ptr)0);
    mpfr_sub(diff, got, want, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);
    mpfr_abs(scale, want, MPFR_RNDN);
    mpfr_mul_d(scale, scale, reltol, MPFR_RNDN);
    bool ok = mpfr_lessequal_p(diff, scale) != 0;
    if (!ok) {
        char gs[64], ws[64];
        mpfr_snprintf(gs, sizeof(gs), "%.20Rg", got);
        mpfr_snprintf(ws, sizeof(ws), "%.20Rg", want);
        fprintf(stderr, "FAIL: %s\n  got:  %s\n  want: %s\n  reltol %g\n",
                what, gs, ws, reltol);
    }
    mpfr_clears(diff, scale, (mpfr_ptr)0);
    ASSERT_MSG(ok, "%s", what);
}

typedef int (*MpfrUnary)(mpfr_t, const mpfr_t, mpfr_rnd_t);

/* Every function below has a directly correctly-rounded MPFR counterpart,
 * so the oracle is exact rather than another implementation of the same
 * series. Composed inverses (ArcCot, ArcSinh, ...) and the special
 * functions MPFR does not carry are covered in test_numeric_stress.c by
 * the machine-vs-high-precision invariant instead. */
typedef struct {
    const char* head;
    MpfrUnary   ref;
    bool        needs_positive;   /* Log / Sqrt / LogGamma / PolyGamma */
} UnaryCase;

static const UnaryCase kUnary[] = {
    { "Sin",      mpfr_sin,      false },
    { "Cos",      mpfr_cos,      false },
    { "Tan",      mpfr_tan,      false },
    { "Cot",      mpfr_cot,      false },
    { "Sec",      mpfr_sec,      false },
    { "Csc",      mpfr_csc,      false },
    { "Tanh",     mpfr_tanh,     false },
    { "Coth",     mpfr_coth,     false },
    { "ArcTan",   mpfr_atan,     false },
    { "Erf",      mpfr_erf,      false },
    { "Erfc",     mpfr_erfc,     false },
    { "Log",      mpfr_log,      true  },
    { "Sqrt",     mpfr_sqrt,     true  },
    { "LogGamma", mpfr_lngamma,  true  },
    { "PolyGamma", mpfr_digamma, true  },
};
static const size_t kUnaryCount = sizeof(kUnary) / sizeof(kUnary[0]);

/* Exact arguments a double cannot hold. Written out in full decimal so the
 * oracle and the expression under test start from the same characters. */
static const char* kLossyArgs[] = {
    "9007199254740993",                          /* 2^53 + 1  */
    "3141592653589793238",                       /* the bug report */
    "10000000000000000000000000",                /* 10^25 */
    "10000000000000000000000000000000000000000", /* 10^40 */
    "100000000000000000000/3",                   /* 10^20 / 3 */
    "22/7",
};
static const size_t kLossyCount = sizeof(kLossyArgs) / sizeof(kLossyArgs[0]);

/* Arguments that ARE exactly representable. These must keep taking the
 * plain machine path — they are the control that the fix did not simply
 * route everything through MPFR. */
static const char* kExactArgs[] = {
    "9007199254740992",           /* 2^53 */
    "100000000000000000000",      /* 10^20 = 2^20 * 5^20, 5^20 < 2^53 */
    "1048576",
    "3/4",
};
static const size_t kExactCount = sizeof(kExactArgs) / sizeof(kExactArgs[0]);

/* A handful of ulps: the value under test is correctly rounded from a
 * working precision at least 64 bits above machine, the oracle from 512,
 * so they agree to the last bit except at a rounding boundary. */
#define MACHINE_RELTOL 1e-14

static void check_unary_against_oracle(const UnaryCase* c, const char* arg) {
    mpfr_t x, want, got;
    mpfr_inits2(ORACLE_BITS, x, want, got, (mpfr_ptr)0);
    exact_of(arg, x);
    if (c->needs_positive && mpfr_sgn(x) <= 0) {
        mpfr_clears(x, want, got, (mpfr_ptr)0);
        return;
    }
    c->ref(want, x, MPFR_RNDN);

    /* Skip arguments where the true value is outside IEEE range: those are
     * the overflow class, covered separately below. */
    if (!mpfr_number_p(want) || mpfr_zero_p(want)) {
        mpfr_clears(x, want, got, (mpfr_ptr)0);
        return;
    }

    char input[256];
    snprintf(input, sizeof(input), "N[%s[%s]]", c->head, arg);
    eval_to_mpfr(input, got);
    assert_close_mpfr(input, got, want, MACHINE_RELTOL);

    mpfr_clears(x, want, got, (mpfr_ptr)0);
}

/* ------------------------------------------------------------------------
 *  The reported transcript
 * ---------------------------------------------------------------------- */

static void test_reported_transcript(void) {
    /* In[6]: the real literal was always right — the parser promotes it to
     * a 62-bit MPFR from its text, never through a double. */
    assert_prints_prefix("Sin[3141592653589793238.]", "-0.4463151633593201122");

    /* In[8]: the bug. Machine N used to answer -0.641653, which is
     * Sin[3141592653589793280] — the nearest double to the argument. */
    assert_prints("N[Sin[3141592653589793238]]", "-0.446315");

    /* In[9]: the high-precision form was already right and must stay so. */
    assert_prints_prefix("N[Sin[3141592653589793238], 30]",
                         "-0.446315163359320112201603619323");

    /* The same defect at the other end of the precision range: a *low*
     * two-argument request used to round the leaf to the requested 34 bits
     * before Sin ever saw it, and answered -0.32885101403. */
    assert_prints_prefix("N[Sin[3141592653589793238], 10]", "-0.446315163");
    assert_prints_prefix("N[Sin[3141592653589793238], 5]",  "-0.44631");

    /* Sin[2^53] and Sin[2^53 + 1] must differ. They did not: both
     * arguments round to the same double. */
    Expr* a = eval_str("N[Sin[2^53]]");
    Expr* b = eval_str("N[Sin[2^53 + 1]]");
    ASSERT(a->type == EXPR_REAL && b->type == EXPR_REAL);
    ASSERT_MSG(a->data.real != b->data.real,
               "Sin[2^53] and Sin[2^53+1] collapsed to the same value (%.17g)",
               a->data.real);
    expr_free(a);
    expr_free(b);
}

/* ------------------------------------------------------------------------
 *  Elementary functions against an MPFR oracle
 * ---------------------------------------------------------------------- */

static void test_unary_lossy_args(void) {
    for (size_t i = 0; i < kUnaryCount; ++i)
        for (size_t j = 0; j < kLossyCount; ++j)
            check_unary_against_oracle(&kUnary[i], kLossyArgs[j]);
}

static void test_unary_exact_args(void) {
    for (size_t i = 0; i < kUnaryCount; ++i)
        for (size_t j = 0; j < kExactCount; ++j)
            check_unary_against_oracle(&kUnary[i], kExactArgs[j]);
}

/* The two-argument form has to hold the argument exactly no matter how few
 * digits were requested. N[Sin[10^40], 5] asked for 17 bits and used to
 * round 10^40 to 17 bits first. */
static void test_two_arg_low_precision(void) {
    static const int digits[] = { 5, 10, 15, 20, 30 };
    for (size_t j = 0; j < kLossyCount; ++j) {
        mpfr_t x, want, got;
        mpfr_inits2(ORACLE_BITS, x, want, got, (mpfr_ptr)0);
        exact_of(kLossyArgs[j], x);
        mpfr_sin(want, x, MPFR_RNDN);
        for (size_t d = 0; d < sizeof(digits) / sizeof(digits[0]); ++d) {
            char input[256];
            snprintf(input, sizeof(input), "N[Sin[%s], %d]",
                     kLossyArgs[j], digits[d]);
            eval_to_mpfr(input, got);
            /* One digit of slack for the requested-precision rounding. */
            assert_close_mpfr(input, got, want, pow(10.0, -(digits[d] - 1)));
        }
        mpfr_clears(x, want, got, (mpfr_ptr)0);
    }
}

/* Two-argument N must not manufacture precision either: the result carries
 * the digits that were asked for, not the working precision used to get
 * them. */
static void test_two_arg_precision_is_preserved(void) {
    assert_prints("Precision[N[Sin[10^25], 30]]", "30.103");
    assert_prints("Precision[N[Sin[10^25], 50]]", "50.272");
    assert_prints("Precision[N[Sin[10^25]]]", "MachinePrecision");
}

/* ------------------------------------------------------------------------
 *  Non-regression: the fast path must be untouched
 * ---------------------------------------------------------------------- */

static void test_fast_path_unchanged(void) {
    /* Bare leaves and small rationals never had a problem and must not
     * start paying for arbitrary precision. */
    assert_prints("N[1/3]",        "0.333333");
    assert_prints("N[2/7]",        "0.285714");
    assert_prints("N[Sin[1/3]]",   "0.327195");
    assert_prints("N[Sin[1]]",     "0.841471");
    assert_prints("N[Pi]",         "3.14159");
    assert_prints("N[E]",          "2.71828");
    assert_prints("N[Sqrt[2]]",    "1.41421");
    assert_prints("N[Sin[Pi/6]]",  "0.5");
    assert_prints("N[Sin[0]]",     "0.0");
    assert_prints("N[Sin[x]]",     "Sin[x]");
    assert_prints("N[Hold[1/3]]",  "Hold[1/3]");
    assert_prints("N[{1/3, 2/7, Pi}]", "{0.333333, 0.285714, 3.14159}");
    assert_prints("N[3 + 4 I]",    "3.0 + 4.0*I");
    assert_prints("N[Sqrt[-4]]",   "0.0 + 2.0*I");

    /* Exactly-representable large arguments were already correct and must
     * stay bit-identical. */
    assert_prints("N[Sin[10^20]]", "-0.645251");
    assert_prints("N[Tan[10^20]]", "-0.844602");
    assert_prints("N[10^20]",      "1e+20");
    assert_prints("N[2^53]",       "9.0072e+15");

    /* Bare N[expr] targets machine precision even for an already-approximate
     * argument: N[N[Pi, 100]] collapses to a machine number, and Sin of a
     * 100-digit Pi (tiny, ~1e-101) comes back at machine precision. */
    assert_prints("N[N[Pi, 100]]", "3.14159");
    assert_prints_prefix("N[Sin[N[Pi, 100]]]", "-8.28552e-101");
    assert_prints("Precision[N[Pi]]", "MachinePrecision");

    /* An inexact machine Real *is* its own exact binary value, so Sin of it
     * is libm's answer for that double -- deliberately NOT the answer for
     * the exact integer 10^25 alongside. Both are correct; they are
     * different questions, and that distinction is the whole point.
     *
     * The two spellings of the machine number must also agree: GMP's
     * mpz_get_d truncates toward zero where an IEEE conversion rounds to
     * nearest, which used to leave N[10^25, MachinePrecision] and 1.0*^25
     * naming doubles 2^31 apart. */
    assert_prints("Sin[1.0*^25]", "-0.305258");
    assert_prints("Sin[N[10^25, MachinePrecision]]", "-0.305258");
    assert_prints("InputForm[N[10^25, MachinePrecision] - 1.0*^25]", "0.0");
    assert_prints("N[Sin[10^25]]", "-0.74479");
}

/* ------------------------------------------------------------------------
 *  Machine numbers have an arbitrary exponent
 * ---------------------------------------------------------------------- */

static void test_overflow_keeps_a_machine_number(void) {
    struct { const char* input; const char* expr_for_oracle; } cases[] = {
        { "N[Exp[1000]]",   "exp"   },
        { "N[Exp[-1000]]",  "expneg"},
        { "N[Sinh[1000]]",  "sinh"  },
        { "N[Cosh[1000]]",  "cosh"  },
        { "N[Exp[710]]",    "exp710"},
    };
    mpfr_t want, got, x;
    mpfr_inits2(ORACLE_BITS, want, got, x, (mpfr_ptr)0);
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        const char* tag = cases[i].expr_for_oracle;
        if      (strcmp(tag, "exp")    == 0) { mpfr_set_si(x, 1000, MPFR_RNDN);  mpfr_exp(want, x, MPFR_RNDN); }
        else if (strcmp(tag, "expneg") == 0) { mpfr_set_si(x, -1000, MPFR_RNDN); mpfr_exp(want, x, MPFR_RNDN); }
        else if (strcmp(tag, "sinh")   == 0) { mpfr_set_si(x, 1000, MPFR_RNDN);  mpfr_sinh(want, x, MPFR_RNDN); }
        else if (strcmp(tag, "cosh")   == 0) { mpfr_set_si(x, 1000, MPFR_RNDN);  mpfr_cosh(want, x, MPFR_RNDN); }
        else                                 { mpfr_set_si(x, 710, MPFR_RNDN);   mpfr_exp(want, x, MPFR_RNDN); }
        eval_to_mpfr(cases[i].input, got);
        ASSERT_MSG(mpfr_number_p(got) && !mpfr_zero_p(got),
                   "%s: degenerated to inf/0", cases[i].input);
        assert_close_mpfr(cases[i].input, got, want, MACHINE_RELTOL);
    }
    mpfr_clears(want, got, x, (mpfr_ptr)0);

    /* Same 53-bit mantissa as the machine number N[1001!] has always
     * produced — this is machine precision with a wider exponent, not a
     * silent upgrade to arbitrary precision. */
    assert_prints("Precision[N[Exp[1000]]]", "15.9546");
    assert_prints("Precision[N[1001!]]",     "15.9546");

    /* Beyond MPFR's own exponent range (~10^323228458) there is nothing
     * left to recover, and the honest answer stays infinite rather than
     * costing a retry that cannot succeed. */
    assert_prints("N[Exp[10^25]]", "inf.0");
}

/* A machine number outside IEEE range must still compare. Coercing it to a
 * double saturates 4e2570 to +inf and 5e-435 to 0, so `compare_numeric`
 * used to leave the whole comparison unevaluated — N[1001!] > 10^6 has
 * never worked, and the overflow fix above makes that reachable far more
 * often (any Abs[...] < tol on a tiny difference). */
static void test_out_of_range_numbers_compare(void) {
    assert_prints("N[Exp[-1000]] < 10^-6",  "True");
    assert_prints("N[Exp[-1000]] < 1",      "True");
    assert_prints("N[Exp[-1000]] == 0",     "False");
    assert_prints("N[Exp[1000]] > 10^6",    "True");
    assert_prints("N[1001!] > 10^6",        "True");
    assert_prints("Abs[N[Exp[-1000]]] < 10^-6", "True");
    /* The exact case that surfaced it: the 13-term partial product differs
     * from its limit by ~3^-8191, which a double flushes to zero. */
    assert_prints("Abs[N[3/2 - Product[1 + (1/3)^(2^k), {k, 0, 12}]]] < 10^-6",
                  "True");

    /* In-range values keep the existing tolerant double comparison. */
    assert_prints("N[Pi, 30] > 3",          "True");
    assert_prints("N[Pi, 30] == N[Pi, 30]", "True");
    assert_prints("N[Pi, 30] > N[E, 30]",   "True");
    assert_prints("2.5`30 < 3",             "True");
}

/* ------------------------------------------------------------------------
 *  Factorial magnitude guard
 * ---------------------------------------------------------------------- */

/* N[Gamma[2^53+1]] used to abort the process: Gamma routes exact integers
 * through Factorial, whose int64 branch called mpz_fac_ui with no size
 * check, and GMP aborts rather than failing. */
static void test_factorial_magnitude_guard(void) {
    assert_prints("Head[Factorial[10^12]]",   "Factorial");
    assert_prints("Head[Factorial[2^53]]",    "Factorial");
    assert_prints("Head[Factorial[100000000]]", "Factorial");

    /* Still exact where it fits. */
    assert_prints("20!",  "2432902008176640000");
    assert_prints("Head[100!]", "Integer");
    assert_prints("1000!/999!", "1000");
    assert_prints("Head[Factorial[1000000]]", "Integer");

    /* The path that took the binary down. Any finite answer is acceptable;
     * surviving is the point. */
    Expr* r = eval_str("N[Gamma[2^53 + 1]]");
    ASSERT(r != NULL);
    expr_free(r);
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_MPFR
    /* test_utils.h arms alarm(60). This suite drives ~15 functions across
     * ~10 hard arguments at 512-bit oracle precision and needs more than
     * that on a loaded machine — still bounded so a hang cannot become a
     * hung CI job. */
    alarm(300);

    TEST(test_reported_transcript);
    TEST(test_unary_lossy_args);
    TEST(test_unary_exact_args);
    TEST(test_two_arg_low_precision);
    TEST(test_two_arg_precision_is_preserved);
    TEST(test_fast_path_unchanged);
    TEST(test_overflow_keeps_a_machine_number);
    TEST(test_out_of_range_numbers_compare);
    TEST(test_factorial_magnitude_guard);
#else
    printf("USE_MPFR is off: nothing to check (the fix is an "
           "arbitrary-precision working path).\n");
#endif

    printf("All numeric_largearg tests passed.\n");
    return 0;
}
