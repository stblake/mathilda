/*
 * test_contfrac.c — unit tests for ContinuedFraction[x] / [x, n].
 *
 * Covers exact rationals, quadratic surds (Sqrt[d]), inexact machine and
 * MPFR reals, adaptive numeric expansion of exact symbolic constants, the
 * Listable attribute, and unevaluated edge cases.  Precision-sensitive
 * cases (inexact reals, numeric-of-exact) are guarded by USE_MPFR so the
 * MPFR-free build still passes.
 */
#include "core.h"
#include "symtab.h"
#include "test_utils.h"

/* The shared assert_eval_eq() helper relies on libc assert(), which the
 * Release test build (-DNDEBUG) compiles away — failures would print but not
 * fail the run.  cf_check() parses + evaluates `input` and aborts via exit(1)
 * on mismatch, matching the ASSERT* macros in test_utils.h. */
static void cf_check(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    struct Expr* ev = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(ev);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(ev);
}

/* ------------------------------------------------------------------ */
/* Exact rationals.                                                    */
/* ------------------------------------------------------------------ */
static void test_integers(void) {
    cf_check("ContinuedFraction[5]", "{5}");
    cf_check("ContinuedFraction[0]", "{0}");
    cf_check("ContinuedFraction[-7]", "{-7}");
    cf_check("ContinuedFraction[5, 3]", "{5}");  /* count ignored, finite */
}

static void test_rationals(void) {
    /* documented example: 47/17 = 2 + 1/(1 + 1/(3 + 1/4)) */
    cf_check("ContinuedFraction[47/17]", "{2, 1, 3, 4}");
    cf_check("ContinuedFraction[5/3]", "{1, 1, 2}");
    cf_check("ContinuedFraction[1/2]", "{0, 2}");
    /* negatives: -5/3 = -2 + 1/3 */
    cf_check("ContinuedFraction[-5/3]", "{-2, 3}");
    /* canonical terminating form ends with a term >= 2 (never ...,k-1,1) */
    cf_check("ContinuedFraction[7/3]", "{2, 3}");
}

static void test_rational_count(void) {
    /* finite rational: count may yield fewer than n terms */
    cf_check("ContinuedFraction[47/17, 10]", "{2, 1, 3, 4}");
    cf_check("ContinuedFraction[47/17, 2]", "{2, 1}");
    cf_check("ContinuedFraction[47/17, 1]", "{2}");
}

static void test_big_rational(void) {
    /* numerator/denominator overflow int64 -> exercises GMP path */
    cf_check(
        "ContinuedFraction[123456789012345678901234567890/987654321]",
        "{124999998873437499901, 1, 1, 2, 1, 1, 4, 1, 2, 1, 1503, 2, 1, 2, 1, 5, 1, 3}");
}

/* ------------------------------------------------------------------ */
/* Quadratic surds Sqrt[d].                                            */
/* ------------------------------------------------------------------ */
static void test_sqrt_periodic(void) {
    /* documented example */
    cf_check("ContinuedFraction[Sqrt[13]]", "{3, {1, 1, 1, 1, 6}}");
    cf_check("ContinuedFraction[Sqrt[2]]", "{1, {2}}");
    cf_check("ContinuedFraction[Sqrt[3]]", "{1, {1, 2}}");
    cf_check("ContinuedFraction[Sqrt[7]]", "{2, {1, 1, 1, 4}}");
    cf_check("ContinuedFraction[Sqrt[19]]", "{4, {2, 1, 3, 1, 2, 8}}");
    cf_check("ContinuedFraction[Sqrt[61]]",
                   "{7, {1, 4, 3, 1, 2, 2, 1, 3, 4, 1, 14}}");
}

static void test_sqrt_count(void) {
    /* documented example: 20 terms of Sqrt[13] */
    cf_check("ContinuedFraction[Sqrt[13], 20]",
                   "{3, 1, 1, 1, 1, 6, 1, 1, 1, 1, 6, 1, 1, 1, 1, 6, 1, 1, 1, 1}");
    cf_check("ContinuedFraction[Sqrt[2], 10]",
                   "{1, 2, 2, 2, 2, 2, 2, 2, 2, 2}");
    cf_check("ContinuedFraction[Sqrt[13], 1]", "{3}");
    /* huge non-square: direct n-term generation must stay cheap and exact */
    cf_check(
        "ContinuedFraction[Sqrt[1000000000000000000000000000099], 4]",
        "{1000000000000000, 20202020202020, 4, 1}");
}

static void test_sqrt_too_long_unevaluated(void) {
    /* period astronomically long -> declined (left unevaluated) */
    cf_check(
        "ContinuedFraction[Sqrt[1000000000000000000000000000099]]",
        "ContinuedFraction[Sqrt[1000000000000000000000000000099]]");
}

/* ------------------------------------------------------------------ */
/* Listable threading.                                                 */
/* ------------------------------------------------------------------ */
static void test_listable(void) {
    cf_check("ContinuedFraction[{2/3, 3/4}]", "{{0, 1, 2}, {0, 1, 3}}");
    cf_check("ContinuedFraction[{Sqrt[2], Sqrt[3]}]",
                   "{{1, {2}}, {1, {1, 2}}}");
}

/* ------------------------------------------------------------------ */
/* Unevaluated / invalid cases.                                        */
/* ------------------------------------------------------------------ */
static void test_unevaluated(void) {
    /* exact non-rational, non-quadratic, no count -> stays put */
    cf_check("ContinuedFraction[Pi]", "ContinuedFraction[Pi]");
    /* non-positive / non-integer count -> stays put */
    cf_check("ContinuedFraction[Pi, 0]", "ContinuedFraction[Pi, 0]");
    cf_check("ContinuedFraction[Pi, -3]", "ContinuedFraction[Pi, -3]");
    /* wrong arity */
    cf_check("ContinuedFraction[]", "ContinuedFraction[]");
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Inexact reals — terms limited by the precision of the input.        */
/* ------------------------------------------------------------------ */
static void test_inexact_machine(void) {
    /* documented: ContinuedFraction[N[Pi]] runs out of machine precision */
    cf_check("ContinuedFraction[N[Pi]]",
                   "{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}");
    /* a count larger than the precision allows still stops early */
    cf_check("ContinuedFraction[N[Pi], 100]",
                   "{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}");
    /* exact float terminates */
    cf_check("ContinuedFraction[0.5]", "{0, 2}");
}

static void test_inexact_mpfr(void) {
    /* documented: 20-digit Pi yields more terms than machine precision */
    cf_check("ContinuedFraction[N[Pi, 20]]",
                   "{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2}");
}

/* ------------------------------------------------------------------ */
/* Adaptive numeric expansion of exact symbolic constants.             */
/* ------------------------------------------------------------------ */
static void test_numeric_of_exact(void) {
    /* documented: 20 terms of Pi */
    cf_check("ContinuedFraction[Pi, 20]",
                   "{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}");
    /* documented: the very regular CF of Sqrt[E] */
    cf_check("ContinuedFraction[Sqrt[E], 25]",
        "{1, 1, 1, 1, 5, 1, 1, 9, 1, 1, 13, 1, 1, 17, 1, 1, 21, 1, 1, 25, 1, 1, 29, 1, 1}");
    /* documented almost-integer: needs adaptive high precision */
    cf_check("ContinuedFraction[Exp[Pi Sqrt[163]], 10]",
        "{262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}");
    /* e = [2; 1, 2, 1, 1, 4, 1, 1, 6, ...] */
    cf_check("ContinuedFraction[E, 11]",
                   "{2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1}");
    /* GoldenRatio = [1; 1, 1, 1, ...] */
    cf_check("ContinuedFraction[GoldenRatio, 8]",
                   "{1, 1, 1, 1, 1, 1, 1, 1}");
}
#endif /* USE_MPFR */

/* ================================================================== */
/* FromContinuedFraction — inverse of ContinuedFraction.               */
/* ================================================================== */

/* --- Simple finite lists: exact rationals. --- */
static void test_fcf_rationals(void) {
    cf_check("FromContinuedFraction[{2, 1, 3, 4}]", "47/17");   /* doc example */
    cf_check("FromContinuedFraction[{5}]", "5");
    cf_check("FromContinuedFraction[{0, 2}]", "1/2");
    cf_check("FromContinuedFraction[{1, 1, 2}]", "5/3");
    cf_check("FromContinuedFraction[{-2, 3}]", "-5/3");
    cf_check("FromContinuedFraction[{2, 3}]", "7/3");
    /* a single integer term reproduces it */
    cf_check("FromContinuedFraction[{-7}]", "-7");
}

/* --- Round-trips against ContinuedFraction (exact rationals). --- */
static void test_fcf_rational_roundtrip(void) {
    cf_check("FromContinuedFraction[ContinuedFraction[47/17]]", "47/17");
    cf_check("FromContinuedFraction[ContinuedFraction[355/113]]", "355/113");
    cf_check("FromContinuedFraction[ContinuedFraction[-22/7]]", "-22/7");
    /* big rational survives the GMP path (round-trip yields lowest terms) */
    cf_check(
        "FromContinuedFraction[ContinuedFraction["
        "123456789012345678901234567890/987654321]]",
        "13717421001371742100137174210/109739369");
}

/* --- Empty / degenerate / unevaluated forms. --- */
static void test_fcf_edge(void) {
    cf_check("FromContinuedFraction[{}]", "0");
    /* not a list -> unevaluated */
    cf_check("FromContinuedFraction[5]", "FromContinuedFraction[5]");
    /* wrong arity -> unevaluated */
    cf_check("FromContinuedFraction[]", "FromContinuedFraction[]");
    cf_check("FromContinuedFraction[{1, 2}, 3]",
             "FromContinuedFraction[{1, 2}, 3]");
    /* a sub-list anywhere but last is invalid -> unevaluated */
    cf_check("FromContinuedFraction[{1, {2, 3}, 4}]",
             "FromContinuedFraction[{1, {2, 3}, 4}]");
    /* nested list inside the period block -> non-integer terms -> declined */
    cf_check("FromContinuedFraction[{1, {2, 3, {4}}}]",
             "FromContinuedFraction[{1, {2, 3, {4}}}]");
    /* empty period block -> declined */
    cf_check("FromContinuedFraction[{2, {}}]",
             "FromContinuedFraction[{2, {}}]");
}

/* --- Symbolic terms give the nested convergent form. --- */
static void test_fcf_symbolic(void) {
    cf_check("FromContinuedFraction[{x}]", "x");
    cf_check("FromContinuedFraction[{a, b}]", "(1 + a b)/b");
    /* documented convergent form (un-expanded) */
    cf_check("FromContinuedFraction[{a, b, c, d}]",
             "(1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)");
    /* Together collapses it to the fully expanded rational */
    cf_check("Together[FromContinuedFraction[{a, b, c, d}]]",
             "(1 + a b + a d + c d + a b c d)/(b + d + b c d)");
}

/* --- Purely periodic blocks -> quadratic irrationals (doc examples). --- */
static void test_fcf_periodic_pure(void) {
    cf_check("FromContinuedFraction[{{1}}]", "1/2 (1 + Sqrt[5])");
    cf_check("FromContinuedFraction[{{1, 2}}]", "1/2 (1 + Sqrt[3])");
    cf_check("FromContinuedFraction[{{1, 2, 3}}]", "1/7 (4 + Sqrt[37])");
    cf_check("FromContinuedFraction[{{1, 2, 3, 4}}]", "1/15 (9 + 2 Sqrt[39])");
    /* purely periodic [2;2,2,...] = 1 + Sqrt[2] */
    cf_check("FromContinuedFraction[{{2}}]", "1 + Sqrt[2]");
}

/* --- Eventually-periodic blocks and Sqrt round-trips. --- */
static void test_fcf_periodic_lead(void) {
    /* documented: the CF of Sqrt[71] reconstructs exactly */
    cf_check("FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]",
             "Sqrt[71]");
    cf_check("FromContinuedFraction[{3, {6}}]", "Sqrt[10]");
    cf_check("FromContinuedFraction[{1, {2}}]", "Sqrt[2]");
    /* round-trips through ContinuedFraction for several surds */
    cf_check("FromContinuedFraction[ContinuedFraction[Sqrt[2]]]", "Sqrt[2]");
    cf_check("FromContinuedFraction[ContinuedFraction[Sqrt[13]]]", "Sqrt[13]");
    cf_check("FromContinuedFraction[ContinuedFraction[Sqrt[61]]]", "Sqrt[61]");
    cf_check("FromContinuedFraction[ContinuedFraction[Sqrt[71]]]", "Sqrt[71]");
}

/* --- FromContinuedFraction is NOT Listable (acts on the list as a whole). */
static void test_fcf_not_listable(void) {
    /* a list of two CF lists would thread if Listable; instead the outer list
     * is read as a single (here invalid, non-integer-prefixed) CF and the
     * trailing sub-list marks a period -> evaluates as one quadratic value. */
    cf_check("FromContinuedFraction[{1, 1, 1, 1, 1, 1}]", "13/8");
}

#ifdef USE_MPFR
/* --- Rational approximation: matches the documented N value. --- */
static void test_fcf_pi_approx(void) {
    cf_check("FromContinuedFraction[ContinuedFraction[Pi, 3]]", "333/106");
    cf_check("FromContinuedFraction[ContinuedFraction[Pi, 6]]", "104348/33215");
}
#endif

int main(void) {
    symtab_init();
    core_init();

    TEST(test_integers);
    TEST(test_rationals);
    TEST(test_rational_count);
    TEST(test_big_rational);
    TEST(test_sqrt_periodic);
    TEST(test_sqrt_count);
    TEST(test_sqrt_too_long_unevaluated);
    TEST(test_listable);
    TEST(test_unevaluated);
#ifdef USE_MPFR
    TEST(test_inexact_machine);
    TEST(test_inexact_mpfr);
    TEST(test_numeric_of_exact);
#endif

    /* FromContinuedFraction (inverse). */
    TEST(test_fcf_rationals);
    TEST(test_fcf_rational_roundtrip);
    TEST(test_fcf_edge);
    TEST(test_fcf_symbolic);
    TEST(test_fcf_periodic_pure);
    TEST(test_fcf_periodic_lead);
    TEST(test_fcf_not_listable);
#ifdef USE_MPFR
    TEST(test_fcf_pi_approx);
#endif

    printf("All contfrac_tests passed.\n");
    return 0;
}
