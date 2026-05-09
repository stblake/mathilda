/*
 * test_radical_canonical.c -- Coverage for the radical canonicalization
 * passes that fold integer-base Power factors with rational coefficients
 * into a single canonical Power form.
 *
 * Two cooperating rules are exercised here:
 *
 *  1. Power.c integer-part extraction. For Power[n, p/q] with n a positive
 *     integer, q > 1, and p >= q, write p = a*q + b (0 <= b < q) so
 *     n^(p/q) = n^a * n^(b/q). Combined with the existing perfect-q-th-
 *     power reduction (n = m^q * r), the coefficient is n^a * m^b and the
 *     residual is r^(b/q).
 *
 *  2. Times.c radical canonicalization. After base-grouping in
 *     builtin_times, each Power[b, q] group with b a positive integer and
 *     q rational pulls factors of b out of the rational coefficient
 *     num_prod -- always from the denominator (so den ends up coprime to
 *     b), and from the numerator only enough to lift a still-negative
 *     exponent up to >= 0 (so we never hand Power.c an exponent > 1 that
 *     would re-extract and loop).
 *
 * Together these turn
 *      Sqrt[2]/2          ->  1/Sqrt[2]              (Power[2, -1/2])
 *      2^(1/3)/2          ->  1/2^(2/3)              (Power[2, -2/3])
 *      3^(3/2)            ->  3 Sqrt[3]
 *      3^(7/2)            ->  27 Sqrt[3]
 *      8/Sqrt[2]          ->  4 Sqrt[2]
 *      2/Sqrt[2]          ->  Sqrt[2]
 * into their canonical Mathematica-compatible forms.
 *
 * The tests check both the user-facing display string and FullForm where
 * the structural representation matters; regression cases ensure that
 * neighbouring rules (radical fusion across different bases, Sqrt
 * reduction, Plus collection of like-radicals) keep working alongside
 * the new canonicalization.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>

/* Strong-assert wrapper. The shared test_utils.h `assert_eval_eq` uses
 * libc assert(), which CMake's Release builds silence via NDEBUG -- a
 * mismatch only prints FAIL and the run still exits 0, masking failures
 * from CTest. Route through this wrapper so each mismatch is a hard
 * exit(1). Pattern lifted from test_radical_simplify.c. */
static void check_eval_eq(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL: parse failure for: %s\n", input);
        exit(1);
    }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
        free(str);
        expr_free(evaluated);
        exit(1);
    }
    free(str);
    expr_free(evaluated);
}
static void check_eval_fullform_eq(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: parse failure for: %s\n", input); exit(1); }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string_fullform(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s (FullForm)\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
        free(str); expr_free(evaluated); exit(1);
    }
    free(str); expr_free(evaluated);
}
#define assert_eval_eq(in, exp, ff) \
    do { if (ff) check_eval_fullform_eq((in), (exp)); else check_eval_eq((in), (exp)); } while (0)

/* ------------------------------------------------------------------ */
/* User-supplied motivating examples                                    */
/* ------------------------------------------------------------------ */

static void test_user_sqrt_over_two(void) {
    /* Sqrt[2]/2 -> 1/Sqrt[2]. The rational coefficient 1/2's denominator
     * absorbs into the Sqrt[2] exponent, giving Power[2, -1/2]. */
    assert_eval_eq("Sqrt[2]/2", "1/Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2]/2", "Power[2, Rational[-1, 2]]", 1);
    /* Same value via different routes. */
    assert_eval_eq("1/2 * Sqrt[2]",        "1/Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2] * 1/2",        "1/Sqrt[2]", 0);
    assert_eval_eq("Power[2, 1/2] / 2",    "1/Sqrt[2]", 0);
    assert_eval_eq("Power[2, 1/2 - 1]",    "1/Sqrt[2]", 0);
    assert_eval_eq("2^(1/2)/2",            "1/Sqrt[2]", 0);
}

static void test_user_cube_root_over_two(void) {
    /* 2^(1/3)/2 -> 1/2^(2/3). p = 1, ed = 3, k_d = 1 -> new exp 1/3 - 1
     * = -2/3. */
    assert_eval_eq("2^(1/3)/2",            "1/2^(2/3)", 0);
    assert_eval_eq("2^(1/3)/2",            "Power[2, Rational[-2, 3]]", 1);
    assert_eval_eq("Power[2, 1/3] / 2",    "1/2^(2/3)", 0);
    assert_eval_eq("Power[2, 1/3 - 1]",    "1/2^(2/3)", 0);
}

static void test_user_radical_sum_collects(void) {
    /* 3^(3/2) + 3^(5/2) + 3^(7/2) -> 39 Sqrt[3].
     *   3^(3/2) = 3 Sqrt[3], 3^(5/2) = 9 Sqrt[3], 3^(7/2) = 27 Sqrt[3]. */
    assert_eval_eq("3^(3/2) + 3^(5/2) + 3^(7/2)", "39 Sqrt[3]", 0);
    /* Order shouldn't matter -- Plus is Orderless. */
    assert_eval_eq("3^(7/2) + 3^(5/2) + 3^(3/2)", "39 Sqrt[3]", 0);
    /* Mixed with other terms. */
    assert_eval_eq("3^(3/2) + 3^(5/2)",                "12 Sqrt[3]", 0);
    assert_eval_eq("3^(3/2) + 3^(5/2) + 1",            "1 + 12 Sqrt[3]", 0);
    assert_eval_eq("x + 3^(3/2) + 3^(5/2)",            "12 Sqrt[3] + x", 0);
}

/* ------------------------------------------------------------------ */
/* Power.c integer-part extraction                                     */
/* ------------------------------------------------------------------ */

static void test_power_extracts_unit_residue(void) {
    /* m == 1 (n is q-th-power-free): n^(p/q) -> n^a * Power[n, b/q]. */
    assert_eval_eq("3^(3/2)", "3 Sqrt[3]", 0);
    assert_eval_eq("3^(5/2)", "9 Sqrt[3]", 0);
    assert_eval_eq("3^(7/2)", "27 Sqrt[3]", 0);
    assert_eval_eq("5^(3/2)", "5 Sqrt[5]", 0);
    assert_eval_eq("7^(5/2)", "49 Sqrt[7]", 0);
    /* Cube roots. */
    assert_eval_eq("2^(4/3)", "2 2^(1/3)", 0);
    assert_eval_eq("2^(5/3)", "2 2^(2/3)", 0);
    assert_eval_eq("2^(7/3)", "4 2^(1/3)", 0);
    assert_eval_eq("3^(4/3)", "3 3^(1/3)", 0);
}

static void test_power_extracts_with_perfect_factor(void) {
    /* m > 1 (n has q-th-power factors): n^(p/q) = n^a * m^b * r^(b/q).
     *
     * 12 = 4 * 3 = 2^2 * 3, so for q=2: m=2, r=3.
     * 12^(3/2): a=1, b=1. coeff = 12^1 * 2^1 = 24. residue = Sqrt[3]. */
    assert_eval_eq("12^(3/2)", "24 Sqrt[3]", 0);
    /* 12^(5/2): a=2, b=1. coeff = 144 * 2 = 288. residue = Sqrt[3]. */
    assert_eval_eq("12^(5/2)", "288 Sqrt[3]", 0);
    /* 18 = 9 * 2 = 3^2 * 2, q=2: m=3, r=2.
     * 18^(3/2): a=1, b=1. coeff = 18 * 3 = 54. residue = Sqrt[2]. */
    assert_eval_eq("18^(3/2)", "54 Sqrt[2]", 0);
    /* 24 = 8 * 3 = 2^3 * 3, q=3: m=2, r=3.
     * 24^(4/3): a=1, b=1. coeff = 24 * 2 = 48. residue = 3^(1/3). */
    assert_eval_eq("24^(4/3)", "48 3^(1/3)", 0);
}

static void test_power_perfect_power_collapses(void) {
    /* When the residual r = 1 the result reduces to a pure integer. */
    assert_eval_eq("4^(3/2)",  "8", 0);
    assert_eval_eq("9^(3/2)",  "27", 0);
    assert_eval_eq("8^(2/3)",  "4", 0);
    assert_eval_eq("27^(2/3)", "9", 0);
    assert_eval_eq("27^(4/3)", "81", 0);
    assert_eval_eq("16^(3/4)", "8", 0);
    assert_eval_eq("32^(2/5)", "4", 0);
    /* Negative integer parts of the exponent likewise collapse. */
    assert_eval_eq("4^(-3/2)", "1/8", 0);
    assert_eval_eq("8^(-2/3)", "1/4", 0);
}

static void test_power_bigint_coefficient(void) {
    /* Integer-part extraction must promote into BigInt cleanly when the
     * coefficient overflows int64. 2^(101/2) = 2^50 * Sqrt[2]; 2^50 fits
     * but 2^(201/2) = 2^100 * Sqrt[2] overflows int64. */
    assert_eval_eq("2^(101/2)", "1125899906842624 Sqrt[2]", 0);
    assert_eval_eq("2^(201/2)", "1267650600228229401496703205376 Sqrt[2]", 0);
    /* 3^(7/2) = 27 Sqrt[3] (regression: ensure a small case still folds
     * correctly after the BigInt path is enabled). */
    assert_eval_eq("3^(7/2)", "27 Sqrt[3]", 0);
}

static void test_power_negative_in_residue(void) {
    /* p < 0 with |p| < q: a = 0, b_rem = p (negative). Yields a rational
     * coefficient and a residual Power[r, b_rem/q] with negative exp. */
    /* 8^(-1/2) = 2^(-1) * 2^(-1/2) = 1/2 * 1/Sqrt[2] = 1/2^(3/2). */
    assert_eval_eq("8^(-1/2)", "1/2^(3/2)", 0);
    /* 12^(-1/2) = 12 = 2^2 * 3, m=2, r=3. coeff = 1/2, residue = 1/Sqrt[3].
     * Times then has num_prod = 1/2, group (3, -1/2): no factors of 3 in
     * num_prod. Final: Times[1/2, Power[3, -1/2]]. The default printer
     * emits this as `1/2/Sqrt[3]`. Mathematically equal to 1/(2 Sqrt[3]). */
    assert_eval_eq("12^(-1/2)", "1/2/Sqrt[3]", 0);
    /* FullForm pins the structure unambiguously. */
    assert_eval_eq("12^(-1/2)",
                   "Times[Rational[1, 2], Power[3, Rational[-1, 2]]]", 1);
    /* 27^(-2/3) = 1/9 (perfect power, residue = 1). */
    assert_eval_eq("27^(-2/3)", "1/9", 0);
}

static void test_power_irreducible_unchanged(void) {
    /* When neither integer-part extraction nor perfect-power reduction
     * applies the Power form must survive intact. */
    assert_eval_eq("Sqrt[2]",     "Sqrt[2]",     0);
    assert_eval_eq("Sqrt[3]",     "Sqrt[3]",     0);
    assert_eval_eq("Sqrt[5]",     "Sqrt[5]",     0);
    assert_eval_eq("2^(1/3)",     "2^(1/3)",     0);
    assert_eval_eq("3^(2/5)",     "3^(2/5)",     0);
    assert_eval_eq("2^(1/2)",     "Sqrt[2]",     0);
    assert_eval_eq("Power[2, -1/2]", "1/Sqrt[2]", 0);
}

static void test_power_negative_base_unchanged(void) {
    /* Negative integer base: integer-part extraction must NOT fire (sign
     * branch handling lives in builtin_power). For q == 2 the -1 still
     * factors out as I; other q values are left for Power's symbolic
     * handler. */
    /* (-1)^(1/2) = I */
    assert_eval_eq("(-1)^(1/2)", "I", 0);
    /* (-3)^(1/2) = I Sqrt[3]. (printer parenthesises the I) */
    assert_eval_eq("Power[-3, 1/2]", "(I) Sqrt[3]", 0);
    /* (-2)^(1/3) is left to the parser/printer round-trip path. */
    assert_eval_eq("Power[-2, 1/3]", "(-2)^(1/3)", 0);
}

/* ------------------------------------------------------------------ */
/* Times.c radical canonicalization                                    */
/* ------------------------------------------------------------------ */

static void test_times_pulls_from_denominator(void) {
    /* The motivating cases: factors of the Power's base sitting in the
     * coefficient's denominator absorb into the exponent. */
    assert_eval_eq("Sqrt[2]/2",      "1/Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2]/4",      "1/2^(3/2)", 0);
    assert_eval_eq("Sqrt[2]/8",      "1/2^(5/2)", 0);
    assert_eval_eq("Sqrt[3]/3",      "1/Sqrt[3]", 0);
    assert_eval_eq("Sqrt[3]/9",      "1/3^(3/2)", 0);
    assert_eval_eq("Sqrt[3]/27",     "1/3^(5/2)", 0);
    assert_eval_eq("Sqrt[5]/5",      "1/Sqrt[5]", 0);
    assert_eval_eq("2^(1/3)/2",      "1/2^(2/3)", 0);
    assert_eval_eq("2^(1/3)/4",      "1/2^(5/3)", 0);
    assert_eval_eq("3^(1/4)/3",      "1/3^(3/4)", 0);
}

static void test_times_pulls_from_numerator_when_negative_exp(void) {
    /* When the Power exponent is negative the canonicalizer pulls factors
     * of b out of the numerator, but only enough to make the new exponent
     * >= 0 (avoiding the loop with Power.c's integer-part extraction). */
    /* 2/Sqrt[2] = Sqrt[2]: pull one factor of 2, exp -1/2 -> 1/2. */
    assert_eval_eq("2/Sqrt[2]",      "Sqrt[2]", 0);
    /* 4/Sqrt[2]: pull one factor (exp goes -1/2 -> 1/2), num 4 -> 2.
     * Stop -- pulling more would drive exp >= 1 and Power.c would
     * re-emit. Final: Times[2, Sqrt[2]]. */
    assert_eval_eq("4/Sqrt[2]",      "2 Sqrt[2]", 0);
    /* 8/Sqrt[2] -> 4 Sqrt[2]. */
    assert_eval_eq("8/Sqrt[2]",      "4 Sqrt[2]", 0);
    /* 16/Sqrt[2] -> 8 Sqrt[2]. */
    assert_eval_eq("16/Sqrt[2]",     "8 Sqrt[2]", 0);
    /* 32/Sqrt[2] -> 16 Sqrt[2]. */
    assert_eval_eq("32/Sqrt[2]",     "16 Sqrt[2]", 0);
    /* Cube-root counterpart: 4/2^(1/3) = 2^2 * 2^(-1/3) = 2^(5/3) =
     * 2 * 2^(2/3). */
    assert_eval_eq("4/2^(1/3)",      "2 2^(2/3)", 0);
    /* 8/2^(1/3) = 2^3 * 2^(-1/3) = 2^(8/3) = 4 * 2^(2/3). */
    assert_eval_eq("8/2^(1/3)",      "4 2^(2/3)", 0);
}

static void test_times_does_not_pull_when_positive_exp(void) {
    /* If the Power exponent is already >= 0, factors of b in num_prod's
     * numerator stay out -- pulling them would push the exponent past 1
     * and Power's integer extraction would emit them right back. */
    assert_eval_eq("2 Sqrt[2]",      "2 Sqrt[2]", 0);
    assert_eval_eq("4 Sqrt[2]",      "4 Sqrt[2]", 0);
    assert_eval_eq("8 Sqrt[2]",      "8 Sqrt[2]", 0);
    assert_eval_eq("3 Sqrt[3]",      "3 Sqrt[3]", 0);
    assert_eval_eq("9 Sqrt[3]",      "9 Sqrt[3]", 0);
    /* Cube roots. */
    assert_eval_eq("2 * 2^(1/3)",    "2 2^(1/3)", 0);
    assert_eval_eq("4 * 2^(1/3)",    "4 2^(1/3)", 0);
    assert_eval_eq("3 * 3^(1/4)",    "3 3^(1/4)", 0);
}

static void test_times_unrelated_base_stays(void) {
    /* Coefficient factors that are coprime to the radical base do not
     * touch the Power. */
    assert_eval_eq("3/Sqrt[2]",      "3/Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2]/3",      "1/3 Sqrt[2]", 0);
    assert_eval_eq("5 Sqrt[3]",      "5 Sqrt[3]", 0);
    assert_eval_eq("Sqrt[6]/5",      "1/5 Sqrt[6]", 0);
    /* The 6 in num doesn't pull from a Sqrt[2] or Sqrt[3] group because
     * we only pull factors of the WHOLE base, not its prime factors. */
    assert_eval_eq("6 Sqrt[2]",      "6 Sqrt[2]", 0);
    assert_eval_eq("6 Sqrt[3]",      "6 Sqrt[3]", 0);
}

static void test_times_rational_coefficient(void) {
    /* Mixed numerator and denominator factors on the same base. */
    /* (3/2) * Sqrt[2] = 3 / Sqrt[2] * (?) ... actually:
     * pull 1 factor of 2 from den. num=3, den=1, exp=-1/2.
     * Then exp < 0 and num=3 has no factors of 2. Stop. -> 3 Power[2,-1/2]
     * = 3/Sqrt[2]. */
    assert_eval_eq("3/2 * Sqrt[2]",      "3/Sqrt[2]", 0);
    /* (4/3) * Sqrt[3]: pull 1 factor of 3 from den. num=4, den=1, exp=-1/2.
     * exp < 0; num=4 coprime to 3. -> 4 Power[3, -1/2] = 4/Sqrt[3]. */
    assert_eval_eq("4/3 * Sqrt[3]",      "4/Sqrt[3]", 0);
    /* (2/3) * Sqrt[2]: den=3 coprime to 2, num=2 -- but exp=1/2 is
     * positive so we don't pull from num. -> 2/3 Sqrt[2]. */
    assert_eval_eq("2/3 * Sqrt[2]",      "2/3 Sqrt[2]", 0);
}

/* ------------------------------------------------------------------ */
/* Plus collection of canonicalized radicals                           */
/* ------------------------------------------------------------------ */

static void test_plus_collects_like_radicals(void) {
    /* After canonicalization the integer-coefficient * Sqrt[r] terms have
     * a common radical factor and Plus collects them. */
    assert_eval_eq("Sqrt[2] + Sqrt[2]",          "2 Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2] + 3 Sqrt[2]",        "4 Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2] + Sqrt[8]",          "3 Sqrt[2]", 0);     /* 1 + 2 */
    assert_eval_eq("Sqrt[2] + Sqrt[18]",         "4 Sqrt[2]", 0);     /* 1 + 3 */
    assert_eval_eq("Sqrt[2] + Sqrt[8] + Sqrt[18]", "6 Sqrt[2]", 0);
    /* Higher-power cases via integer extraction. */
    assert_eval_eq("3^(3/2) + 3^(3/2)",           "6 Sqrt[3]", 0);
    assert_eval_eq("2^(5/3) + 2^(5/3)",           "4 2^(2/3)", 0);
    assert_eval_eq("2^(5/3) + 2^(2/3)",           "3 2^(2/3)", 0);
    /* Coefficient-only collection. */
    assert_eval_eq("Sqrt[2]/2 + Sqrt[2]/2",       "Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2]/2 - Sqrt[2]/2",       "0", 0);
    /* The motivating case from the user's report. */
    assert_eval_eq("Sqrt[2]/2 - 1/Sqrt[2]",       "0", 0);
}

static void test_plus_distinct_radicals_remain(void) {
    /* Different radical bases must NOT merge -- they're independent
     * algebraic generators. */
    assert_eval_eq("Sqrt[2] + Sqrt[3]",     "Sqrt[2] + Sqrt[3]", 0);
    /* Plus is Orderless and the canonical ordering puts the cube root
     * before the square root by exponent. */
    assert_eval_eq("Sqrt[2] + 2^(1/3)",     "2^(1/3) + Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2] - Sqrt[3]",     "Sqrt[2] - Sqrt[3]", 0);
}

/* ------------------------------------------------------------------ */
/* Compatibility with existing radical machinery                       */
/* ------------------------------------------------------------------ */

static void test_radical_fusion_compat(void) {
    /* The new Times canonicalization runs BEFORE the existing radical
     * fusion (a^q * b^(-q) -> (a/b)^q). The fusion rule still fires for
     * different positive integer bases. */
    assert_eval_eq("Sqrt[6]/Sqrt[2]",        "Sqrt[3]", 0);
    assert_eval_eq("Sqrt[15]/Sqrt[5]",       "Sqrt[3]", 0);
    assert_eval_eq("Sqrt[10]/Sqrt[3]",       "Sqrt[10/3]", 0);
    assert_eval_eq("Power[6, 1/3] / Power[2, 1/3]", "3^(1/3)", 0);
}

static void test_sqrt_reduction_compat(void) {
    /* Existing perfect-square extraction inside Sqrt continues to fire. */
    assert_eval_eq("Sqrt[8]",   "2 Sqrt[2]", 0);
    assert_eval_eq("Sqrt[12]",  "2 Sqrt[3]", 0);
    assert_eval_eq("Sqrt[18]",  "3 Sqrt[2]", 0);
    assert_eval_eq("Sqrt[50]",  "5 Sqrt[2]", 0);
    assert_eval_eq("Sqrt[200]", "10 Sqrt[2]", 0);
    /* Same-base cancellation from Power grouping. */
    assert_eval_eq("Sqrt[2] * Sqrt[2]", "2", 0);
    assert_eval_eq("Sqrt[3] * Sqrt[3]", "3", 0);
}

static void test_radical_with_symbolic(void) {
    /* Symbolic factors pass through; only the numeric-coefficient *
     * Power[b, q] interaction is rewritten. After Times canonicalization
     * x Sqrt[2]/2 becomes x * Power[2, -1/2] = x/Sqrt[2]. */
    assert_eval_eq("x Sqrt[2]/2",      "x/Sqrt[2]", 0);
    assert_eval_eq("x / (2 Sqrt[2])",  "x/2^(3/2)", 0);
    /* Mathematical identity preservation through Simplify. */
    assert_eval_eq("Simplify[x Sqrt[2]/2 - x/Sqrt[2]]",     "0", 0);
    assert_eval_eq("Simplify[2^(1/3) x / 2 - x/2^(2/3)]",   "0", 0);
    /* Plus collection across canonicalized symbolic terms: both terms
     * have the same Power[2,-1/2] base so they merge to 2 x / Sqrt[2],
     * and Times then absorbs the 2 / Sqrt[2] -> Sqrt[2]. */
    assert_eval_eq("x Sqrt[2]/2 + x/Sqrt[2]",  "Sqrt[2] x", 0);
}

/* ------------------------------------------------------------------ */
/* Round-trip / fixed-point stability                                  */
/* ------------------------------------------------------------------ */

static void test_canonical_form_idempotent(void) {
    /* Re-evaluating a canonical form must yield the same canonical form.
     * If the rules did loop, eval would re-extract or re-fuse and the
     * second pass would diverge. */
    assert_eval_eq("1/Sqrt[2]",       "1/Sqrt[2]", 0);
    assert_eval_eq("1/2^(2/3)",       "1/2^(2/3)", 0);
    assert_eval_eq("1/2^(3/2)",       "1/2^(3/2)", 0);
    assert_eval_eq("3 Sqrt[3]",       "3 Sqrt[3]", 0);
    assert_eval_eq("39 Sqrt[3]",      "39 Sqrt[3]", 0);
    assert_eval_eq("4 Sqrt[2]",       "4 Sqrt[2]", 0);
}

static void test_canonical_form_negative_exponent(void) {
    /* Negative-exponent Power forms stay Power[b, q] (Power.c does NOT
     * extract negative integer parts) so the Times canonicalizer
     * remains the only normaliser for these and the result is stable. */
    assert_eval_eq("Power[2, -1/2]",  "1/Sqrt[2]", 0);
    assert_eval_eq("Power[2, -3/2]",  "1/2^(3/2)", 0);
    assert_eval_eq("Power[2, -5/2]",  "1/2^(5/2)", 0);
    assert_eval_eq("Power[2, -2/3]",  "1/2^(2/3)", 0);
    assert_eval_eq("Power[3, -3/2]",  "1/3^(3/2)", 0);
    /* FullForm equivalents. */
    assert_eval_eq("Power[2, -1/2]",  "Power[2, Rational[-1, 2]]", 1);
    assert_eval_eq("Power[2, -3/2]",  "Power[2, Rational[-3, 2]]", 1);
}

static void test_round_trip_arithmetic(void) {
    /* Algebraic identities that exercise both passes simultaneously. */
    /* Sqrt[2] * Sqrt[2] = 2. */
    assert_eval_eq("Sqrt[2] Sqrt[2]",       "2", 0);
    /* (Sqrt[2])^3 = 2 Sqrt[2]. */
    assert_eval_eq("Sqrt[2]^3",             "2 Sqrt[2]", 0);
    /* (Sqrt[2])^4 = 4. */
    assert_eval_eq("Sqrt[2]^4",             "4", 0);
    /* (Sqrt[2])^(-1) = 1/Sqrt[2]. */
    assert_eval_eq("Sqrt[2]^(-1)",          "1/Sqrt[2]", 0);
    /* (1/Sqrt[2])^2 = 1/2. */
    assert_eval_eq("(1/Sqrt[2])^2",         "1/2", 0);
    /* (1/Sqrt[2])^3 = 1/2^(3/2). */
    assert_eval_eq("(1/Sqrt[2])^3",         "1/2^(3/2)", 0);
    /* Mixed: Sqrt[2] * (1/Sqrt[2]) = 1. */
    assert_eval_eq("Sqrt[2] * (1/Sqrt[2])", "1", 0);
    /* Sqrt[2]/Sqrt[2] = 1. */
    assert_eval_eq("Sqrt[2]/Sqrt[2]",       "1", 0);
}

/* ------------------------------------------------------------------ */
/* Stress / regression                                                 */
/* ------------------------------------------------------------------ */

static void test_zero_and_one_exponents_unchanged(void) {
    /* Power[b, 0] and Power[b, 1] short-circuit before reaching the
     * radical-canonical paths; ensure they still do. */
    assert_eval_eq("2^0",      "1", 0);
    assert_eval_eq("2^1",      "2", 0);
    assert_eval_eq("0^0",      "1", 0);
    assert_eval_eq("Power[2, 1]", "2", 0);
    assert_eval_eq("Power[3, 0]", "1", 0);
}

static void test_real_coefficient_skips(void) {
    /* Inexact Real coefficients trigger numeric contagion (Sqrt[2]
     * numericalises) so they bypass the canonicalization entirely. */
    assert_eval_startswith("0.5 * Sqrt[2]", "0.7071");
    assert_eval_startswith("Sqrt[2]/2.0",   "0.7071");
    /* MPFR-precision flag-style evaluation also goes numeric. */
    assert_eval_startswith("N[Sqrt[2]/2]",  "0.7071");
}

static void test_chained_canonicalization(void) {
    /* Multiple radical groups share num_prod via separate pull passes:
     * each group claims its own prime factors. */
    /* Sqrt[2]/2 * Sqrt[3]/3 = (1/Sqrt[2]) * (1/Sqrt[3]) = 1/Sqrt[6]
     * after the existing radical-fusion same-exp consolidation -- but
     * if that does not fire we still get an equivalent normal form.
     * Pin the value through a Simplify zero-test instead of a pinned
     * string. */
    assert_eval_eq("Simplify[(Sqrt[2]/2)(Sqrt[3]/3) Sqrt[6] - 1]", "0", 0);
    /* Nested:  (Sqrt[2]/2)^3 = 1/2^(3/2). */
    assert_eval_eq("(Sqrt[2]/2)^3",         "1/2^(3/2)", 0);
    /* Long sum of like radicals, all routed through canonicalization to
     * Power[2, -1/2] / Power[2, 1/2] forms.  Sqrt[8]/4 = 2 Sqrt[2]/4 =
     * 1/Sqrt[2]; total: Sqrt[2] + 2/Sqrt[2] = 2 Sqrt[2]. */
    assert_eval_eq("Sqrt[2] + 1/Sqrt[2] + Sqrt[8]/4",  "2 Sqrt[2]", 0);
    assert_eval_eq("Simplify[Sqrt[2] + 1/Sqrt[2] + Sqrt[8]/4 - 2 Sqrt[2]]",
                   "0", 0);
}

static void test_negative_coefficient(void) {
    /* Negative numerator is preserved -- canonicalization only moves
     * |b|^k factors, not signs. */
    assert_eval_eq("-Sqrt[2]/2",   "-1/Sqrt[2]", 0);
    assert_eval_eq("-2/Sqrt[2]",   "-Sqrt[2]", 0);
    assert_eval_eq("Sqrt[2]/(-2)", "-1/Sqrt[2]", 0);
    /* Negative coefficient out front (parser binds unary minus tightest;
     * (-1) * 3^(3/2) is the canonical sign-extracted form). */
    assert_eval_eq("-(3^(3/2))",                "-3 Sqrt[3]", 0);
    assert_eval_eq("(-1) 3^(3/2)",              "-3 Sqrt[3]", 0);
    assert_eval_eq("(-3)^(3/2) - (-3)^(3/2)",   "0", 0);
}

/* ------------------------------------------------------------------ */
/* Same-exponent product collapse (Times generalized fusion, k = +1)   */
/* ------------------------------------------------------------------ */

static void test_same_exponent_sqrt(void) {
    /* a^q * b^q -> (a*b)^q  for positive numeric bases */
    assert_eval_eq("Sqrt[2] Sqrt[3]",           "Sqrt[6]",      0);
    assert_eval_eq("Sqrt[2] Sqrt[5]",           "Sqrt[10]",     0);
    assert_eval_eq("Sqrt[3] Sqrt[5]",           "Sqrt[15]",     0);
    /* Three-way chain via repeated pairwise fusion. */
    assert_eval_eq("Sqrt[2] Sqrt[3] Sqrt[5]",   "Sqrt[30]",     0);
    assert_eval_eq("Sqrt[2] Sqrt[3] Sqrt[5] Sqrt[7]", "Sqrt[210]", 0);
}

static void test_same_exponent_cube_root(void) {
    /* Cube-root product collapse. */
    assert_eval_eq("2^(1/3) 3^(1/3)",           "6^(1/3)",      0);
    assert_eval_eq("2^(1/3) 5^(1/3)",           "10^(1/3)",     0);
    assert_eval_eq("3^(1/3) 5^(1/3) 7^(1/3)",   "105^(1/3)",    0);
    /* The product collapses to a perfect cube and Power's perfect-cube
     * reduction then fires -- 2^(1/3) * 4^(1/3) = 8^(1/3) = 2. */
    assert_eval_eq("2^(1/3) 4^(1/3)",           "2",            0);
    /* 4^(1/3) auto-unifies to 2^(2/3) before fusion can pair it with
     * 2^(1/3); the generalized-fusion k=+2 path then collapses
     * 2^(1/3) * 2^(2/3) -> 2^1 = 2. Same result either way. */
}

static void test_same_exponent_higher_q(void) {
    /* Quartic and quintic radicals collapse the same way. */
    assert_eval_eq("2^(1/4) 3^(1/4)",           "6^(1/4)",      0);
    assert_eval_eq("2^(1/5) 3^(1/5)",           "6^(1/5)",      0);
    /* k = +2 with q = 3: a^(1/3) * b^(2/3) -> (a*b^2)^(1/3) */
    assert_eval_eq("2^(1/3) 3^(2/3)",           "18^(1/3)",     0);
    assert_eval_eq("2^(2/3) 3^(1/3)",           "12^(1/3)",     0);
    /* Same exponent with rational coefficient. The (1/2) gets absorbed
     * into Sqrt[2]'s exponent before fusion, then the generalized
     * fusion sees Power[2, -1/2] * Power[3, 1/2] and produces the
     * combined-base form Sqrt[3/2] = Sqrt[6]/2. */
    assert_eval_eq("(1/2) Sqrt[2] Sqrt[3]",     "Sqrt[3/2]",    0);
}

/* ------------------------------------------------------------------ */
/* Generalized fusion for integer-ratio exponents (k = -2, +2, ...)    */
/* ------------------------------------------------------------------ */

static void test_generalized_fusion_negative_k(void) {
    /* k = -1: classic radical fusion (Sqrt[6]/Sqrt[2] -> Sqrt[3]) is
     * preserved; the regression that motivated the generalized version. */
    assert_eval_eq("Sqrt[6]/Sqrt[2]",           "Sqrt[3]",      0);
    assert_eval_eq("Sqrt[10]/Sqrt[5]",          "Sqrt[2]",      0);
    /* k = -2: a^(1/3) * b^(-2/3) -> (a/b^2)^(1/3). After perfect-power
     * unification 4^(-1/3) becomes 2^(-2/3), so this exact form arises
     * from 12^(1/3) * 4^(-1/3) too. */
    assert_eval_eq("12^(1/3) * 4^(-1/3)",       "3^(1/3)",      0);
    assert_eval_eq("12^(1/3) * 2^(-2/3)",       "3^(1/3)",      0);
    /* k = -3: a^(1/4) * b^(-3/4) -> (a/b^3)^(1/4). */
    assert_eval_eq("16^(1/4) * 2^(-3/4)",       "2^(1/4)",      0);
}

static void test_generalized_fusion_positive_k(void) {
    /* k = +2: a^(1/3) * b^(2/3) -> (a*b^2)^(1/3). */
    assert_eval_eq("2^(1/3) 8^(2/3)",           "4 2^(1/3)",    0);
    /* 2^(1/3) * (2^3)^(2/3) -> 2^(1/3) * 2^2 -> 4 * 2^(1/3) -- after
     * unification of 8^(2/3) into 2^2 the same outcome falls out
     * via simple integer multiplication. */
    /* k = +3: a^(1/4) * b^(3/4) -> (a*b^3)^(1/4). */
    assert_eval_eq("3^(1/4) 2^(3/4)",           "24^(1/4)",     0);
}

static void test_generalized_fusion_rational_ratio_skipped(void) {
    /* When the exponents are NOT in integer ratio, fusion must NOT
     * fire -- 2^(2/5) * 3^(3/5) shouldn't combine. (3/5 is not an
     * integer multiple of 2/5.) */
    assert_eval_eq("2^(2/5) 3^(3/5)",           "2^(2/5) 3^(3/5)", 0);
    /* Different denominators don't fuse either. */
    assert_eval_eq("2^(1/2) 3^(1/3)",           "Sqrt[2] 3^(1/3)", 0);
}

/* ------------------------------------------------------------------ */
/* Perfect-power base unification (Power.c)                            */
/* ------------------------------------------------------------------ */

static void test_perfect_power_base_unification(void) {
    /* The motivating case: 4^(2/3) -> 2^(4/3) -> 2 * 2^(1/3). */
    assert_eval_eq("4^(2/3)",                   "2 2^(1/3)",    0);
    /* 9 = 3^2, residual exponent 2/3 on a smaller base. */
    assert_eval_eq("9^(1/3)",                   "3^(2/3)",      0);
    assert_eval_eq("9^(2/3)",                   "3 3^(1/3)",    0);
    /* 8 = 2^3 is already minimal-base in some forms; perfect-cube
     * reduction collapses these to integers. */
    assert_eval_eq("8^(2/3)",                   "4",            0);
    assert_eval_eq("8^(5/3)",                   "32",           0);
    /* 16 = 2^4. Different exponents trigger different combinations. */
    assert_eval_eq("16^(1/4)",                  "2",            0);
    assert_eval_eq("16^(3/4)",                  "8",            0);
    assert_eval_eq("16^(2/3)",                  "4 2^(2/3)",    0);
    /* 32 = 2^5: lattice of cases. */
    assert_eval_eq("32^(1/5)",                  "2",            0);
    assert_eval_eq("32^(2/5)",                  "4",            0);
    /* 1024 = 2^10. */
    assert_eval_eq("1024^(1/2)",                "32",           0);
    assert_eval_eq("1024^(1/5)",                "4",            0);
    /* 27 = 3^3. */
    assert_eval_eq("27^(1/3)",                  "3",            0);
    assert_eval_eq("27^(2/3)",                  "9",            0);
    assert_eval_eq("27^(4/3)",                  "81",           0);
}

static void test_perfect_power_minimal_base(void) {
    /* The unification picks the SMALLEST base. 64 = 2^6, NOT 4^3 or 8^2,
     * so 64^(1/2) -> 8 (not 4 or 8 stays as 8) and 64^(1/3) -> 4. */
    assert_eval_eq("64^(1/2)",                  "8",            0);
    assert_eval_eq("64^(1/3)",                  "4",            0);
    assert_eval_eq("64^(1/6)",                  "2",            0);
    assert_eval_eq("64^(5/6)",                  "32",           0);
    /* 729 = 3^6 = 27^2 = 9^3; smallest base is 3. */
    assert_eval_eq("729^(1/3)",                 "9",            0);
    assert_eval_eq("729^(1/6)",                 "3",            0);
}

static void test_perfect_power_non_perfect_unchanged(void) {
    /* 6, 10, 12, 30, 100 (perfect square but caught by sqrt path),
     * 2, 3, 5 (primes). These should NOT have their bases rewritten by
     * the unification step. */
    assert_eval_eq("6^(1/3)",                   "6^(1/3)",      0);
    assert_eval_eq("6^(2/3)",                   "6^(2/3)",      0);
    assert_eval_eq("6^(2/5)",                   "6^(2/5)",      0);
    assert_eval_eq("10^(1/3)",                  "10^(1/3)",     0);
    assert_eval_eq("12^(1/3)",                  "12^(1/3)",     0);
    assert_eval_eq("30^(2/7)",                  "30^(2/7)",     0);
    /* 100 = 10^2 IS a perfect power; should unify to 10^(2*p/q). */
    assert_eval_eq("100^(1/3)",                 "10^(2/3)",     0);
    assert_eval_eq("100^(1/4)",                 "Sqrt[10]",     0);
}

static void test_perfect_power_negative_exp(void) {
    /* Negative-exponent unification: Power[4, -2/3] -> Power[2, -4/3].
     * Then integer-part extraction kicks in: -4/3 = -1 + (-1)/3 in the
     * sense p<0,a=-1 ... actually our existing block keeps a=0, b=-2
     * and produces a rational coefficient. The end form is consistent
     * with what Sqrt[2]/2 already produces. */
    assert_eval_eq("4^(-1/3)",                  "1/2^(2/3)",    0);
    assert_eval_eq("9^(-1/3)",                  "1/3^(2/3)",    0);
    assert_eval_eq("8^(-1/3)",                  "1/2",          0);
    assert_eval_eq("16^(-1/4)",                 "1/2",          0);
    assert_eval_eq("16^(-3/4)",                 "1/8",          0);
}

/* ------------------------------------------------------------------ */
/* Sign normalization for negative bases with odd q                    */
/* ------------------------------------------------------------------ */

static void test_sign_normalization_perfect_cube(void) {
    /* Mathematica form for negative perfect-cube bases:
     * (-8)^(p/q) -> coeff * (-1)^(b/q) when residue is 1. */
    assert_eval_eq("(-8)^(1/3)",                "2 (-1)^(1/3)", 0);
    assert_eval_eq("(-8)^(2/3)",                "4 (-1)^(2/3)", 0);
    assert_eval_eq("(-27)^(1/3)",               "3 (-1)^(1/3)", 0);
    assert_eval_eq("(-27)^(2/3)",               "9 (-1)^(2/3)", 0);
    assert_eval_eq("(-64)^(1/3)",               "4 (-1)^(1/3)", 0);
    assert_eval_eq("(-1000)^(1/3)",             "10 (-1)^(1/3)", 0);
}

static void test_sign_normalization_integer_part(void) {
    /* Higher-power exponents pick up a sign from the integer part:
     * (-8)^(5/3) = (-8)^1 * (-8)^(2/3) = -8 * 4 (-1)^(2/3) = -32 (-1)^(2/3).
     * a_int = 1 (odd), so coefficient is negated. */
    assert_eval_eq("(-8)^(5/3)",                "-32 (-1)^(2/3)", 0);
    /* a_int = 2 (even), no sign flip. */
    assert_eval_eq("(-8)^(7/3)",                "128 (-1)^(1/3)", 0);
    /* (-27)^(7/3) = (-27)^2 * (-27)^(1/3) = 729 * 3 * (-1)^(1/3). */
    assert_eval_eq("(-27)^(7/3)",               "2187 (-1)^(1/3)", 0);
    /* a_int = 1 (odd) for (-27)^(4/3). */
    assert_eval_eq("(-27)^(4/3)",               "-81 (-1)^(1/3)", 0);
}

static void test_sign_normalization_with_residue(void) {
    /* (-r) stays under the radical when r > 1. For (-72)^(1/3),
     * 72 = 2^3 * 9 so factor 2 out, residue (-9). */
    assert_eval_eq("(-72)^(1/3)",               "2 (-9)^(1/3)", 0);
    assert_eval_eq("(-72)^(2/3)",               "4 (-9)^(2/3)", 0);
    /* (-200)^(1/3): 200 = 2^3 * 25, factor 2, residue (-25). */
    assert_eval_eq("(-200)^(1/3)",              "2 (-25)^(1/3)", 0);
    /* (-24)^(1/3): 24 = 2^3 * 3, factor 2, residue (-3). */
    assert_eval_eq("(-24)^(1/3)",               "2 (-3)^(1/3)", 0);
}

static void test_sign_normalization_no_extraction(void) {
    /* When |n| has no perfect q-th-power factor, the negative-base form
     * stays as-is. (-2)^(1/3), (-3)^(1/3), (-6)^(1/3) — m = 1, no factor
     * to pull out, leave unevaluated. */
    assert_eval_eq("(-2)^(1/3)",                "(-2)^(1/3)", 0);
    assert_eval_eq("(-3)^(1/3)",                "(-3)^(1/3)", 0);
    assert_eval_eq("(-6)^(1/3)",                "(-6)^(1/3)", 0);
    assert_eval_eq("(-12)^(1/3)",               "(-12)^(1/3)", 0);
}

static void test_sign_normalization_even_q_unchanged(void) {
    /* Even q with negative base (other than q == 2) is intentionally
     * left unevaluated; the principal-branch result has no clean
     * canonical form yet. q == 2 still goes through the I-extraction
     * path, q == 4 / 6 / ... stays. */
    assert_eval_eq("(-4)^(1/2)",                "2*I",          0);
    assert_eval_eq("(-9)^(1/2)",                "3*I",          0);
    assert_eval_eq("(-2)^(1/4)",                "(-2)^(1/4)",   0);
    assert_eval_eq("(-2)^(1/6)",                "(-2)^(1/6)",   0);
    assert_eval_eq("(-16)^(1/4)",               "(-16)^(1/4)",  0);
}

/* ------------------------------------------------------------------ */
/* Integration: combinations of all three new canonical forms          */
/* ------------------------------------------------------------------ */

static void test_combined_canonical_forms(void) {
    /* Same-exponent product into a perfect-cube collapse. */
    assert_eval_eq("2^(1/3) 32^(1/3)",          "4",            0);
    /* Same-exponent product where unification then re-fuses. */
    assert_eval_eq("4^(1/3) 2^(1/3)",           "2",            0);
    /* Sign normalization combined with same-exponent collapse. */
    assert_eval_eq("(-8)^(1/3) 27^(1/3)",       "6 (-1)^(1/3)", 0);
    /* Perfect-power unification then radical fusion across coefficients.
     * 4^(2/3) -> 2 * 2^(1/3); divided by 2 cancels the 2: result 2^(1/3). */
    assert_eval_eq("4^(2/3)/2",                 "2^(1/3)",      0);
    /* Idempotence check: re-evaluating the canonical form is a no-op. */
    assert_eval_eq("FullForm[4^(2/3)]",
        "Times[2, Power[2, Rational[1, 3]]]", 0);
    assert_eval_eq("FullForm[(-8)^(1/3)]",
        "Times[2, Power[-1, Rational[1, 3]]]", 0);
    assert_eval_eq("FullForm[(-72)^(1/3)]",
        "Times[2, Power[-9, Rational[1, 3]]]", 0);
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    /* User-facing motivating examples. */
    TEST(test_user_sqrt_over_two);
    TEST(test_user_cube_root_over_two);
    TEST(test_user_radical_sum_collects);

    /* Power.c integer-part extraction. */
    TEST(test_power_extracts_unit_residue);
    TEST(test_power_extracts_with_perfect_factor);
    TEST(test_power_perfect_power_collapses);
    TEST(test_power_bigint_coefficient);
    TEST(test_power_negative_in_residue);
    TEST(test_power_irreducible_unchanged);
    TEST(test_power_negative_base_unchanged);

    /* Times.c radical canonicalization. */
    TEST(test_times_pulls_from_denominator);
    TEST(test_times_pulls_from_numerator_when_negative_exp);
    TEST(test_times_does_not_pull_when_positive_exp);
    TEST(test_times_unrelated_base_stays);
    TEST(test_times_rational_coefficient);

    /* Plus collection. */
    TEST(test_plus_collects_like_radicals);
    TEST(test_plus_distinct_radicals_remain);

    /* Compatibility with existing radical machinery. */
    TEST(test_radical_fusion_compat);
    TEST(test_sqrt_reduction_compat);
    TEST(test_radical_with_symbolic);

    /* Stability and round-trip. */
    TEST(test_canonical_form_idempotent);
    TEST(test_canonical_form_negative_exponent);
    TEST(test_round_trip_arithmetic);

    /* Stress / regression. */
    TEST(test_zero_and_one_exponents_unchanged);
    TEST(test_real_coefficient_skips);
    TEST(test_chained_canonicalization);
    TEST(test_negative_coefficient);

    /* Same-exponent product collapse (Times generalized fusion, k = +1). */
    TEST(test_same_exponent_sqrt);
    TEST(test_same_exponent_cube_root);
    TEST(test_same_exponent_higher_q);

    /* Generalized fusion for k = -2, +2, ... */
    TEST(test_generalized_fusion_negative_k);
    TEST(test_generalized_fusion_positive_k);
    TEST(test_generalized_fusion_rational_ratio_skipped);

    /* Perfect-power base unification. */
    TEST(test_perfect_power_base_unification);
    TEST(test_perfect_power_minimal_base);
    TEST(test_perfect_power_non_perfect_unchanged);
    TEST(test_perfect_power_negative_exp);

    /* Sign normalization for negative base, odd q. */
    TEST(test_sign_normalization_perfect_cube);
    TEST(test_sign_normalization_integer_part);
    TEST(test_sign_normalization_with_residue);
    TEST(test_sign_normalization_no_extraction);
    TEST(test_sign_normalization_even_q_unchanged);

    /* Combinations of the three new canonical forms. */
    TEST(test_combined_canonical_forms);

    printf("All radical_canonical tests passed.\n");
    return 0;
}
