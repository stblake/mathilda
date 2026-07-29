/* Integer radical extraction: Power[n, p/q] for exact integer n.
 *
 * `Sqrt[n]` splits n into a rational coefficient and a q-th-power-free
 * residue, which needs n's prime factorisation. That factorisation used to
 * be trial division all the way to sqrt(n) — up to 3e9 iterations for an
 * int64 — so a single 18-19 digit non-square radicand took tens of seconds:
 *
 *     Sqrt[3141592653589793238]     24 s        (FactorInteger[same]: 3.6 ms)
 *     Sqrt[1000000000000000003]     20 s
 *
 * The residue past a trial-division limit now goes to the same general
 * factoriser FactorInteger uses (Pollard rho / SQUFOF / ECM), which answers
 * in milliseconds. This file pins both halves of that: every extracted form
 * is unchanged, and the cost no longer scales with sqrt(n).
 *
 * Run binary directly: ./power_integer_radical_tests
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_to_string(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT_MSG(parsed != NULL, "parse failed: %s", input);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    ASSERT_MSG(r != NULL, "evaluate returned NULL: %s", input);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

static void check(const char* input, const char* expected) {
    char* s = eval_to_string(input);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "%s\n  expected: %s\n  actual:   %s", input, expected, s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Extraction is unchanged
 * ---------------------------------------------------------------------- */

/* The small radicands the trial-division prefix still handles directly.
 * These are the regression guard for the fast path. */
static void test_small_radicands(void) {
    check("Sqrt[4]",        "2");
    check("Sqrt[8]",        "2 Sqrt[2]");
    check("Sqrt[12]",       "2 Sqrt[3]");
    check("Sqrt[18]",       "3 Sqrt[2]");
    check("Sqrt[72]",       "6 Sqrt[2]");
    check("Sqrt[2^62]",     "2147483648");
    check("Sqrt[-12]",      "(2*I) Sqrt[3]");
    check("Sqrt[1/12]",     "1/2/Sqrt[3]");
    check("12^(1/3)",       "2^(2/3) 3^(1/3)");
    check("54^(1/3)",       "3 2^(1/3)");
    check("(16/81)^(1/4)",  "2/3");
}

/* Radicands past the trial-division limit, where the general factoriser now
 * does the work. 3141592653589793238 = 2*3*11*10513*311743*14523877 — no
 * repeated prime, so nothing comes out and the radical stays whole. */
static void test_large_radicands(void) {
    check("Sqrt[3141592653589793238]", "Sqrt[3141592653589793238]");
    check("Sqrt[1000000000000000003]", "Sqrt[1000000000000000003]");
    check("Sqrt[999999999999999989]",  "Sqrt[999999999999999989]");
    check("Sqrt[2^61 - 1]",            "Sqrt[2305843009213693951]");
    check("Sqrt[10^18 + 1]",           "Sqrt[1000000000000000001]");

    /* A square factor well above the trial-division limit must still come
     * out — this is what a bare cutoff with no fallback would have broken. */
    check("Sqrt[1000003^2 * 1000033]", "1000003 Sqrt[1000033]");
    check("Sqrt[9 * 1000000000000000003]", "3 Sqrt[1000000000000000003]");
    check("(8 * 1000000000000000003)^(1/3)", "2 1000000000000000003^(1/3)");

    /* Perfect powers still collapse completely. */
    check("Sqrt[4611686018427387904]", "2147483648");
    check("Sqrt[1000000000000000000000000]", "1000000000000");
}

/* Squaring the extracted form must give the radicand back, whatever the
 * split was — an independent check on the coefficient/residue pair. */
static void test_round_trip(void) {
    static const char* radicands[] = {
        "3141592653589793238",
        "1000000000000000003",
        "1000003^2 * 1000033",
        "9 * 1000000000000000003",
        "999999999999999989",
        "72",
        "12",
    };
    for (size_t i = 0; i < sizeof(radicands) / sizeof(radicands[0]); ++i) {
        char in[128], expect[128];
        snprintf(in, sizeof(in), "Simplify[Sqrt[%s]^2]", radicands[i]);
        snprintf(expect, sizeof(expect), "%s", radicands[i]);
        char* got = eval_to_string(in);
        char* want = eval_to_string(expect);
        ASSERT_MSG(strcmp(got, want) == 0,
                   "%s\n  expected: %s\n  actual:   %s", in, want, got);
        free(got);
        free(want);
    }
}

/* ------------------------------------------------------------------------
 *  Cost no longer scales with sqrt(n)
 * ---------------------------------------------------------------------- */

/* Twelve 19-digit radicands. Before the fix each cost ~31 s (measured), so
 * this loop took ~6 minutes; it now runs in well under a tenth of a second.
 * The bound is absolute rather than normalized because the margin is three
 * orders of magnitude either way — no plausible machine sits between "a few
 * hundred milliseconds" and "several minutes". */
static void test_large_radicands_are_not_trial_divided(void) {
    const double budget_seconds = 5.0;
    clock_t t0 = clock();
    for (int k = 1; k <= 12; ++k) {
        char in[128];
        snprintf(in, sizeof(in), "Sqrt[%lld]",
                 (long long)(1000000000000000003LL + 2 * k));
        char* s = eval_to_string(in);
        free(s);
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < budget_seconds,
               "12 large-integer Sqrt[] calls took %.2f s (budget %.1f s) — "
               "the radicand factorisation has regressed to trial division",
               secs, budget_seconds);
    printf("  12 large radicands in %.3f s\n", secs);
}

/* Same check one exponent up, so a regression confined to the k-th-power
 * path (factor_out_kth_power) is caught too. */
static void test_large_cube_roots_are_not_trial_divided(void) {
    const double budget_seconds = 5.0;
    clock_t t0 = clock();
    for (int k = 1; k <= 12; ++k) {
        char in[128];
        snprintf(in, sizeof(in), "(%lld)^(1/3)",
                 (long long)(1000000000000000003LL + 2 * k));
        char* s = eval_to_string(in);
        free(s);
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < budget_seconds,
               "12 large-integer cube roots took %.2f s (budget %.1f s)",
               secs, budget_seconds);
    printf("  12 large cube roots in %.3f s\n", secs);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_small_radicands);
    TEST(test_large_radicands);
    TEST(test_round_trip);
    TEST(test_large_radicands_are_not_trial_divided);
    TEST(test_large_cube_roots_are_not_trial_divided);

    printf("All power integer-radical tests passed!\n");
    return 0;
}
