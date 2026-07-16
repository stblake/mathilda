#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Assert that Apart[input] prints exactly as `expected`. */
static void assert_apart_form(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        free(res_str); expr_free(e); expr_free(res);
        ASSERT(false);
    }
    printf("PASS: %s -> %s\n", input, expected);
    free(res_str); expr_free(e); expr_free(res);
}

/* Correctness by recombination: Together[Apart[f] - f] must be exactly 0.
 * This is immune to output-form canonicalisation (unlike a string compare):
 * it directly proves the partial-fraction decomposition is mathematically
 * equal to the original rational function.  `var` may be NULL (1-arg Apart
 * with automatic variable selection) or an explicit variable string. */
static void assert_apart_recombines(const char* f, const char* var) {
    char buf[4096];
    if (var)
        snprintf(buf, sizeof(buf),
                 "Together[(Apart[%s, %s]) - (%s)]", f, var, f);
    else
        snprintf(buf, sizeof(buf),
                 "Together[(Apart[%s]) - (%s)]", f, f);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, "0") != 0) {
        printf("FAIL recombination: Apart[%s%s%s]\n"
               "  Together[Apart[f] - f] = %s  (expected 0)\n",
               f, var ? ", " : "", var ? var : "", res_str);
        free(res_str); expr_free(e); expr_free(res);
        ASSERT(false);
    }
    printf("PASS recombination: Apart[%s%s%s]\n",
           f, var ? ", " : "", var ? var : "");
    free(res_str); expr_free(e); expr_free(res);
}

/* Baseline pure-rational-over-Q forms (all route through the FLINT fast path).
 * These output forms are stable — exact string compare.  Multivariate cases,
 * whose printed form drifts with Factor/Together canonicalisation, are checked
 * by recombination instead (test_apart_multivariate_decline). */
void test_apart_basic(void) {
    struct { const char* input; const char* expected; } tests[] = {
        {"Apart[1/((1+x)(5+x))]", "-1/4/(5 + x) + 1/4/(1 + x)"},
        {"Apart[(x^5-2)/((1+x+x^2)(2+x)(1-x))]", "2 - x - 34/9/(2 + x) + 1/9/(-1 + x) + (-1 - 1/3 x)/(1 + x + x^2)"},
        {"Apart[1/(1-x^3)]", "-1/3/(-1 + x) + (2/3 + 1/3 x)/(1 + x + x^2)"},
        {"Apart[16x/((1+x)^2*(5+x))]", "-5/(5 + x) - 4/(1 + x)^2 + 5/(1 + x)"},
    };
    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++)
        assert_apart_form(tests[i].input, tests[i].expected);
}

/* Exact output-form locks for the FLINT-native pure-rational-over-Q path:
 * distinct/repeated linear, irreducible + repeated irreducible quadratic,
 * improper fractions, non-monic and rational-coefficient factors.  Output must
 * match the classical RowReduce Apart byte-for-byte. */
void test_apart_flint_forms(void) {
    struct { const char* input; const char* expected; } tests[] = {
        /* irreducible quadratic: (a + b x)/(x^2 + 1) numerator preserved */
        {"Apart[(2x+1)/((x^2+1)(x-3)), x]",
         "7/10/(-3 + x) + (-1/10 - 7/10 x)/(1 + x^2)"},
        /* repeated irreducible quadratic + a linear factor */
        {"Apart[(x^3+1)/((x^2+x+1)^2 (x+5)), x]",
         "-124/441/(5 + x) + (8/21 - 2/21 x)/(1 + x + x^2)^2 + (-55/441 + 124/441 x)/(1 + x + x^2)"},
        /* improper fraction: polynomial (quotient) part present */
        {"Apart[(x^5+1)/((x+1)(x+2)(x+3)), x]",
         "25 - 6 x + x^2 - 121/(3 + x) + 31/(2 + x)"},
        /* repeated linear */
        {"Apart[1/((x-1)^2 (x+2)), x]",
         "-1/9/(-1 + x) + 1/9/(2 + x) + 1/3/(-1 + x)^2"},
        /* rational leading content in the denominator */
        {"Apart[(x+1)/(2 (x-1)(x-3)), x]",
         "1/(-3 + x) - 1/2/(-1 + x)"},
        /* non-monic linear factors */
        {"Apart[1/((2x+1)(3x-2)), x]",
         "-2/7/(1 + 2 x) + 3/7/(-2 + 3 x)"},
        /* two repeated linear factors, higher multiplicity */
        {"Apart[(x^2+1)/((x-1)^3 (x+2)^2), x]",
         "-5/27/(2 + x)^2 - 1/27/(2 + x) + 1/27/(-1 + x) + 2/27/(-1 + x)^2 + 2/9/(-1 + x)^3"},
        /* two linear + one irreducible quadratic (x^4 - 1) */
        {"Apart[1/(x^4-1), x]",
         "-1/4/(1 + x) + 1/4/(-1 + x) - 1/2/(1 + x^2)"},
        /* repeated irreducible quadratic + linear, non-trivial numerators */
        {"Apart[(5x^2-3)/((x^2+2)^2 (x-1)), x]",
         "2/9/(-1 + x) + (13/3 + 13/3 x)/(2 + x^2)^2 + (-2/9 - 2/9 x)/(2 + x^2)"},
    };
    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++)
        assert_apart_form(tests[i].input, tests[i].expected);
}

/* Broad correctness sweep of the FLINT path via recombination — every shape a
 * proper/improper rational function over Q can take. */
void test_apart_recombination(void) {
    const char* cases[] = {
        /* distinct linear */
        "1/((x+1)(x+2)(x+3))",
        "(x^2+x+1)/((x-1)(x-2)(x-3)(x-4))",
        "(3x-7)/((x+5)(x-11)(2x+1))",
        /* repeated linear, various multiplicities */
        "1/((x-1)^2 (x+2))",
        "(x^2+1)/((x-1)^3 (x+2)^2)",
        "1/((x+1)^5)",
        "(x^3-2x+5)/((x-2)^4 (x+3)^3)",
        /* irreducible quadratics (ArcTan-form numerators) */
        "(2x+1)/((x^2+1)(x-3))",
        "1/(x^3+1)",
        "1/(x^4+1)",
        "(x+2)/((x^2+x+1)(x^2+2))",
        /* repeated irreducible quadratics */
        "(x^3+1)/((x^2+x+1)^2 (x+5))",
        "(5x^2-3)/((x^2+2)^2 (x-1))",
        "1/((x^2+1)^3)",
        /* improper fractions (nonzero polynomial part) */
        "(x^5+1)/((x+1)(x+2)(x+3))",
        "(x^7)/((x^2+1)(x-2))",
        "(x^6-3x^2+2)/((x-1)^2 (x+4))",
        /* rational / large integer coefficients, non-monic factors */
        "(x+1)/(2 (x-1)(x-3))",
        "1/((2x+1)(3x-2))",
        "(7x-3)/((5x+2)^2 (x-9))",
        "(x^2+100000)/((x-123)(x+456))",
        /* everything at once: linear + repeated linear + quadratic + repeated */
        "(x^2-2)/((x+1)^3 (x^2+2)^2 (x-4))",
        /* higher degree — exercises the packed FLINT loop */
        "(x-99)/((x+2)^8 (x-3)^6)",
        "1/((x-1)^10 (x+1)^8)",
        "(x^3+x+1)/((x^2+1)^6 (x-2)^4)",
    };
    for (int i = 0; i < (int)(sizeof(cases) / sizeof(cases[0])); i++)
        assert_apart_recombines(cases[i], "x");
}

/* The multivariate / algebraic decline path (coefficients carrying another
 * symbol) must fall back to the classical elimination and still recombine. */
void test_apart_multivariate_decline(void) {
    assert_apart_recombines("1/((x+a)(x+b)(x+c))", "x");
    assert_apart_recombines("(a x + 1)/((x+1)(x+a))", "x");
    assert_apart_recombines("(x+y)/((x+1)(y+1)(x-y))", "x");
    assert_apart_recombines("(x+y)/((x+1)(y+1)(x-y))", "y");
    assert_apart_recombines("1/((x-k)^2 (x+k))", "x");
    /* radical coefficient: also declines the fmpq_poly conversion */
    assert_apart_recombines("1/((x - Sqrt[2])(x + Sqrt[3]))", "x");
}

/* Structural edge cases. */
void test_apart_edge_cases(void) {
    /* constant (degree-0) denominator in var -> pure Expand, no fractions */
    assert_apart_form("Apart[(x^2+2x+1)/4, x]", "1/4 + 1/2 x + 1/4 x^2");
    /* already a polynomial: unchanged (Expanded) */
    assert_apart_form("Apart[x^2 + 3 x + 2, x]", "2 + 3 x + x^2");
    /* single simple pole */
    assert_apart_form("Apart[1/(x-5), x]", "1/(-5 + x)");
    /* numerator degree exactly one less than denominator (proper, no quotient) */
    assert_apart_recombines("(2x+3)/((x+1)(x+4))", "x");
    /* negative leading coefficients throughout */
    assert_apart_recombines("(-x-1)/((-x+2)(-x-5))", "x");
    /* zero numerator -> 0 */
    assert_apart_form("Apart[0/((x+1)(x+2)), x]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running parfrac tests...\n");
    TEST(test_apart_basic);
    TEST(test_apart_flint_forms);
    TEST(test_apart_recombination);
    TEST(test_apart_multivariate_decline);
    TEST(test_apart_edge_cases);
    printf("All parfrac tests passed!\n");
    return 0;
}
