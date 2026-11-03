/* test_mathar_regressions.c
 *
 * Regression tests for the wrong answers and memory errors the Mathar integral
 * corpus (tests/Mathar*IntegralsCorpus.m, arXiv:1207.5845) exposed.  Each test
 * names the corpus entry that found it and pins the root-cause fix, not just
 * the integral:
 *
 *   series_arctan_sign      Series/Limit of ArcTan[k x] at infinity took +Pi/2
 *                           for every non-literal k (src/calculus/series.c);
 *                           a negative radical k made Newton-Leibniz telescope
 *                           to a wrong definite value (M6968).
 *   squarefree_bigint       FactorSquareFree enrolled bigint coefficients as
 *                           polynomial variables and lost the square
 *                           (src/poly/poly.c is_number); the rational
 *                           integrator then returned 0 (M5695).
 *   squarefree_verified     FactorSquareFree returned a product unequal to its
 *                           input when its gcds were not true divisors
 *                           (src/poly/facpoly_squarefree.inc); through
 *                           Simplify, a wrong residue closed form (singular
 *                           corpus case 49).
 *   pmt_foreign_kernel      ParallelMixedTower treated a non-elementary kernel
 *                           of an x-dependent argument as a constant and leaked
 *                           its tower variable into the antiderivative
 *                           (src/calculus/integrate.c, M13537).
 *   int64_min_rational      Rational[INT64_MIN, d] trapped on n / -1 in
 *                           make_rational (src/arithmetic.c, M7801).
 *   tc_unwind_block_values  a TimeConstrained unwind freed Block values still in
 *                           use (src/modular.c; eight corpus crashes).
 *   residue_region          the symbolic residue families classified poles at
 *                           one sample point although the radical root
 *                           expressions swap across the region
 *                           (src/calculus/integrate_residue.c, M6968).
 *
 * Values that Mathilda leaves in an equivalent unsimplified form are pinned
 * numerically via Chop[N[value - reference]] == 0.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

static void test_series_arctan_sign(void) {
    /* Numeric negative radical leading coefficient. */
    check_eq("Limit[2 ArcTan[x/(1 - Sqrt[3])], x -> Infinity]", "-Pi");
    check_eq("Limit[Sqrt[2] ArcTan[x/(1 - Sqrt[3])], x -> Infinity]", "-Pi/Sqrt[2]");
    /* The constant is (Pi/2) Sign[(1 - Sqrt[3]) w]: the reference's
     * -(Pi/2) w Sqrt[w^-2] generically, -Pi/2 for w > 0 (v0.6162). */
    check_eq("Series[ArcTan[1/((1 - Sqrt[3]) w)], {w, 0, 1}] === SeriesData[w, 0, {-Pi Sqrt[w^(-2)] w/2, -1 + Sqrt[3]}, 0, 2, 1]", "True");
    check_eq("Series[ArcTan[1/((1 - Sqrt[3]) w)], {w, 0, 1}, Assumptions -> w > 0]",
             "-1/2 Pi + (-1 + Sqrt[3]) w + O[w]^2");
    /* Positive coefficient: unchanged. */
    check_eq("Limit[2 ArcTan[x/(Sqrt[3] - 1)], x -> Infinity]", "Pi");
    /* Symbolic coefficient: the sign factor, never an assumed sign. */
    /* (w -> 0 in the constant presumes the side: w > 0, stated.) */
    check_eq("Chop[N[(Normal[Series[ArcTan[c/w], {w, 0, 1}, Assumptions -> w > 0]] /. w -> 0) /. c -> -2]]", "-1.5708");
    check_eq("Chop[N[(Normal[Series[ArcTan[c/w], {w, 0, 1}, Assumptions -> w > 0]] /. w -> 0) /. c -> 3]]", "1.5708");
    /* The definite integral that exposed it (M6968 instance). */
    check_eq("Chop[N[Integrate[1/(1 + 2 x^2 + x^4/2), {x, 0, Infinity}] - Sqrt[1 - 1/Sqrt[2]] Pi/2]]", "0");
    check_eq("Chop[N[Integrate[1/(1 + 5/3 x^2 + 3/8 x^4), {x, 0, Infinity}] - 0.9237722717237228]]", "0");
}

static void test_squarefree_bigint(void) {
    check_eq("FactorSquareFree[Expand[3^45 (1 + t)^2]]", "2954312706550833698643 (1 + t)^2");
    check_eq("FactorSquareFree[Expand[687194767360000 (16 + 625 t^2)^2]]",
             "687194767360000 (16 + 625 t^2)^2");
    /* The rational integrals that returned 0 (M5695 instance). */
    check_eq("Integrate[x/(625/4096 + x^4), x]", "32/25 ArcTan[64/25 x^2]");
    check_eq("Integrate[x/(81/4096 + x^4), x]", "32/9 ArcTan[64/9 x^2]");
}

static void test_squarefree_verified(void) {
    /* FactorSquareFree must reproduce its input whatever it factors. */
    check_eq("Expand[FactorSquareFree[x^2 - 2 x Arg[I a] + Arg[I a]^2 + Pi x]"
             " - (x^2 - 2 x Arg[I a] + Arg[I a]^2 + Pi x)]", "0");
    /* Singular-definite corpus case 49: Pi^3/(8 a) + Pi Log[a]^2/(2 a). */
    check_eq("Chop[N[(Integrate[Log[x]^2/(a^2 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]"
             " - (Pi^3/(8 a) + Pi Log[a]^2/(2 a))) /. a -> 2]]", "0");
    check_eq("Chop[N[(Integrate[Log[x]^2/(a^2 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]"
             " - (Pi^3/(8 a) + Pi Log[a]^2/(2 a))) /. a -> 1/3]]", "0");
}

static void test_residue_region(void) {
    /* 1/(1 + 2 a x^2 + b x^4): the poles swap across a^2 = b inside a, b > 0.
     * Whatever closed form comes back must be right on BOTH sides. */
    const char* r = "rr = Integrate[1/(1 + 2 a x^2 + b x^4), {x, 0, Infinity}, "
                    "Assumptions -> a > 0 && b > 0]";
    Expr* p = parse_expression(r);
    Expr* v = evaluate(p);
    expr_free(p);
    bool unevaluated = v && v->type == EXPR_FUNCTION &&
        v->data.function.head->type == EXPR_SYMBOL &&
        strcmp(v->data.function.head->data.symbol.name, "Integrate") == 0;
    expr_free(v);
    if (!unevaluated) {
        check_eq("Chop[N[(rr /. {a -> 5/6, b -> 3/8}) - 0.9237722717237228]]", "0");
        check_eq("Chop[N[(rr /. {a -> 1/5, b -> 3}) - 0.7990898378557233]]", "0");
    }
}

static void test_pmt_foreign_kernel(void) {
    /* Not elementary: must stay unevaluated, never a constant-folded g[r..]. */
    check_eq("Integrate[g[Cos[t]] Sin[n (z - t)], t]", "Integrate[g[Cos[t]] Sin[n (-t + z)], t]");
    check_eq("FreeQ[Integrate[LegendreP[1 + n, -m, Cos[t]] Sin[n (z - t)], t], "
             "s_Symbol /; StringMatchQ[SymbolName[s], \"r\" ~~ DigitCharacter ..]]", "True");
    /* Elementary integrands still reach the stage. */
    check_eq("D[Integrate[Sqrt[x + Sqrt[x]], x], x] - Sqrt[x + Sqrt[x]] // Simplify", "0");
}

static void test_int64_min_rational(void) {
    check_eq("Rational[-9223372036854775808, -2]", "4611686018427387904");
    check_eq("Rational[-9223372036854775808, 3]", "-9223372036854775808/3");
    check_eq("Rational[6, -9223372036854775808]", "-3/4611686018427387904");
}

/* A TimeConstrained unwind inside ParallelMixedSpecial freed a Block-local
 * memo Association still in use (segfault, M6534's instance reproduced it in
 * ~20 s).  Any result is fine; the process must survive. */
static void test_tc_unwind_block_values(void) {
    check_eq("Head[Integrate[Log[Sin[5/4] + Sin[(3*x)/8]], {x, 0, 5/8}]] =!= Null", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_series_arctan_sign);
    TEST(test_squarefree_bigint);
    TEST(test_squarefree_verified);
    TEST(test_residue_region);
    TEST(test_pmt_foreign_kernel);
    TEST(test_int64_min_rational);
    TEST(test_tc_unwind_block_values);

    printf("All Mathar regression tests passed!\n");
    return 0;
}
