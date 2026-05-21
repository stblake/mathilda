/*
 * test_deriv.c -- Correctness tests for the native D, Dt, Derivative.
 *
 * Exercises the full dispatch surface of src/deriv.c: the partial-
 * derivative fast paths for Plus/Times/Power/List, every elementary
 * unary function derivative, the chain rule on unknown single and
 * multi-argument functions, Derivative[n][f][g] reduction, higher-order
 * derivatives via {x, n}, mixed derivatives via D[f, x, y, ...], and
 * the full complement of Dt (total) identities.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Compare evaluate(input) against expected, in FullForm so intermediate
 * orderings don't trip us up. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Derivative mismatch for %s:\n    expected: %s\n    got:      %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* Symbolic-equality check: reduce (lhs) - (rhs) via Expand[Together[...]]
 * and require the final form to be 0. Handy when the canonical form
 * from the printer differs from our hard-coded spelling. */
static void check_zero(const char* expr) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Expand[Together[%s]]", expr);
    Expr* e = parse_expression(buf);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", buf);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, "0") == 0,
               "Expected zero for %s, got %s", expr, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* --- Base cases & structural rules ------------------------------------ */
static void test_base_cases(void) {
    check("D[5, x]", "0");
    check("D[Pi, x]", "0");
    check("D[x, x]", "1");
    check("D[y, x]", "0");
    check("D[x, y]", "0");
}

/* --- Sum, product and power rules ------------------------------------- */
static void test_algebraic_rules(void) {
    check("D[a x + b y, x]", "a");
    check("D[x^2, x]", "Times[2, x]");
    check("D[x^3, x]", "Times[3, Power[x, 2]]");
    /* product rule on three factors */
    check_zero("D[x y z, x] - y z");
    check_zero("D[x y z, y] - x z");
    /* f^g with both varying -- x^x */
    check_zero("D[x^x, x] - x^x (1 + Log[x])");
}

/* --- Elementary unary chain rules ------------------------------------- */
static void test_elementary_chain(void) {
    check_zero("D[Sin[x], x] - Cos[x]");
    check_zero("D[Cos[x], x] + Sin[x]");
    check_zero("D[Tan[x], x] - Sec[x]^2");
    check_zero("D[Cot[x], x] + Csc[x]^2");
    check_zero("D[Sec[x], x] - Sec[x] Tan[x]");
    check_zero("D[Csc[x], x] + Csc[x] Cot[x]");

    check_zero("D[ArcSin[x], x] - 1/Sqrt[1 - x^2]");
    check_zero("D[ArcCos[x], x] + 1/Sqrt[1 - x^2]");
    check_zero("D[ArcTan[x], x] - 1/(1 + x^2)");
    check_zero("D[ArcCot[x], x] + 1/(1 + x^2)");

    check_zero("D[Sinh[x], x] - Cosh[x]");
    check_zero("D[Cosh[x], x] - Sinh[x]");
    check_zero("D[Tanh[x], x] - Sech[x]^2");

    check_zero("D[Exp[x], x] - Exp[x]");
    check_zero("D[Log[x], x] - 1/x");
    check_zero("D[Sqrt[x], x] - 1/(2 Sqrt[x])");

    /* chain rule with a non-trivial inner argument */
    check_zero("D[Sin[x^2], x] - 2 x Cos[x^2]");
    check_zero("D[Exp[a x], x] - a Exp[a x]");
    check_zero("D[Log[b, x], x] - 1/(x Log[b])");
}

/* --- Unknown-function chain rule (Derivative introduction) ------------ */
static void test_chain_unknown(void) {
    check("D[f[x], x]", "Derivative[1][f][x]");
    check_zero("D[f[g[x]], x] - Derivative[1][f][g[x]] Derivative[1][g][x]");
    check("D[f[x, y], x]", "Derivative[1, 0][f][x, y]");
    check("D[f[x, y], y]", "Derivative[0, 1][f][x, y]");
}

/* --- Chain on existing Derivative expressions ------------------------- */
static void test_chain_on_derivative(void) {
    check("D[Derivative[2][f][x], x]", "Derivative[3][f][x]");
    check_zero("D[Derivative[1][f][x^2], x] - 2 x Derivative[2][f][x^2]");
}

/* --- Higher-order derivatives and mixed partials ---------------------- */
static void test_higher_order(void) {
    check("D[x^3, {x, 2}]", "Times[6, x]");
    check("D[x^3, {x, 3}]", "6");
    check("D[x^3, {x, 0}]", "Power[x, 3]");
    check_zero("D[Sin[a x], {x, 2}] + a^2 Sin[a x]");
    /* mixed partial should commute */
    check_zero("D[x^2 y^3, x, y] - D[x^2 y^3, y, x]");
    check_zero("D[x^2 y^3, x, y] - 6 x y^2");
}

/* --- Threading across lists ------------------------------------------- */
static void test_list_threading(void) {
    check("D[{x, x^2, Sin[x]}, x]", "List[1, Times[2, x], Cos[x]]");
}

/* --- Total derivative (Dt) -------------------------------------------- */
static void test_dt(void) {
    check("Dt[5]", "0");
    check("Dt[Pi]", "0");
    check("Dt[E]", "0");
    check("Dt[I]", "0");
    check("Dt[x]", "Dt[x]");
    /* Known elementary functions still fire. */
    check_zero("Dt[Sin[x]] - Cos[x] Dt[x]");
    check_zero("Dt[y^2] - 2 y Dt[y]");
    /* Dt with a variable argument reduces to D. */
    check("Dt[x^2, x]", "Times[2, x]");
    check("Dt[y^2, x]", "0");
    check("Dt[x^2, {x, 2}]", "2");
}

/* --- Distribution over Equal -------------------------------------------- */
static void test_equal_distribution(void) {
    /* MMA: D[a == b, x] -> 0 == 0 -> True (both sides x-free).
     * (a == a evaluates to True *before* reaching D, so the canonical
     * test uses two distinct symbols.) */
    check("D[a == b, x]", "True");
    /* Plus inside Equal -- the classic implicit-circle setup. */
    check("D[x^2 + y^2 == 1, x]", "Equal[Times[2, x], 0]");
    /* Both sides depend on x. */
    check("D[Sin[x] == Cos[x], x]",
          "Equal[Cos[x], Times[-1, Sin[x]]]");
    /* Nested Equal distributes through both layers. (Mathilda parses
     * a == b == c left-associatively as Equal[Equal[a, b], c].) */
    check("D[Sin[x] == Cos[x] == Tan[x], x]",
          "Equal[Equal[Cos[x], Times[-1, Sin[x]]], Power[Sec[x], 2]]");
}

/* --- NonConstants option ------------------------------------------------ */
static void test_nonconstants(void) {
    /* The implicit-derivative base case: a symbol declared NonConstant
     * produces a literally-unevaluated D[sym, var, NonConstants -> {...}]. */
    check("D[y, x, NonConstants -> {y}]",
          "D[y, x, Rule[NonConstants, List[y]]]");
    /* Single-symbol shorthand canonicalises to a List. */
    check("D[y, x, NonConstants -> y]",
          "D[y, x, Rule[NonConstants, List[y]]]");
    /* Constant short-circuit must NOT fire on symbols that mention a
     * declared-non-constant symbol. */
    check_zero("D[2 y, x, NonConstants -> {y}] - 2 D[y, x, NonConstants -> {y}]");
    /* Symbols not in the list still vanish. */
    check("D[z, x, NonConstants -> {y}]", "0");
    /* x itself differentiates to 1, regardless of NonConstants. */
    check("D[x, x, NonConstants -> {y}]", "1");
    /* The full implicit-circle worked example from the MMA docs. */
    check_zero("D[x^2 + y^2 == 1, x, NonConstants -> y] - "
               "(2 x + 2 y D[y, x, NonConstants -> {y}] == 0)");
    /* Multiple non-constant symbols. */
    check_zero("D[x + y + z, x, NonConstants -> {y, z}] - "
               "(1 + D[y, x, NonConstants -> {y, z}] + "
               "D[z, x, NonConstants -> {y, z}])");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_base_cases);
    TEST(test_algebraic_rules);
    TEST(test_elementary_chain);
    TEST(test_chain_unknown);
    TEST(test_chain_on_derivative);
    TEST(test_higher_order);
    TEST(test_list_threading);
    TEST(test_dt);
    TEST(test_equal_distribution);
    TEST(test_nonconstants);

    printf("All derivative tests passed.\n");
    return 0;
}
