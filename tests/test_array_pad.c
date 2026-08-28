#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "ndarray.h"
#include "pack.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Compare an evaluated result (delisted so a packed buffer reads as a List) in
 * OutputForm against `expected`. */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = test_delist(evaluate(e));
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "ArrayPad %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Constant padding, amount forms ---------- */

static void test_ap_scalar(void) {
    run("ArrayPad[{1,2,3},1]", "{0, 1, 2, 3, 0}");
}
static void test_ap_matrix_scalar(void) {
    run("ArrayPad[{{1,2},{3,4}},2]",
        "{{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 1, 2, 0, 0}, "
        "{0, 0, 3, 4, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}");
}
static void test_ap_pair(void) {
    run("ArrayPad[{1,2,3},{2,4}]", "{0, 0, 1, 2, 3, 0, 0, 0, 0}");
    run("ArrayPad[{1,2,3},{0,2}]", "{1, 2, 3, 0, 0}");
}
static void test_ap_constant_value(void) {
    run("ArrayPad[{1,2,3},2,x]", "{x, x, 1, 2, 3, x, x}");
}

/* ---------- Per-level amounts ---------- */

static void test_ap_perlevel_symmetric(void) {
    run("ArrayPad[{{1,2},{3,4}},{{1},{5}}]",
        "{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "
        "{0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0}, "
        "{0, 0, 0, 0, 0, 3, 4, 0, 0, 0, 0, 0}, "
        "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}");
}
static void test_ap_perlevel_rows(void) {
    run("ArrayPad[{{1,2},{3,4}}, {{1},{0}}]", "{{0, 0}, {1, 2}, {3, 4}, {0, 0}}");
}
static void test_ap_perlevel_cols(void) {
    run("ArrayPad[{{1,2},{3,4}}, {{0},{1}}]", "{{0, 1, 2, 0}, {0, 3, 4, 0}}");
}
static void test_ap_perlevel_mixed(void) {
    run("ArrayPad[{{1,2},{3,4}},{{-1,2},{1,3}}]",
        "{{0, 3, 4, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}");
}
static void test_ap_first_level_only(void) {
    run("ArrayPad[{{1,2},{3,4}},{{1}}]", "{0, {1, 2}, {3, 4}, 0}");
}

/* ---------- Negative amounts (removal) ---------- */

static void test_ap_negative(void) {
    run("ArrayPad[Range[10],-2]", "{3, 4, 5, 6, 7, 8}");
}

/* ---------- Named schemes ---------- */

static void test_ap_fixed(void) {
    run("ArrayPad[{1,2,3},2,\"Fixed\"]", "{1, 1, 1, 2, 3, 3, 3}");
}
static void test_ap_periodic(void) {
    run("ArrayPad[{a,b,c},4,\"Periodic\"]", "{c, a, b, c, a, b, c, a, b, c, a}");
}
static void test_ap_reversed(void) {
    run("ArrayPad[{a,b,c},4,\"Reversed\"]", "{c, c, b, a, a, b, c, c, b, a, a}");
}
static void test_ap_reversed_negation(void) {
    run("ArrayPad[{a,b,c},4,\"ReversedNegation\"]",
        "{c, -c, -b, -a, a, b, c, -c, -b, -a, a}");
}
static void test_ap_reflected(void) {
    run("ArrayPad[{a,b,c},4,\"Reflected\"]", "{a, b, c, b, a, b, c, b, a, b, c}");
}
static void test_ap_reflected_differences(void) {
    run("ArrayPad[{a,b,c},3,\"ReflectedDifferences\"]",
        "{2 a + b - 2 c, 2 a - c, 2 a - b, a, b, c, -b + 2 c, -a + 2 c, -2 a + b + 2 c}");
}
static void test_ap_reversed_differences(void) {
    run("ArrayPad[{a,b,c},4,\"ReversedDifferences\"]",
        "{2 a - c, 2 a - c, 2 a - b, a, a, b, c, c, -b + 2 c, -a + 2 c, -a + 2 c}");
}

/* ---------- Extrapolated ---------- */

static void test_ap_extrapolated_linear_symbolic(void) {
    run("ArrayPad[{a,b,c},3,\"Extrapolated\"]",
        "{4 a - 3 b, 3 a - 2 b, 2 a - b, a, b, c, -b + 2 c, -2 b + 3 c, -3 b + 4 c}");
}
static void test_ap_extrapolated_order2(void) {
    run("ArrayPad[{1,5,7,8},{0,3},\"Extrapolated\",InterpolationOrder->2]",
        "{1, 5, 7, 8, 8, 7, 5}");
}
static void test_ap_extrapolated_default(void) {
    run("ArrayPad[{1,5,7,8},{0,3},\"Extrapolated\"]", "{1, 5, 7, 8, 9, 10, 11}");
}
static void test_ap_extrapolated_infinity(void) {
    run("ArrayPad[{1,5,7,8},{0,3},\"Extrapolated\",InterpolationOrder->Infinity]",
        "{1, 5, 7, 8, 9, 11, 15}");
}
static void test_ap_extrapolated_orders(void) {
    run("Table[ArrayPad[{1,4,9},4,\"Extrapolated\", InterpolationOrder->o],{o,0,2}]",
        "{{1, 1, 1, 1, 1, 4, 9, 9, 9, 9, 9}, "
        "{-11, -8, -5, -2, 1, 4, 9, 14, 19, 24, 29}, "
        "{9, 4, 1, 0, 1, 4, 9, 16, 25, 36, 49}}");
}

/* ---------- Empty dimensions ---------- */

static void test_ap_empty_constant(void) {
    run("ArrayPad[{{},{}},1,x]", "{{x, x}, {x, x}, {x, x}, {x, x}}");
}
static void test_ap_empty_fixed(void) {
    run("ArrayPad[{{},{}},1,\"Fixed\"]", "{{}, {}, {}, {}}");
}

/* ---------- mindimsize ---------- */

static void test_ap_mindimsize(void) {
    /* value-difference scheme, axis of length < 2 -> unevaluated. (expr_to_string
     * keeps the quotes on the scheme string; OutputForm/Print drops them, as in
     * Wolfram's transcript.) */
    run("ArrayPad[{{1,2,4}},{1,2},\"ReflectedDifferences\"]",
        "ArrayPad[{{1, 2, 4}}, {1, 2}, \"ReflectedDifferences\"]");
    /* explicit per-dim amounts avoid it */
    run("ArrayPad[{{1,2,4}},{{0,0},{1,2}},\"ReflectedDifferences\"]",
        "{{0, 1, 2, 4, 6, 7}}");
}

/* ---------- Packed input ---------- */

static void test_ap_packed_int(void) {
    run("ArrayPad[Range[5],2]", "{0, 0, 1, 2, 3, 4, 5, 0, 0}");
    run("ArrayPad[Range[3],2,7]", "{7, 7, 1, 2, 3, 7, 7}");
    run("ArrayPad[Range[10],-2]", "{3, 4, 5, 6, 7, 8}");
}
static void test_ap_packed_real(void) {
    run("ArrayPad[Range[1.,5.],{0,2},9.]", "{1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 9.0}");
    run("ArrayPad[Range[1.,5.],2]", "{0, 0, 1.0, 2.0, 3.0, 4.0, 5.0, 0, 0}");
}

/* ---------- Attributes / documentation / options ---------- */

static void test_ap_protected(void) {
    run("MemberQ[Attributes[ArrayPad], Protected]", "True");
}
static void test_ap_options(void) {
    run("Options[ArrayPad]", "{InterpolationOrder -> 1}");
}
static void test_ap_docstring(void) {
    SymbolDef* def = symtab_get_def("ArrayPad");
    ASSERT_MSG(def && def->docstring && def->docstring[0],
               "ArrayPad should have a non-empty docstring");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_ap_arity(void) {
    run("ArrayPad[]", "ArrayPad[]");
    run("ArrayPad[{1,2,3}]", "ArrayPad[{1, 2, 3}]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_ap_scalar);
    TEST(test_ap_matrix_scalar);
    TEST(test_ap_pair);
    TEST(test_ap_constant_value);
    TEST(test_ap_perlevel_symmetric);
    TEST(test_ap_perlevel_rows);
    TEST(test_ap_perlevel_cols);
    TEST(test_ap_perlevel_mixed);
    TEST(test_ap_first_level_only);
    TEST(test_ap_negative);
    TEST(test_ap_fixed);
    TEST(test_ap_periodic);
    TEST(test_ap_reversed);
    TEST(test_ap_reversed_negation);
    TEST(test_ap_reflected);
    TEST(test_ap_reflected_differences);
    TEST(test_ap_reversed_differences);
    TEST(test_ap_extrapolated_linear_symbolic);
    TEST(test_ap_extrapolated_order2);
    TEST(test_ap_extrapolated_default);
    TEST(test_ap_extrapolated_infinity);
    TEST(test_ap_extrapolated_orders);
    TEST(test_ap_empty_constant);
    TEST(test_ap_empty_fixed);
    TEST(test_ap_mindimsize);
    TEST(test_ap_packed_int);
    TEST(test_ap_packed_real);
    TEST(test_ap_protected);
    TEST(test_ap_options);
    TEST(test_ap_docstring);
    TEST(test_ap_arity);

    printf("All ArrayPad tests passed!\n");
    return 0;
}
