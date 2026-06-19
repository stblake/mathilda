#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("Rescale test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Two-argument numeric form: (x - min)/(max - min) ---- */
void test_rescale_two_arg_numeric() {
    check("Rescale[2.5, {-10, 10}]", "0.625");
    check("Rescale[12.5, {-10, 10}]", "1.125");      /* may exceed [0,1] */
    check("Rescale[-3/2, {-2, 2}]", "1/8");          /* exact rational */
    check("Rescale[3, {0, 6}]", "1/2");
    check("Rescale[0, {0, 10}]", "0");
    check("Rescale[10, {0, 10}]", "1");
}

/* ---- Endpoints rescale to exactly 0 and 1 ---- */
void test_rescale_endpoints() {
    check("Rescale[a, {a, b}]", "0");
    check("Rescale[a, {b, a}]", "1");
    check("Rescale[b, {a, b}]", "1");
}

/* ---- Real numbers track input precision ---- */
void test_rescale_real() {
    check("Rescale[Pi, {0, 2.5}]", "1.25664");
    /* Precision of the output tracks the precision of the input. */
    check("Rescale[8, {-9, 7.111111111111111111111}]",
          "1.055172413793103448276");
}

/* ---- Complex inputs ---- */
void test_rescale_complex() {
    check("Rescale[1 + 2 I, {0, 5}]", "1/5 + 2/5*I");
    check("Rescale[1 + 2 I, {0, 1 + I}]", "3/2 + 1/2*I");
}

/* ---- High-precision evaluation via N ---- */
void test_rescale_high_precision() {
    check("N[Rescale[1/11, {1/7, 5}], 50]",
          "-0.0106951871657754010695187165775401069518716577540107");
    check("N[Rescale[1/8, {1/17, 1/11}], 100]", "2.0625");
}

/* ---- Three-argument form: y0 + (y1 - y0)(x - min)/(max - min) ---- */
void test_rescale_three_arg_numeric() {
    check("Rescale[3, {-9, 7}, {11, 28}]", "95/4");
    check("Rescale[a, {a, b}, {c, d}]", "c");   /* lower endpoint -> y0 */
    check("Rescale[b, {a, b}, {c, d}]", "d");   /* upper endpoint -> y1 */
}

/* ---- Degenerate target range collapses to the constant ---- */
void test_rescale_degenerate_target() {
    check("Rescale[x, {a, b}, {c, c}]", "c");
    check("Rescale[5, {0, 10}, {7, 7}]", "7");
}

/* ---- Symbolic forms (kept factored; verified equivalent below) ---- */
void test_rescale_symbolic() {
    check("Rescale[0, {a, b}]", "-a/(-a + b)");
    check("Rescale[0, {a, b}, {c, d}]", "c - (a (-c + d))/(-a + b)");
    check("Rescale[x, {a, b}, {c, d}]", "c + ((-c + d) (-a + x))/(-a + b)");
}

/* ---- Symbolic forms are mathematically equal to the closed forms ---- */
void test_rescale_symbolic_equivalence() {
    check("Simplify[Rescale[x, {a, b}] - (x - a)/(b - a)]", "0");
    check("Simplify[Rescale[0, {a, b}] - (-(a/(-a + b)))]", "0");
    check("Simplify[Rescale[x, {a, b}, {c, d}] - "
          "(c + (d - c) (x - a)/(b - a))]", "0");
    check("Simplify[Rescale[0, {a, b}, {c, d}] - "
          "(-((-b c + a d)/(-a + b)))]", "0");
}

/* ---- One-argument form rescales a list over its own Min..Max ---- */
void test_rescale_one_arg_list() {
    check("Rescale[{-2, 0, 2}]", "{0, 1/2, 1}");
    check("Rescale[{-.7, .5, 1.2, 5.6, 1.8}]",
          "{0.0, 0.190476, 0.301587, 1.0, 0.396825}");
    /* Constant list: Min == Max so every element divides 0/0. */
    check("Rescale[{1, 1, 1}]", "{Indeterminate, Indeterminate, Indeterminate}");
}

/* ---- Threading over a list first argument ---- */
void test_rescale_threading() {
    check("Rescale[{-2, 0, 2}, {-5, 5}]", "{3/10, 1/2, 7/10}");
    check("Rescale[{-2, 0, 2}, {-5, 5}, {-1, 1}]", "{-2/5, 0, 2/5}");
    check("Rescale[{1, 2, 3, 4, 5, 6}, {0, a}]",
          "{1/a, 2/a, 3/a, 4/a, 5/a, 6/a}");
    /* Threads recursively through nested lists. */
    check("Rescale[{{1, 2}, {3, 4}}, {0, 4}]", "{{1/4, 1/2}, {3/4, 1}}");
}

/* ---- Wrong argument counts leave the call unevaluated (after argb) ---- */
void test_rescale_arg_errors() {
    check("Rescale[]", "Rescale[]");
    check("Rescale[1, 3, 4, 5, 6]", "Rescale[1, 3, 4, 5, 6]");
}

/* ---- Malformed range arguments are left unevaluated ---- */
void test_rescale_bad_range() {
    check("Rescale[x, {a, b, c}]", "Rescale[x, {a, b, c}]");
    check("Rescale[x, a]", "Rescale[x, a]");
    check("Rescale[x, {a, b}, {c, d, e}]", "Rescale[x, {a, b}, {c, d, e}]");
}

/* ---- Attributes: NumericFunction and Protected ---- */
void test_rescale_attributes() {
    check("MemberQ[Attributes[Rescale], NumericFunction]", "True");
    check("MemberQ[Attributes[Rescale], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_rescale_two_arg_numeric);
    TEST(test_rescale_endpoints);
    TEST(test_rescale_real);
    TEST(test_rescale_complex);
    TEST(test_rescale_high_precision);
    TEST(test_rescale_three_arg_numeric);
    TEST(test_rescale_degenerate_target);
    TEST(test_rescale_symbolic);
    TEST(test_rescale_symbolic_equivalence);
    TEST(test_rescale_one_arg_list);
    TEST(test_rescale_threading);
    TEST(test_rescale_arg_errors);
    TEST(test_rescale_bad_range);
    TEST(test_rescale_attributes);

    printf("All Rescale tests passed!\n");
    return 0;
}
