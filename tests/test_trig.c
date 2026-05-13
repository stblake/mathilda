
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_trig_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[0]", "0"},
        {"Sin[Pi/2]", "1"},
        {"Sin[Pi/6]", "Rational[1, 2]"},
        {"Sin[Pi/4]", "Power[2, Rational[-1, 2]]"},
        {"Sin[Pi/3]", "Times[Rational[1, 2], Power[3, Rational[1, 2]]]"},
        {"Cos[0]", "1"},
        {"Cos[Pi/2]", "0"},
        {"Cos[Pi/3]", "Rational[1, 2]"},
        {"Tan[Pi/4]", "1"},
        {"Cot[Pi/4]", "1"},
        {"Sec[0]", "1"},
        {"Csc[Pi/2]", "1"},
        // Exact values for d=12
        {"Sin[Pi/12]", "Times[Rational[1, 4], Plus[Power[6, Rational[1, 2]], Times[-1, Power[2, Rational[1, 2]]]]]"},
        {"Cos[Pi/12]", "Times[Rational[1, 4], Plus[Power[2, Rational[1, 2]], Power[6, Rational[1, 2]]]]"},
        // Exact values for d=10, 5
        {"Sin[Pi/10]", "Times[Rational[1, 4], Plus[-1, Power[5, Rational[1, 2]]]]"},
        {"Cos[Pi/5]", "Times[Rational[1, 4], Plus[1, Power[5, Rational[1, 2]]]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_inverse() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"ArcSin[0]", "0"},
        {"ArcSin[1]", "Times[Rational[1, 2], Pi]"},
        {"ArcSin[1/2]", "Times[Rational[1, 6], Pi]"},
        {"ArcCos[0]", "Times[Rational[1, 2], Pi]"},
        {"ArcCos[1]", "0"},
        {"ArcTan[0]", "0"},
        {"ArcTan[1]", "Times[Rational[1, 4], Pi]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_forward_of_inverse() {
    /* f[f_inv[x]] = x for each direct trig / inverse-trig pair. Uses
     * InputForm-style strings since the simplification strips the inverse
     * and leaves the inner argument untouched. */
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[ArcSin[x]]", "x"},
        {"Cos[ArcCos[y]]", "y"},
        {"Tan[ArcTan[a]]", "a"},
        {"Cot[ArcCot[b]]", "b"},
        {"Sec[ArcSec[c]]", "c"},
        {"Csc[ArcCsc[d]]", "d"},
        /* Composite arguments survive as-is */
        {"Sin[ArcSin[x^2 + 1]]", "1 + x^2"},
        {"Cos[ArcCos[2 y]]", "2 y"},
        /* ArcTan[x, y] two-arg form must NOT be stripped: our rule guards
         * on arg_count==1, so Tan[ArcTan[x, y]] stays unevaluated here
         * rather than collapsing via the (wrong) single-argument rule. */
        {"Tan[ArcTan[3, 4]]", "Tan[ArcTan[3, 4]]"},
        /* Opposite direction is NOT folded: ArcSin[Sin[x]] stays put */
        {"ArcSin[Sin[x]]", "ArcSin[Sin[x]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0,
                   "Forward-of-inverse %s: expected %s, got %s",
                   cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

int main() {
    symtab_init();
    core_init();

    TEST(test_trig_forward);
    TEST(test_trig_inverse);
    TEST(test_trig_forward_of_inverse);

    return 0;
}
