#include "expand.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

void test_expand() {
    run_test("Expand[(x+3)(x+2)]", "Plus[6, Times[5, x], Power[x, 2]]");
    run_test("Expand[(x+y)^2 (x-y)^2]", "Plus[Power[x, 4], Times[-2, Times[Power[x, 2], Power[y, 2]]], Power[y, 4]]");
    run_test("Expand[(x+y)^5 - (x^5+y^5)]", "Plus[Times[5, Times[Power[x, 4], y]], Times[10, Times[Power[x, 3], Power[y, 2]]], Times[10, Times[Power[x, 2], Power[y, 3]]], Times[5, Times[x, Power[y, 4]]]]");
    run_test("Expand[(f[Subscript[x, 1]] + f[Subscript[x, 2]])^3]", "Plus[Power[f[Subscript[x, 1]], 3], Power[f[Subscript[x, 2]], 3], Times[3, Times[Power[f[Subscript[x, 1]], 2], f[Subscript[x, 2]]]], Times[3, Times[Power[f[Subscript[x, 2]], 2], f[Subscript[x, 1]]]]]");
    run_test("Expand[(Sin[x] + Cos[x])^2]", "Plus[Power[Cos[x], 2], Power[Sin[x], 2], Times[2, Times[Cos[x], Sin[x]]]]");
    run_test("Expand[Sqrt[(1+x)^2]]", "Power[Power[Plus[1, x], 2], Rational[1, 2]]");
    run_test("Expand[1 < (x+y)^2 < 2]", "Inequality[1, Less, Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]], Less, 2]");
    
    // Pattern-based Expand
    run_test("Expand[(x+1)^2 + (y+1)^2, x]", "Plus[1, Times[2, x], Power[x, 2], Power[Plus[1, y], 2]]");
    run_test("Expand[(x+1)^2 + (y+1)^2, y]", "Plus[1, Power[Plus[1, x], 2], Times[2, y], Power[y, 2]]");
    run_test("Expand[(x+1)^2 + (y+1)^2, _Symbol]", "Plus[2, Times[2, x], Power[x, 2], Times[2, y], Power[y, 2]]");
    run_test("Expand[(x^2+1)^2 + (y+1)^2, x^2]", "Plus[1, Times[2, Power[x, 2]], Power[x, 4], Power[Plus[1, y], 2]]");
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_expand);
    
    printf("All expand tests passed!\n");
    return 0;
}
