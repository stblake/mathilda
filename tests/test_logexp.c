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

void test_logexp_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Log[E]", "1"},
        {"Log[1]", "0"},
        {"Log[0]", "-Infinity"},
        {"Log[Infinity]", "Infinity"},
        {"Log[2, 8]", "3"},
        {"Log[b, b]", "1"},
        {"Exp[0]", "1"},
        {"Exp[-Infinity]", "0"},
        {"Exp[Infinity]", "Infinity"},
        {"Log[10, 10]", "1"},
        {"Log[E, x]", "Log[x]"},
        {"Exp[0.0]", "1.0"},
        {"Log[-1.0]", "0.0 + 3.14159*I"},
        {"Log[-1]", "(I) Pi"},
        {"Log[-5]", "Log[5] + (I) Pi"},
        {"Exp[Log[x]]", "x"},
        {"Exp[b Log[a]]", "a^b"},
        {"Exp[2 Log[x]]", "x^2"},
        {"Exp[3 Log[2]]", "8"},
        {"3^Log[3, x]", "x"},
        {"b^Log[b, a]", "a"},
        {"Power[base, b Log[base, a]]", "a^b"},
        {"10^(3 Log[10, x])", "x^3"},
        {"Log[E^4]", "4"},
        {"Log[E^(1/3)]", "1/3"},
        {"Log[E^(-2)]", "-2"},
        {"Log[E^x]", "Log[E^x]"},
        {"Log[2, 2^(1/3)]", "1/3"},
        {"Log[E, E^4]", "4"},
        {"Table[Exp[I * n * Pi / 2], {n, 0, 4}]", "{1, I, -1, -I, 1}"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_logexp_forward);
    
    return 0;
}
