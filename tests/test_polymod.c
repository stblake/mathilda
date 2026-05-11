#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

void test_polymod_basic(void) {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"PolynomialMod[3x^2+2x+1,2]", "1 + x^2"},
        {"PolynomialMod[3x^2+2x+1,x^2+1]", "-2 + 2 x"},
        {"PolynomialMod[3x^3+21x^2-7x+55,9]", "1 + 2 x + 3 x^2 + 3 x^3"},
        {"PolynomialMod[35x^3+21x^2*y^2-17x*y^3+55z-123,19]", "10 + 16 x^3 + 2 x^2 y^2 + 2 x y^3 + 17 z"},
        {"PolynomialMod[3x^3+21x^2-7x+55,2x^2-7]", "257/2 + 7/2 x"},
        {"PolynomialMod[3x^3+21x^2-7x+55,x^2+I]", "55 - 21*I + (-7 - 3*I) x"},
        {"PolynomialMod[3x^3+21x^2-7x+55,{2x^2-7,9}]", "8 + x^2 + x^3"},
        {"PolynomialMod[3x^3+21x^2*y^2-7x*y^3+55,{2x^2-7,x*y-3, 9}]", "1 + 7 x + x^3 + 4 y^2"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("FAIL: %s\n  expected: %s\n  got:      %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(false);
        } else {
            printf("PASS: %s -> %s\n", tests[i].input, tests[i].expected);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

int main() {
    symtab_init();
    core_init();
    
    printf("Running polymod tests...\n");
    TEST(test_polymod_basic);
    printf("All polymod tests passed!\n");
    return 0;
}