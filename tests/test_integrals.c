#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include <stdio.h>
#include <string.h>

void test_integrals() {
    symtab_init();
    core_init();
    
    Expr* load_cmd = parse_expression("Get[\"src/internal/CRCMathTablesIntegrals.m\"]");
    Expr* res = evaluate(load_cmd);
    
    if (res->type == EXPR_SYMBOL && strcmp(res->data.symbol, "$Failed") == 0) {
        expr_free(res);
        expr_free(load_cmd);
        load_cmd = parse_expression("Get[\"../../src/internal/CRCMathTablesIntegrals.m\"]");
        res = evaluate(load_cmd);
    }
    
    expr_free(res);
    expr_free(load_cmd);
    
    printf("Formula 59\\n");
    assert_eval_eq("Integrate[(1 + 2x)^3/(3 + 4x)^2, x]", "1/4 x^2 + 1/(96 + 128 x) + 3/32 Log[3 + 4 x]", 0);
    
    printf("Formula 60\n");
    assert_eval_eq("Integrate[1/(3 + 4x^2), x]", "(1/2 ArcTan[(2 x)/Sqrt[3]])/Sqrt[3]", 0);
    
    printf("Formula 61\n");
    assert_eval_eq("Integrate[1/(3 - 4x^2), x]", "-1/4 Log[-1/2 Sqrt[3] + x]/Sqrt[3] + 1/4 Log[1/2 Sqrt[3] + x]/Sqrt[3]", 0);
    
    printf("Formula 62\n");
    assert_eval_eq("Integrate[1/(9 + 16x^2), x]", "1/12 ArcTan[4/3 x]", 0);
    
    printf("Formula 63\n");
    assert_eval_eq("Integrate[x/(3 + 4x^2), x]", "1/8 Log[3 + 4 x^2]", 0);
    
    printf("Formula 64\n");
    assert_eval_eq("Integrate[x^2/(3 + 4x^2), x]", "1/4 x - 1/8 Sqrt[3] ArcTan[(2 x)/Sqrt[3]]", 0);
    printf("Done!\n");
}

int main() {
    TEST(test_integrals);
    return 0;
}
