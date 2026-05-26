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

void test_piecewise() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        // Integer inputs
        {"Round[5]", "5"},
        {"Floor[5]", "5"},
        {"Ceiling[5]", "5"},
        {"IntegerPart[5]", "5"},
        {"FractionalPart[5]", "0"},
        
        // Real inputs
        {"Floor[2.4]", "2"},
        {"Ceiling[2.4]", "3"},
        {"IntegerPart[2.4]", "2"},
        {"FractionalPart[2.4]", "0.4"},
        {"Round[2.4]", "2"},
        {"Round[2.5]", "2"},
        {"Round[3.5]", "4"},
        {"Floor[-2.4]", "-3"},
        {"Ceiling[-2.4]", "-2"},
        {"IntegerPart[-2.4]", "-2"},
        {"FractionalPart[-2.4]", "-0.4"},
        {"Round[-2.5]", "-2"},
        {"Round[-3.5]", "-4"},
        
        // Arbitrary-precision (MPFR) Real inputs — must reduce to exact integers
        // for Floor/Ceiling (the sibling MPFR coverage follows in later commits).
        {"Floor[N[3.7, 50]]", "3"},
        {"Floor[N[-3.7, 50]]", "-4"},
        {"Floor[N[5, 50]]", "5"},
        {"Ceiling[N[3.2, 50]]", "4"},
        {"Ceiling[N[-3.2, 50]]", "-3"},
        {"Ceiling[N[5, 50]]", "5"},

        // Rational inputs
        {"Floor[5/2]", "2"},
        {"Ceiling[5/2]", "3"},
        {"IntegerPart[5/2]", "2"},
        {"FractionalPart[5/2]", "Rational[1, 2]"},
        {"Round[5/2]", "2"},
        {"Round[7/2]", "4"},
        {"Floor[-5/2]", "-3"},
        {"Ceiling[-5/2]", "-2"},
        {"IntegerPart[-5/2]", "-2"},
        {"FractionalPart[-5/2]", "Rational[-1, 2]"},
        {"Round[-5/2]", "-2"},
        {"Round[-7/2]", "-4"},
        
        // Complex inputs
        {"Floor[Complex[2.4, -2.4]]", "Complex[2, -3]"},
        {"Round[Complex[2.5, -3.5]]", "Complex[2, -4]"},
        {"IntegerPart[Complex[5/2, -5/2]]", "Complex[2, -2]"},
        {"FractionalPart[Complex[5/2, -5/2]]", "Complex[Rational[1, 2], Rational[-1, 2]]"},
        {"FractionalPart[Complex[2, 3]]", "0"}, // Real 0, Im 0 -> integer 0
        
        // Sign extraction (symbolic)
        {"Floor[-x]", "Times[-1, Ceiling[x]]"},
        {"Ceiling[-x]", "Times[-1, Floor[x]]"},
        {"Round[-x]", "Times[-1, Round[x]]"},
        {"Floor[-2 x]", "Times[-1, Ceiling[Times[2, x]]]"},
        {"Ceiling[Times[-3, x, y]]", "Times[-1, Floor[Times[3, x, y]]]"},
        {"Floor[-x] + Ceiling[x]", "0"},

        // Idempotency / composition (symbolic)
        {"Floor[Floor[x]]", "Floor[x]"},
        {"Ceiling[Ceiling[x]]", "Ceiling[x]"},
        {"Round[Round[x]]", "Round[x]"},
        {"Floor[Ceiling[x]]", "Ceiling[x]"},
        {"Ceiling[Floor[x]]", "Floor[x]"},
        {"Floor[Round[x]]", "Round[x]"},
        {"Ceiling[Round[x]]", "Round[x]"},
        {"Round[Floor[x]]", "Floor[x]"},
        {"Round[Ceiling[x]]", "Ceiling[x]"},
        {"Ceiling[Floor[Ceiling[x]]]", "Ceiling[x]"},

        // IntegerPart / FractionalPart do NOT participate in these rules
        {"IntegerPart[-x]", "IntegerPart[Times[-1, x]]"},
        {"FractionalPart[Floor[x]]", "FractionalPart[Floor[x]]"},

        // Two argument forms
        {"Floor[7, 3]", "6"},
        {"Ceiling[7, 3]", "9"},
        {"Round[7, 3]", "6"},
        {"Round[2.4, 0.5]", "2.5"},
        
        // Infinity inputs
        {"Floor[Infinity]", "Infinity"},
        {"Ceiling[Infinity]", "Infinity"},
        {"Round[Infinity]", "Infinity"},
        {"IntegerPart[Infinity]", "Infinity"},
        {"FractionalPart[Infinity]", "0"},
        {"Floor[-Infinity]", "Times[-1, Infinity]"},
        {"Ceiling[-Infinity]", "Times[-1, Infinity]"},
        {"Round[-Infinity]", "Times[-1, Infinity]"},
        {"IntegerPart[-Infinity]", "Times[-1, Infinity]"},
        {"FractionalPart[-Infinity]", "0"},
        
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_piecewise);
    
    return 0;
}
