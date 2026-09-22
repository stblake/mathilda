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

void test_purefunc_basic() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Function[u, 3+u][x]", "Plus[3, x]"},
        {"Function[3+#][x]", "Plus[3, x]"},
        {"(3+#)&[x]", "Plus[3, x]"},
        {"Function[{u, v}, u^2 + v^4][x, y]", "Plus[Power[x, 2], Power[y, 4]]"},
        {"(#1^2 + #2^4)&[x, y]", "Plus[Power[x, 2], Power[y, 4]]"},
        
        // Unbound slots remain unevaluated
        {"(3+#)&[]", "Plus[3, Slot[1]]"},
        
        // Sequence slots
        {"f[##]&[a, b, c, d]", "f[a, b, c, d]"},
        {"f[X, ##, Y, ##]&[a, b, c, d]", "f[X, a, b, c, d, Y, a, b, c, d]"},
        {"f[##2]&[a, b, c, d]", "f[b, c, d]"}, // SlotSequence[2]
        
        // Mapping a pure function over a list
        {"g[#, #^2]& /@ {x, y, z}", "List[g[x, Power[x, 2]], g[y, Power[y, 2]], g[z, Power[z, 2]]]"},
        
        // Array from a pure function
        {"Array[1+#^2&, 3]", "List[2, 5, 10]"},
        
        // Listable pure functions
        {"Function[{u}, g[u], Listable][{a, b, c}]", "List[g[a], g[b], g[c]]"},
        {"Function[{u}, g[u]][{a, b, c}]", "g[List[a, b, c]]"},
        
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        if (strcmp(s, cases[i].expected) != 0) {
            printf("Forward %s: expected %s, got %s\n", cases[i].input, cases[i].expected, s);
        }
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

/* Function[params, body, attrs]: attributes are respected during the
 * argument-evaluation phase, and lexical substitution means held references
 * to a parameter (e.g. Unevaluated[x]) see the substituted expression. */
void test_purefunc_attributes() {
    /* HoldAll holds the argument, so x is bound to Plus[1,1,1] unevaluated.
     * Length[Unevaluated[Plus[1,1,1]]] gives 3. */
    assert_eval_eq("Function[{x}, Length[Unevaluated[x]], {HoldAll}][1+1+1]", "3", 0);
    /* Single attribute (not wrapped in List) is also accepted. */
    assert_eval_eq("Function[{x}, Length[Unevaluated[x]], HoldAll][1+1+1]", "3", 0);

    /* No Hold attribute: the argument is evaluated first (1+1+1 -> 3) before
     * substitution, so Length[Unevaluated[3]] = Length[3] = 0. */
    assert_eval_eq("Function[{x}, Length[Unevaluated[x]]][1+1+1]", "0", 0);

    /* Default Function has no Hold attribute, matching Mathematica:
     * (Hold[#]&)[1+2] evaluates the arg to 3 before substitution, so the
     * result is Hold[3] (not Hold[1+2]). */
    assert_eval_eq("(Hold[#]&)[1+2]", "Hold[3]", 0);

    /* HoldFirst holds arg 1 but evaluates arg 2. */
    assert_eval_eq("Function[{a, b}, {Length[Unevaluated[a]], Length[Unevaluated[b]]}, HoldFirst][1+2+3, 4+5+6]",
                   "{3, 0}", 0);

    /* Listable threads the pure function over a list argument. */
    assert_eval_eq("Function[{u}, g[u], Listable][{a, b, c}]",
                   "{g[a], g[b], g[c]}", 0);
}

/* Lexical shadowing: the inner Function's parameter shadows the outer. */
void test_purefunc_shadowing() {
    assert_eval_eq("Function[{x}, x + Function[{x}, x][10]][1]", "11", 0);
    assert_eval_eq("Function[{x}, x + (# &)[99]][7]", "106", 0);
}

/* Function is a lexical closure (Wolfram-Language semantics): an outer named
 * parameter is visible inside a nested pure function unless that nested
 * Function rebinds the SAME name. Slots (#) are still claimed by the innermost
 * enclosing Function of any form. */
void test_purefunc_closure() {
    /* Outer named parameter reaches inside a nested slot function -- the case
     * that used to leave the parameter free (Part[pl,1] survived). */
    assert_eval_eq("Function[pl, Map[(pl[[1]] + #) &, {10, 20}]][{5, 9}]",
                   "{15, 25}", 0);
    /* The partial-application idiom: Function[x, (x + #)&][a] -> (a + #)&. */
    assert_eval_eq("Function[x, (x + #) &][10][20]", "30", 0);
    /* Named parameter reaches inside a nested NAMED function it does not
     * shadow (previously returned Function[x, x^n] with n free). */
    assert_eval_eq("Function[n, Function[x, x^n]][2][3]", "9", 0);
    assert_eval_eq("Function[x, Function[y, x + y]][3][4]", "7", 0);
    /* Three levels deep. */
    assert_eval_eq("Function[a, Function[b, Function[c, a + b + c]]][1][2][3]",
                   "6", 0);
    /* Shadowing: an inner Function that rebinds the name blocks substitution. */
    assert_eval_eq("Function[x, Function[x, x]][1][2]", "2", 0);
    /* Partial shadowing across a multi-parameter frame: x rebound, y not. */
    assert_eval_eq("Function[{x, y}, Function[{x}, x + y]][10, 20][5]", "25", 0);
    /* Slot shielding is preserved: the inner slot function keeps its own #. */
    assert_eval_eq("Function[Function[# + 1]][5][10]", "11", 0);
    /* A10: an outer SLOT-form function's # must flow into a nested NAMED
     * function's body (previously left #1 free -> "2 + #1"). */
    assert_eval_eq("Function[Function[y, # + y]][1][2]", "3", 0);
    assert_eval_eq("Function[Function[{y}, #1 * y]][3][4]", "12", 0);
    /* Nested Map over an association-row parameter (the shape that surfaced the
     * bug in the ParallelMixed integrator's `splittable` gate). */
    assert_eval_eq("Map[Function[r, Map[r[[1]] + # &, {1, 2}]], {{10}, {20}}]",
                   "{{11, 12}, {21, 22}}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_purefunc_basic);
    TEST(test_purefunc_attributes);
    TEST(test_purefunc_shadowing);
    TEST(test_purefunc_closure);

    return 0;
}
