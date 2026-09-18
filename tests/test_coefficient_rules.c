/* test_coefficient_rules.c — MonomialList, CoefficientRules, FromCoefficientRules
 *
 * Mirrors the test_poly.c convention: parse -> evaluate -> compare FullForm
 * string. Uses ASSERT_MSG (never libc assert; Release defines NDEBUG). Every
 * driver call frees the printed string and both Expr trees.
 */
#include "poly.h"
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

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* ------------------------------------------------------------------ */

static void test_coefficientrules_basic(void) {
    run_test("CoefficientRules[(x+y)^3]",
             "List[Rule[List[3, 0], 1], Rule[List[2, 1], 3], "
             "Rule[List[1, 2], 3], Rule[List[0, 3], 1]]");
    run_test("CoefficientRules[a x^2+b x y+c y^2,{x,y}]",
             "List[Rule[List[2, 0], a], Rule[List[1, 1], b], Rule[List[0, 2], c]]");
    /* not explicitly expanded */
    run_test("CoefficientRules[(x+1)(x+2),{x}]",
             "List[Rule[List[2], 1], Rule[List[1], 3], Rule[List[0], 2]]");
    /* constant, no variables -> empty exponent vector */
    run_test("CoefficientRules[7]", "List[Rule[List[], 7]]");
    /* zero polynomial -> empty */
    run_test("CoefficientRules[0]", "List[]");
    /* symbolic coefficients merge across like monomials */
    run_test("CoefficientRules[a x + b x, {x}]",
             "List[Rule[List[1], Plus[a, b]]]");
    run_test("CoefficientRules[3 x^2 - 3 x^2 + y, {x, y}]",
             "List[Rule[List[0, 1], 1]]");
}

static void test_monomiallist_basic(void) {
    run_test("MonomialList[(x+y)^3]",
             "List[Power[x, 3], Times[3, Power[x, 2], y], "
             "Times[3, x, Power[y, 2]], Power[y, 3]]");
    run_test("MonomialList[a x^2+b x y+c y^2,{x,y}]",
             "List[Times[a, Power[x, 2]], Times[b, x, y], Times[c, Power[y, 2]]]");
    /* Plus @@ MonomialList reconstructs the (expanded) polynomial */
    run_test("Plus @@ MonomialList[a x^2+b x y+c y^2,{x,y}] == Expand[a x^2+b x y+c y^2]",
             "True");
}

static void test_orders(void) {
    /* DegreeReverseLexicographic, as a string and as the equivalent matrix */
    run_test("CoefficientRules[a x y^2+b x^2 z,{x,y,z},\"DegreeReverseLexicographic\"]",
             "List[Rule[List[1, 2, 0], a], Rule[List[2, 0, 1], b]]");
    run_test("CoefficientRules[a x y^2+b x^2 z,{x,y,z},{{1,1,1},{0,0,-1},{0,-1,0}}] "
             "=== CoefficientRules[a x y^2+b x^2 z,{x,y,z},\"DegreeReverseLexicographic\"]",
             "True");
    run_test("MonomialList[x^2 y^2+x^3,{x,y},\"DegreeLexicographic\"]",
             "List[Times[Power[x, 2], Power[y, 2]], Power[x, 3]]");
    run_test("MonomialList[x^2 y^2+x^3,{x,y},{{1,1},{1,0}}] "
             "=== MonomialList[x^2 y^2+x^3,{x,y},\"DegreeLexicographic\"]",
             "True");
    /* NegativeLexicographic == Sort of exponent vectors (ascending) */
    run_test("CoefficientRules[x^2 - 1, x, \"NegativeLexicographic\"]",
             "List[Rule[List[0], -1], Rule[List[2], 1]]");
    /* For two variables DegreeLexicographic and DegreeReverseLexicographic coincide */
    run_test("CoefficientRules[3+10 x+2 x^2-3 x^5 y+9 x^10 y-6 x^4 y^2-2 x y^4-10 y^5,{x,y},"
             "\"DegreeLexicographic\"] === CoefficientRules[3+10 x+2 x^2-3 x^5 y+9 x^10 y"
             "-6 x^4 y^2-2 x y^4-10 y^5,{x,y},\"DegreeReverseLexicographic\"]",
             "True");
    /* NegativeDegreeReverseLexicographic == reversed-vars DegreeLexicographic, reversed */
    run_test("MonomialList[-x^3 y z-3 x^9 y^7 z+6 x^10 y^8 z-5 x y^10 z-x^2 y^5 z^4,"
             "{x,y,z},\"NegativeDegreeReverseLexicographic\"] === Reverse@MonomialList["
             "-x^3 y z-3 x^9 y^7 z+6 x^10 y^8 z-5 x y^10 z-x^2 y^5 z^4,Reverse@{x,y,z},"
             "\"DegreeLexicographic\"]",
             "True");
}

static void test_modulus(void) {
    run_test("CoefficientRules[(x+1)^5,x,Modulus->2]",
             "List[Rule[List[5], 1], Rule[List[4], 1], Rule[List[1], 1], Rule[List[0], 1]]");
    run_test("MonomialList[(x+1)^5,x,Modulus->2]",
             "List[Power[x, 5], Power[x, 4], x, 1]");
    /* Modulus reducing every coefficient to zero */
    run_test("CoefficientRules[(x+1)^2, x, Modulus->1]", "List[]");
}

static void test_fromcoefficientrules(void) {
    run_test("FromCoefficientRules[{{2,0}->a,{1,1}->b,{0,2}->c},{x,y}] "
             "== Expand[a x^2 + b x y + c y^2]", "True");
    run_test("FromCoefficientRules[{{1,0,0}->1,{0,1,1}->1},{y,x,z}]",
             "Plus[y, Times[x, z]]");
    run_test("FromCoefficientRules[{},{x,y}]", "0");
    /* full round trip */
    run_test("FromCoefficientRules[CoefficientRules[(x+2 y)^3,{x,y}],{x,y}] "
             "== Expand[(x+2 y)^3]", "True");
    /* exponent-vector length mismatch -> unevaluated */
    run_test("FromCoefficientRules[{{2}->a},{x,y}]",
             "FromCoefficientRules[List[Rule[List[2], a]], List[x, y]]");
}

static void test_decline(void) {
    /* a bare List is not a polynomial -> unevaluated (symbolic/structural head) */
    run_test("MonomialList[{1,2,3}]", "MonomialList[List[1, 2, 3]]");
    run_test("CoefficientRules[{1,2,3}]", "CoefficientRules[List[1, 2, 3]]");
    /* unknown order string -> unevaluated */
    run_test("CoefficientRules[x, x, \"Nonsense\"]",
             "CoefficientRules[x, x, \"Nonsense\"]");
    /* All is equivalent to the omitted-variable form Variables[poly] (which
     * here treats a, b, c as variables too, so both give the same 5-vector rules) */
    run_test("CoefficientRules[a x^2+b x y+c y^2, All] "
             "=== CoefficientRules[a x^2+b x y+c y^2]", "True");
    run_test("CoefficientRules[a x^2+b x y+c y^2, All]",
             "List[Rule[List[1, 0, 0, 2, 0], 1], Rule[List[0, 1, 0, 1, 1], 1], "
             "Rule[List[0, 0, 1, 0, 2], 1]]");
}

int main(void) {
    setbuf(stdout, NULL);
    printf("Starting coefficient_rules_tests\n");
    symtab_init();
    core_init();

    TEST(test_coefficientrules_basic);
    TEST(test_monomiallist_basic);
    TEST(test_orders);
    TEST(test_modulus);
    TEST(test_fromcoefficientrules);
    TEST(test_decline);

    printf("All MonomialList / CoefficientRules / FromCoefficientRules tests passed!\n");
    return 0;
}
