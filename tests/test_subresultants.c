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

void run_test(const char* input, const char* expected) {
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

/* Numeric principal subresultant coefficients (the Wolfram reference
 * outputs), including defective (degree-gap) chains where leading PSCs
 * vanish. */
void test_subresultants_numeric() {
    run_test("Subresultants[2x^7+3x^3-7x+1,3x^5-17x+21,x]",
             "List[273612691817, 68946901, 1299537, 16641, 0, 9]");
    run_test("Subresultants[(x-1)(x-2)(x-3),(x-1)(x-4)(x-5),x]",
             "List[0, 12, -4, 1]");
    run_test("Subresultants[(x-1)^5(x-2)(x-3),(x-1)^4(x-4)(x-5),x]",
             "List[0, 0, 0, 0, 144, 18, 1]");
    run_test("Subresultants[(x-1)(x-2)^2(x-3)^3,(x-1)(x-2)^2(x-4)(x-5)^2,x]",
             "List[0, 0, 0, -64, 44, -5, 1]");
    /* Three common roots -> first three PSCs vanish. */
    run_test("Subresultants[(x-1)(x-2)(x-3),(x-1)(x-4)(x-5),x]",
             "List[0, 12, -4, 1]");
}

/* PSCs are polynomials in the coefficients of the inputs. */
void test_subresultants_symbolic() {
    run_test("Subresultants[a x^3+b x^2+c x+d,x^3-5b x-7a,x]",
             "List[Plus[Times[-343, Power[a, 6]], Times[245, Times[Power[a, 4], Power[b, 2]]], "
             "Times[-49, Times[Power[a, 2], Power[b, 3]]], Times[147, Times[Power[a, 3], b, c]], "
             "Times[-175, Times[Power[a, 3], Power[b, 2], c]], Times[35, Times[a, Power[b, 3], c]], "
             "Times[-70, Times[Power[a, 2], b, Power[c, 2]]], Times[-7, Times[a, Power[c, 3]]], "
             "Times[-147, Times[Power[a, 4], d]], Times[-35, Times[Power[a, 2], Power[b, 2], d]], "
             "Times[125, Times[Power[a, 2], Power[b, 3], d]], Times[-25, Times[Power[b, 4], d]], "
             "Times[21, Times[a, b, c, d]], Times[50, Times[a, Power[b, 2], c, d]], "
             "Times[5, Times[b, Power[c, 2], d]], Times[-21, Times[Power[a, 2], Power[d, 2]]], "
             "Times[-10, Times[Power[b, 2], Power[d, 2]]], Times[-1, Power[d, 3]]], "
             "Plus[Times[-7, Times[Power[a, 2], b]], Times[25, Times[Power[a, 2], Power[b, 2]]], "
             "Times[-5, Power[b, 3]], Times[10, Times[a, b, c]], Power[c, 2], Times[-1, Times[b, d]]], "
             "Times[-1, b], 1]");
    /* One pair of equal roots: first PSC disappears. */
    run_test("Subresultants[(x-a)(x-b)(x-c),(x-1)(x-2)(x-3),x] /. {a->1}",
             "List[0, Plus[36, Times[-30, b], Times[6, Power[b, 2]], Times[-30, c], "
             "Times[25, Times[b, c]], Times[-5, Times[Power[b, 2], c]], Times[6, Power[c, 2]], "
             "Times[-5, Times[b, Power[c, 2]]], Times[Power[b, 2], Power[c, 2]]], "
             "Plus[-5, b, c], 1]");
    /* Two pairs of equal roots: first two PSCs disappear. */
    run_test("Subresultants[(x-a)(x-b)(x-c),(x-1)(x-2)(x-3),x] /. {a->1,b->2}",
             "List[0, 0, Plus[-3, c], 1]");
}

/* The first element of Subresultants equals Resultant. */
void test_subresultants_resultant() {
    run_test("Subresultants[2x^5+3x^3-7x^2+11x+21,-3x^4+27x^3-11x+9,x]",
             "List[-68252080236, 265691684, 350853, 1485, -3]");
    run_test("Subresultants[2x^5+3x^3-7x^2+11x+21,-3x^4+27x^3-11x+9,x][[1]] "
             "- Resultant[2x^5+3x^3-7x^2+11x+21,-3x^4+27x^3-11x+9,x]", "0");
    run_test("First[Subresultants[a x^3+b x^2+c x+d,x^3-5b x-7a,x]] "
             "- Resultant[a x^3+b x^2+c x+d,x^3-5b x-7a,x]", "0");
}

/* The length is Min[Exponent[p1,var], Exponent[p2,var]] + 1. */
void test_subresultants_length() {
    run_test("Length[Subresultants[x^50+a,x^20+b,x]]", "21");
    run_test("Length[Subresultants[2x^7+3x^3-7x+1,3x^5-17x+21,x]]", "6");
    run_test("Length[Subresultants[(x-1)(x-2)(x-3),(x-1)(x-4)(x-5),x]]", "4");
}

/* Argument order: the swap multiplies element j by (-1)^((n-j)(m-j)); the
 * first element tracks the sign of Resultant. */
void test_subresultants_swap() {
    run_test("Subresultants[x^3+1, x+2, x]", "List[7, 1]");
    run_test("Subresultants[x+2, x^3+1, x]", "List[-7, 1]");      /* n*m odd -> sign flip */
    run_test("Subresultants[x^2+2, x^5+x+1, x]", "List[51, 5, 1]");
    run_test("Subresultants[x^5+x+1, x^2+2, x]", "List[51, 5, 1]");/* n*m even -> same */
}

/* Degenerate inputs: constants and one constant argument. */
void test_subresultants_degenerate() {
    run_test("Subresultants[7, 4, x]", "List[1]");                /* both constant */
    run_test("Subresultants[x+1, 5, x]", "List[5]");              /* Resultant[x+1,5] = 5 */
    run_test("Subresultants[x+1, 5, x][[1]] - Resultant[x+1,5,x]", "0");
    /* Identical polynomials -> GCD is the polynomial itself; only top PSC. */
    run_test("Subresultants[x^3+x+1, x^3+x+1, x]", "List[0, 0, 0, 1]");
    run_test("Subresultants[2x^2+3x+1, 2x^2+3x+1, x]", "List[0, 0, 1]");
    /* Shared factor x^2-1 of degree 2 -> first two PSCs vanish. */
    run_test("Subresultants[x^4-1, x^2-1, x]", "List[0, 0, 1]");
}

/* Algebraic-number coefficients route through the Sylvester-minor
 * determinant fallback (PRS is skipped to avoid chain bloat). */
void test_subresultants_algebraic() {
    run_test("Subresultants[x^2 + Sqrt[2] x + 1, x + 1, x]",
             "List[Plus[2, Times[-1, Power[2, Rational[1, 2]]]], 1]");
    run_test("Subresultants[x^2 - 2, x^2 - Sqrt[2] x, x]",
             "List[0, Times[-1, Power[2, Rational[1, 2]]], 1]");
}

int main() {
    setbuf(stdout, NULL);
    printf("Starting subresultants_tests\n");
    symtab_init();
    core_init();

    TEST(test_subresultants_numeric);
    TEST(test_subresultants_symbolic);
    TEST(test_subresultants_resultant);
    TEST(test_subresultants_length);
    TEST(test_subresultants_swap);
    TEST(test_subresultants_degenerate);
    TEST(test_subresultants_algebraic);

    printf("All subresultants tests passed!\n");
    return 0;
}
