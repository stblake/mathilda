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

/* Numeric subresultant polynomials, including a defective (degree-gap)
 * chain where one entry vanishes and another has degree below its index. */
void test_subrespoly_numeric() {
    run_test("SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]",
             "List[0, Plus[-36, Times[36, x]], "
             "Plus[38, Times[-49, x], Times[11, Power[x, 2]]], "
             "Plus[-16, Times[24, x], Times[-9, Power[x, 2]], Power[x, 3]]]");
    run_test("SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]",
             "List[-183782157189, Plus[-761749829, Times[3208696817, x]], "
             "Plus[-3143546, Times[11222638, x], Times[3838135, Power[x, 2]]], "
             "Plus[-21609, Times[163611, x], Times[-49392, Power[x, 2]], Times[64827, Power[x, 3]]], "
             "0, "
             "Plus[-49, Times[371, x], Times[-112, Power[x, 2]], Times[147, Power[x, 3]]], "
             "Plus[-9, Times[8, x], Times[7, Power[x, 6]]]]");
}

/* Coefficients of the subresultant polynomials are polynomials in the
 * coefficients of the inputs. */
void test_subrespoly_symbolic() {
    run_test("SubresultantPolynomials[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]",
             "List[Plus[Times[4, Times[Power[a, 2], Power[c, 3]]], "
             "Times[2, Times[a, Power[b, 3], d]], "
             "Times[-18, Times[Power[a, 2], b, c, d]], "
             "Times[27, Times[Power[a, 3], Power[d, 2]]]], "
             "Plus[Times[-2, Times[a, b, c]], Times[9, Times[Power[a, 2], d]], "
             "Times[-2, Times[a, Power[b, 2], x]], Times[6, Times[Power[a, 2], c, x]]], "
             "Plus[c, Times[b, x], Times[3, Times[a, Power[x, 2]]]]]");
}

/* The first subresultant polynomial equals the Resultant. */
void test_subrespoly_resultant() {
    run_test("First[SubresultantPolynomials[2x^5+3x^3-7x^2+11x+21,3x^4+27x^3-11x+9,x]] "
             "- Resultant[2x^5+3x^3-7x^2+11x+21,3x^4+27x^3-11x+9,x]", "0");
    run_test("First[SubresultantPolynomials[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]] "
             "- Resultant[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]", "0");
}

/* The list length is Exponent[poly2, var] + 1. */
void test_subrespoly_length() {
    run_test("Length[SubresultantPolynomials[x^50 + a, x^20 + b, x]]", "21");
    run_test("Length[SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]]", "7");
    run_test("Length[SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]]", "4");
}

/* The coefficient of var^(i-1) in S_(i-1) is the i-th principal subresultant
 * coefficient, i.e. Subresultants[poly1, poly2, var]. */
void test_subrespoly_psc() {
    run_test("Table[Coefficient[SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x][[i]], "
             "x, i - 1], {i, 7}] - Subresultants[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]",
             "List[0, 0, 0, 0, 0, 0, 0]");
    run_test("Table[Coefficient[SubresultantPolynomials[-14 + 4 x + 3 x^2 + 2 x^3 + 5 x^8, "
             "-21 + 7 x + 6 x^2 + 5 x^3 + 3 x^7, x][[i]], x, i - 1], {i, 8}] "
             "- Subresultants[-14 + 4 x + 3 x^2 + 2 x^3 + 5 x^8, -21 + 7 x + 6 x^2 + 5 x^3 + 3 x^7, x]",
             "List[0, 0, 0, 0, 0, 0, 0, 0]");
}

/* The last subresultant polynomial is poly2 when deg(poly1) = deg(poly2) + 1. */
void test_subrespoly_last() {
    run_test("Last[SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]] "
             "- (7 x^6 + 8 x - 9)", "0");
    run_test("Last[SubresultantPolynomials[a x^4 + b x + c, a c x^3 + b^2 x + a - b, x]] "
             "- (a c x^3 + b^2 x + a - b)", "0");
}

/* Degenerate inputs: shared factors (leading zeros), equal degrees, and a
 * constant second argument. */
void test_subrespoly_degenerate() {
    /* GCD x^2 - 1 of degree 2 -> first two entries vanish, last is the GCD. */
    run_test("SubresultantPolynomials[x^4 - 1, x^2 - 1, x]",
             "List[0, 0, Plus[-1, Power[x, 2]]]");
    /* Equal degrees: top entry is the empty (0x0) minor, 1. */
    run_test("SubresultantPolynomials[2 x^2 + 3 x + 1, 5 x^2 + x + 7, x]",
             "List[341, Plus[9, Times[-13, x]], 1]");
    run_test("SubresultantPolynomials[2 x^2 + 3 x + 1, 5 x^2 + x + 7, x][[1]] "
             "- Resultant[2 x^2 + 3 x + 1, 5 x^2 + x + 7, x]", "0");
    /* Constant second argument -> Resultant[poly1, c] = c^deg(poly1). */
    run_test("SubresultantPolynomials[x^3 + 2 x + 1, 5, x]", "List[125]");
}

/* Algebraic-number coefficients route through the determinant-polynomial
 * fallback (the PRS path is skipped to avoid chain bloat). */
void test_subrespoly_algebraic() {
    run_test("SubresultantPolynomials[x^2 + Sqrt[2] x + 1, x + 1, x]",
             "List[Plus[2, Times[-1, Power[2, Rational[1, 2]]]], Plus[1, x]]");
}

/* Inexact coefficients or deg(poly1) < deg(poly2): emit SubresultantPolynomials
 * ::npolys and leave the call unevaluated. */
void test_subrespoly_errors() {
    run_test("SubresultantPolynomials[x^2 - 1.2, x - 3.4, x]",
             "SubresultantPolynomials[Plus[-1.2, Power[x, 2]], Plus[-3.4, x], x]");
    run_test("SubresultantPolynomials[x + 1, x^2 + 1, x]",
             "SubresultantPolynomials[Plus[1, x], Plus[1, Power[x, 2]], x]");
}

int main() {
    setbuf(stdout, NULL);
    printf("Starting subresultantpolynomials_tests\n");
    symtab_init();
    core_init();

    TEST(test_subrespoly_numeric);
    TEST(test_subrespoly_symbolic);
    TEST(test_subrespoly_resultant);
    TEST(test_subrespoly_length);
    TEST(test_subrespoly_psc);
    TEST(test_subrespoly_last);
    TEST(test_subrespoly_degenerate);
    TEST(test_subrespoly_algebraic);
    TEST(test_subrespoly_errors);

    printf("All subresultantpolynomials tests passed!\n");
    return 0;
}
