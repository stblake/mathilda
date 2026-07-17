#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("Exponent test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic maximum exponent ---- */
void test_exponent_basic() {
    check("Exponent[1+x^2+a x^3, x]", "3");
    check("Exponent[1, x]", "0");
    check("Exponent[2 x^3, x]", "3");
    check("Exponent[a, x]", "0");           /* free of x */
    check("Exponent[x, x]", "1");
    check("Exponent[x^5, x]", "5");
}

/* ---- The genuine zero polynomial: empty exponent set ---- */
void test_exponent_zero() {
    check("Exponent[0, x]", "-Infinity");   /* Max[] */
    check("Exponent[0, x, Min]", "Infinity"); /* Min[] */
    check("Exponent[0, x, List]", "{}");
}

/* ---- Works on unexpanded input (Exponent expands first) ---- */
void test_exponent_expands() {
    check("Exponent[(x^2+1)^3+1, x]", "6");
    check("Exponent[(x^2+1)^3-1, x, Min]", "2");   /* constant term cancels */
    check("Exponent[(x+1)(x-1), x]", "2");
    check("Exponent[(1+x)^10, x]", "10");
}

/* ---- Rational and symbolic exponents ---- */
void test_exponent_symbolic() {
    check("Exponent[x^(n+1)+2 Sqrt[x]+1, x]", "Max[1/2, 1 + n]");
    check("Exponent[x^(1/2)+x^(3/2), x]", "3/2");
    check("Exponent[Sqrt[x], x]", "1/2");
    check("Exponent[x^(2 n), x]", "2 n");
}

/* ---- Laurent (negative) exponents ---- */
void test_exponent_laurent() {
    check("Exponent[1/x, x]", "-1");
    check("Exponent[1/x^3 + x, x, Min]", "-3");
    check("Exponent[1/x^3 + x, x]", "1");
}

/* ---- h = Min / List ---- */
void test_exponent_h() {
    check("Exponent[x + x^2, x, Min]", "1");
    check("Exponent[1+x^2+a x^3, x, List]", "{0, 2, 3}");
    check("Exponent[x^2 + a x^2 + x, x, List]", "{1, 2}"); /* like terms collapse */
    check("Exponent[x^4 + x^2, x, List]", "{2, 4}");
}

/* ---- Purely syntactic: no zero-coefficient recognition ---- */
void test_exponent_syntactic() {
    /* zero is numerically 0 but not in normal form; Exponent must not
     * recognise the vanishing coefficient. */
    check("Exponent[(Sqrt[2]+Sqrt[3]-Sqrt[5+2Sqrt[6]]) x^2 + x + 1, x]", "2");
    /* Once the coefficient is genuinely reduced away, the degree drops. */
    check("Exponent[x + 1, x]", "1");
}

/* ---- form is a kernel, not necessarily a bare symbol ---- */
void test_exponent_kernel_form() {
    check("Exponent[Sin[x], x]", "0");          /* nested, not a power of x */
    check("Exponent[Sin[x]^3 + Sin[x], Sin[x]]", "3");
    check("Exponent[Sin[x]^3 + Sin[x], Sin[x], Min]", "1");
    check("Exponent[Log[x]^2 + 1, Log[x]]", "2");
}

/* ---- form is a product of terms ---- */
void test_exponent_product_form() {
    check("Exponent[(x y)^2 + x y, x y]", "2");
    check("Exponent[x^2 y^2, x y]", "2");
    check("Exponent[x^2 y, x y]", "1");         /* min(2,1) */
    check("Exponent[a (x y)^3, x y]", "3");
}

/* ---- Listable: threading over a list of forms ---- */
void test_exponent_listable() {
    check("Exponent[1+x^2+a x^3, {x, a}]", "{3, 1}");
    check("Exponent[x^2 y^3, {x, y}]", "{2, 3}");
    check("Exponent[x^2 y^3 + x, {x, y}, Min]", "{1, 0}");
    check("MemberQ[Attributes[Exponent], Listable]", "True");
}

/* ---- Attributes ---- */
void test_exponent_attributes() {
    check("Attributes[Exponent]", "{Listable, Protected}");
    check("MemberQ[Attributes[Exponent], Protected]", "True");
}

/* ---- Argument errors: unevaluated + diagnostic ---- */
void test_exponent_arity() {
    /* argt: 0 or 1 args -> unevaluated (message to stderr). */
    check("Exponent[]", "Exponent[]");
    check("Exponent[x]", "Exponent[x]");
    /* nonopt: extra non-option args beyond position 3 -> unevaluated. */
    check("Exponent[1, x, 3, 4, 5]", "Exponent[1, x, 3, 4, 5]");
}

/* ---- A polynomial's degree via Exponent ---- */
void test_exponent_degree() {
    check("Exponent[x^6 + 3 x^4 + 3 x^2 + 2, x]", "6");
    check("Exponent[3 + 2 x - 5 x^2, x]", "2");
    check("Exponent[c, y]", "0");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_exponent_basic);
    TEST(test_exponent_zero);
    TEST(test_exponent_expands);
    TEST(test_exponent_symbolic);
    TEST(test_exponent_laurent);
    TEST(test_exponent_h);
    TEST(test_exponent_syntactic);
    TEST(test_exponent_kernel_form);
    TEST(test_exponent_product_form);
    TEST(test_exponent_listable);
    TEST(test_exponent_attributes);
    TEST(test_exponent_arity);
    TEST(test_exponent_degree);

    printf("All Exponent tests passed!\n");
    return 0;
}
