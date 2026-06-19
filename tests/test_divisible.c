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
        printf("Divisible test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Machine integers ---- */
void test_divisible_machine_integers() {
    check("Divisible[10, 2]", "True");
    check("Divisible[5, 2]", "False");
    check("Divisible[6, 3]", "True");
    check("Divisible[48, 8]", "True");
    check("Divisible[7, 7]", "True");
    check("Divisible[7, 1]", "True");
    check("Divisible[1, 7]", "False");
    check("Divisible[100, 25]", "True");
    check("Divisible[101, 25]", "False");
}

/* ---- Negative arguments: divisibility ignores sign ---- */
void test_divisible_negative() {
    check("Divisible[10, -2]", "True");
    check("Divisible[-10, 2]", "True");
    check("Divisible[-10, -3]", "False");
    check("Divisible[-9, 3]", "True");
}

/* ---- Zero edge cases (GMP semantics: divisible by 0 iff n == 0) ---- */
void test_divisible_zero() {
    check("Divisible[0, 5]", "True");
    check("Divisible[0, 0]", "True");
    check("Divisible[6, 0]", "False");
}

/* ---- BigInt integers and large-integer divisibility ---- */
void test_divisible_bignum() {
    check("Divisible[10^3000 + 1, 16001]", "True");
    check("Divisible[10^3000 + 2, 16001]", "False");
    check("Divisible[2^200, 2^199]", "True");
    check("Divisible[2^200, 2^201]", "False");
    check("Divisible[3^500, 3]", "True");
    /* BigInt divided by a BigInt that does not divide it. */
    check("Divisible[100000000000000000000000001, 7]", "False");
    check("Divisible[7 * 100000000000000000000000001, 7]", "True");
}

/* ---- Gaussian integers ---- */
void test_divisible_gaussian() {
    check("Divisible[3 + I, 1 - I]", "True");   /* quotient 1 + 2 I */
    check("Divisible[4 + 2 I, 2]", "True");      /* quotient 2 + I */
    check("Divisible[5, 1 + 2 I]", "True");      /* 5 = (1+2I)(1-2I) */
    check("Divisible[3 + I, 2]", "False");       /* quotient 3/2 + I/2 */
    check("Divisible[1 + I, 3]", "False");
}

/* ---- Rationals ---- */
void test_divisible_rationals() {
    check("Divisible[3/2, 1/2]", "True");        /* quotient 3 */
    check("Divisible[2/3, 1/3]", "True");        /* quotient 2 */
    check("Divisible[1/2, 1/3]", "False");       /* quotient 3/2 */
    check("Divisible[5/2, 1/4]", "True");        /* quotient 10 */
    check("Divisible[6, 3/2]", "True");          /* quotient 4 */
}

/* ---- Symbolic forms of numeric quantities ---- */
void test_divisible_symbolic_numeric() {
    check("Divisible[2 Pi, Pi/2]", "True");      /* quotient 4 */
    check("Divisible[4 Pi/3, 2 Pi/3]", "True");  /* quotient 2 */
    check("Divisible[Pi, Pi]", "True");
    check("Divisible[3 Pi, Pi]", "True");
    check("Divisible[Pi, 2 Pi]", "False");       /* quotient 1/2 */
}

/* ---- Numeric quantities that are not divisible ---- */
void test_divisible_numeric_quantities() {
    check("Divisible[Sqrt[6], Sqrt[2]]", "False");  /* quotient Sqrt[3] */
    check("Divisible[Sqrt[8], Sqrt[2]]", "True");   /* quotient 2 */
    check("Divisible[Sqrt[18], Sqrt[2]]", "True");  /* quotient 3 */
    check("Divisible[Sqrt[12], Sqrt[3]]", "True");  /* quotient 2 */
    check("Divisible[Sqrt[5], Sqrt[2]]", "False");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_divisible_listable() {
    check("Divisible[{1, 2, 3, 4, 5, 6}, 2]",
          "{False, True, False, True, False, True}");
    check("Divisible[12, {2, 3, 4, 5, 6}]",
          "{True, True, True, False, True}");
    check("Divisible[{10, 15, 20}, {2, 5, 3}]",
          "{True, True, False}");
    /* Threads recursively through nested lists. */
    check("Divisible[{{4, 5}, {6, 7}}, 2]",
          "{{True, False}, {True, False}}");
}

/* ---- Symbolic, non-numeric arguments are left unevaluated ---- */
void test_divisible_symbolic_unevaluated() {
    check("Divisible[x, 2]", "Divisible[x, 2]");
    check("Divisible[2, x]", "Divisible[2, x]");
    check("Divisible[x, y]", "Divisible[x, y]");
    check("Divisible[Pi, x]", "Divisible[Pi, x]");
    check("Divisible[n + 1, n]", "Divisible[1 + n, n]");
}

/* ---- Wrong argument counts warn and leave the call unevaluated ---- */
void test_divisible_arg_errors() {
    check("Divisible[]", "Divisible[]");
    check("Divisible[7]", "Divisible[7]");
    check("Divisible[6, 5, 4, 3]", "Divisible[6, 5, 4, 3]");
}

/* ---- Cross-checks against related primitives ---- */
void test_divisible_consistency() {
    check("{IntegerQ[10/2], Divisible[10, 2]}", "{True, True}");
    check("{IntegerQ[10/3], Divisible[10, 3]}", "{False, False}");
    /* If n is divisible by m, GCD[n, m] == m (for positive m). */
    check("Divisible[48, 8] && (GCD[48, 8] == 8)", "True");
    /* Divisible[n, m] agrees with Mod[n, m] == 0 over the integers. */
    check("Divisible[42, 7] == (Mod[42, 7] == 0)", "True");
    check("Divisible[43, 7] == (Mod[43, 7] == 0)", "True");
}

/* ---- Attributes: Listable and Protected ---- */
void test_divisible_attributes() {
    check("Attributes[Divisible]", "{Listable, Protected}");
    check("MemberQ[Attributes[Divisible], Listable]", "True");
    check("MemberQ[Attributes[Divisible], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_divisible_machine_integers);
    TEST(test_divisible_negative);
    TEST(test_divisible_zero);
    TEST(test_divisible_bignum);
    TEST(test_divisible_gaussian);
    TEST(test_divisible_rationals);
    TEST(test_divisible_symbolic_numeric);
    TEST(test_divisible_numeric_quantities);
    TEST(test_divisible_listable);
    TEST(test_divisible_symbolic_unevaluated);
    TEST(test_divisible_arg_errors);
    TEST(test_divisible_consistency);
    TEST(test_divisible_attributes);

    printf("All Divisible tests passed!\n");
    return 0;
}
