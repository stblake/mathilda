/*
 * test_factor_options.c
 * ---------------------
 * The options of Factor and its family: Modulus (GF(p), any prime size,
 * univariate and multivariate), GaussianIntegers, option validation and
 * messages, and forwarding through FactorList / FactorSquareFree(List).
 * Messages are checked through Check[], which returns its second argument
 * exactly when a message was issued.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static int failures = 0;

/* Evaluate `input` and compare its FullForm with `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL (parse): %s\n", input);
        failures++;
        return;
    }
    Expr* r = evaluate(e);
    char* got = expr_to_string_fullform(r);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        failures++;
    } else {
        printf("PASS: %s\n", input);
    }
    free(got);
    expr_free(r);
    expr_free(e);
}

static void check_true(const char* input) { check(input, "True"); }

static void check_msg(const char* input) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Check[%s, $Failed]", input);
    check(buf, "$Failed");
}

static void test_modulus(void) {
    check("Factor[3 x^2 + 3, Modulus -> 5]",
          "Times[3, Plus[2, x], Plus[3, x]]");
    check("Factor[2 x^2 + 4 x + 2, Modulus -> 7]", "Times[2, Power[Plus[1, x], 2]]");
    check("Factor[x/2 + 1, Modulus -> 5]", "Times[3, Plus[2, x]]");
    check("Factor[x^2 + 1, Modulus -> -5]", "Times[Plus[2, x], Plus[3, x]]");
    check("Factor[5, Modulus -> 3]", "2");
    check("Factor[4 x^2 + 4, Modulus -> 2]", "0");
    check("Factor[x^4 - y^4, Modulus -> 3]",
          "Times[Plus[x, y], Plus[x, Times[2, y]], Plus[Power[x, 2], Power[y, 2]]]");
    check("Factor[y^2 x^2 - 1, Modulus -> 7]",
          "Times[Plus[1, Times[x, y]], Plus[6, Times[x, y]]]");
    check("Factor[x^4 + 1, Modulus -> 18446744073709551629]",
          "Times[Plus[2370518075556110396, Power[x, 2]], "
          "Plus[16076225998153441233, Power[x, 2]]]");
    check("Factor[(x^2 - 1)/(x + 3), Modulus -> 5]",
          "Times[Plus[1, x], Power[Plus[3, x], -1], Plus[4, x]]");
    check("Factor[x^2 - 1 == 0, Modulus -> 3]",
          "Equal[Times[Plus[1, x], Plus[2, x]], 0]");
    check("Factor[{x^2 + 1, x^3 + 1}, Modulus -> 2]",
          "List[Power[Plus[1, x], 2], Times[Plus[1, x], Plus[1, x, Power[x, 2]]]]");
    check("Factor[3 (2 x + 2)^2 (x^2 + 1), Modulus -> 5]",
          "Times[2, Power[Plus[1, x], 2], Plus[2, x], Plus[3, x]]");
    check("Factor[x^4 + 1, Modulus -> 0]", "Plus[1, Power[x, 4]]");
    check("Factor[x^2 - 1, Modulus -> a]", "Times[Plus[-1, x], Plus[1, x]]");
    /* never expanded over Q: (1 + x)^32768 has enormous binomial coefficients */
    check("Factor[(x + 1)^32768 (x^2 + x + 1), Modulus -> 2]",
          "Times[Power[Plus[1, x], 32768], Plus[1, x, Power[x, 2]]]");
    check("FactorList[3 x^2 + 3, Modulus -> 5]",
          "List[List[3, 1], List[Plus[2, x], 1], List[Plus[3, x], 1]]");
    check("FactorSquareFree[x^3 + x^2, Modulus -> 2]",
          "Times[Power[x, 2], Plus[1, x]]");
    check("FactorSquareFreeList[x^3 + x^2, Modulus -> 2]",
          "List[List[1, 1], List[x, 2], List[Plus[1, x], 1]]");
}

static void test_option_messages(void) {
    check_msg("Factor[x^2 + 1, Modulus -> 9]");                 /* modp */
    check_msg("Factor[x^2 + 1, Modulus -> 1]");                 /* modp */
    check_msg("Factor[x^2 + 1, Modulus -> 2.5]");               /* bmod */
    check_msg("Factor[x^2 + 1, Modulus -> 2, Extension -> Sqrt[2]]");   /* mext */
    check_msg("Factor[x^2 + 1, GaussianIntegers -> True, Modulus -> 3]"); /* facim */
    check_msg("Factor[x^2 + I, Modulus -> 5]");                 /* facim */
    check_msg("Factor[x/5 + 1, Modulus -> 5]");                 /* nmod */
    check_msg("Factor[x^2 + 1, Method -> \"Foo\"]");            /* optx */
    check_msg("Factor[x^2 - 1, 3]");                            /* nonopt */
    check_msg("FactorList[x^2 + 1, Modulus -> 4]");
    check_msg("FactorList[x^2 - 1, Foo -> 1]");
    check_msg("FactorSquareFree[x^2 - 1, GaussianIntegers -> True]");
    check("Head[Quiet[Factor[x^2 + 1, Modulus -> 9]]]", "Factor");
    check("Quiet[FactorList[x^2 + 1, Modulus -> 4]]",
          "FactorList[Plus[1, Power[x, 2]], Rule[Modulus, 4]]");
    /* first occurrence of an option wins */
    check("Factor[x^2 - 1, Modulus -> 3, Modulus -> 5]",
          "Times[Plus[1, x], Plus[2, x]]");
}

/* Options at their defaults are the plain form, with every fast path:
 * Extension -> None used to push x^99 - y^99 onto the classical path
 * (over 100 s). */
static void test_default_options_plain(void) {
    check_true("Length[FactorList[x^99 - y^99, Extension -> None]] == 7");
    check("Factor[x^2 - 1, GaussianIntegers -> False, Trig -> False]",
          "Times[Plus[-1, x], Plus[1, x]]");
    check("Factor[x^2 + 1, GaussianIntegers -> True]",
          "Times[Plus[Complex[0, -1], x], Plus[Complex[0, 1], x]]");
    check("FactorList[x^2 + 1, GaussianIntegers -> True]",
          "List[List[1, 1], List[Plus[Complex[0, -1], x], 1], List[Plus[Complex[0, 1], x], 1]]");
    check("Options[Factor]",
          "List[Rule[Extension, None], Rule[GaussianIntegers, False], "
          "Rule[Modulus, 0], Rule[Trig, False]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor option tests...\n");
    test_modulus();
    test_option_messages();
    test_default_options_plain();

    if (failures) {
        printf("%d factor option test(s) FAILED\n", failures);
        return 1;
    }
    printf("All factor option tests passed!\n");
    return 0;
}
