/* Tests for MangoldtLambda[n], the von Mangoldt function.
 * Semantics follow the reference: exact Log[p] on a prime power p^k, 0
 * otherwise (n = 1, n <= 0 included); MangoldtLambda::exact and unevaluated
 * on an inexact argument; unevaluated, silently, on anything else. */

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

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        printf("MangoldtLambda test failed: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    }
    free(s);
    expr_free(e);
    expr_free(r);
}

static void test_values(void) {
    check("MangoldtLambda[2]", "Log[2]");
    check("MangoldtLambda[8]", "Log[2]");
    check("MangoldtLambda[9]", "Log[3]");
    check("MangoldtLambda[49]", "Log[7]");
    check("MangoldtLambda[1]", "0");
    check("MangoldtLambda[6]", "0");
    check("MangoldtLambda[12]", "0");
    check("MangoldtLambda[0]", "0");
    check("MangoldtLambda[-8]", "0");
    check("Table[MangoldtLambda[n], {n, 1, 10}]",
          "{0, Log[2], Log[3], Log[2], Log[5], 0, Log[7], Log[2], Log[3], 0}");
    check("N[MangoldtLambda[9]]", "1.09861");
}

static void test_bignum(void) {
    check("MangoldtLambda[(10^20 + 39)^3]", "Log[100000000000000000039]");
    check("MangoldtLambda[2^127 - 1]", "Log[170141183460469231731687303715884105727]");
    check("MangoldtLambda[(2^61 - 1) (2^31 - 1)]", "0");
}

/* Sum_{d | n} Lambda(d) = Log n */
static void test_identity(void) {
    check("Simplify[PowerExpand[Total[MangoldtLambda /@ Divisors[360]] - Log[360]]]", "0");
}

static void test_unevaluated(void) {
    check("Quiet[MangoldtLambda[2.]]", "MangoldtLambda[2.0]");
    check("Check[MangoldtLambda[2.], fired]", "fired");
    check("Check[MangoldtLambda[x], fired]", "MangoldtLambda[x]");
    check("MangoldtLambda[9/2]", "MangoldtLambda[9/2]");
    check("MangoldtLambda[2 I]", "MangoldtLambda[2*I]");
}

static void test_attributes(void) {
    check("Attributes[MangoldtLambda]", "{Listable, Protected}");
    check("MangoldtLambda[{2, 3, 4, 6}]", "{Log[2], Log[3], Log[2], 0}");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_values);
    TEST(test_bignum);
    TEST(test_identity);
    TEST(test_unevaluated);
    TEST(test_attributes);
    printf("All MangoldtLambda tests passed!\n");
    return 0;
}
