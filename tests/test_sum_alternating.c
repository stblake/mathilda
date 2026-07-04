/*
 * test_sum_alternating.c -- Sum`Alternating: infinite alternating rational sums.
 *
 * check() compares the InputForm string against an expected string; same()
 * cross-checks a symbolic closed form via Simplify[a - b] == 0.
 */

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
static int checks = 0;

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "PARSE FAIL: %s\n", input); exit(1); }
    Expr* v = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(v);
    expr_free(v);
    return s;
}

static void check(const char* input, const char* expected) {
    checks++;
    char* got = eval_str(input);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        failures++;
    }
    free(got);
}

static void same(const char* a, const char* b) {
    checks++;
    size_t n = strlen(a) + strlen(b) + 32;
    char* buf = malloc(n);
    snprintf(buf, n, "Simplify[(%s) - (%s)]", a, b);
    char* s = eval_str(buf);
    if (strcmp(s, "0") != 0) {
        fprintf(stderr, "FAIL (oracle): %s  vs  %s\n  Simplify[diff] = %s\n", a, b, s);
        failures++;
    }
    free(s); free(buf);
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- Headline elementary closed forms ---- */
    check("Sum[(-1)^(k + 1)/k, {k, 1, Infinity}]", "Log[2]");
    same("Sum[(-1)^k/(2 k + 1), {k, 1, Infinity}]", "Pi/4 - 1");
    check("Sum[(-1)^k/(2 k + 1)^2, {k, 0, Infinity}]", "Catalan");

    /* ---- More of the same families ---- */
    check("Sum[(-1)^k/k, {k, 1, Infinity}]", "-Log[2]");
    same("Sum[(-1)^k/k^2, {k, 1, Infinity}]", "-Pi^2/12");
    same("Sum[(-1)^k/(2 k + 1)^3, {k, 0, Infinity}]", "Pi^3/32");     /* Dirichlet beta(3) */
    same("Sum[(-1)^k/(2 k + 1), {k, 0, Infinity}]", "Pi/4");           /* Leibniz */
    same("Sum[(-1)^(k + 1)/(2 k + 1), {k, 0, Infinity}]", "-Pi/4");    /* negated Leibniz */

    /* numeric cross-check for a general linear pole (a = 1/3) */
    check("Chop[N[Sum[(-1)^k/(3 k + 1), {k, 0, Infinity}] "
          "- Sum[(-1)^k/(3 k + 1), {k, 0, 5000}], 15], 10^-3]", "0");

    /* ---- Non-alternating / divergent / complex-pole inputs stay held ---- */
    check("Sum[1/k^2, {k, 1, Infinity}]", "1/6 Pi^2");                 /* not intercepted */
    check("Sum[(-1)^k k, {k, 1, Infinity}]",
          "Sum[(-1)^k k, {k, 1, Infinity}]");                          /* divergent */
    check("Sum[(-1)^k/(k^2 + 1), {k, 1, Infinity}]",
          "Sum[(-1)^k/(k^2 + 1), {k, 1, Infinity}]");                  /* complex poles */

    if (failures) {
        fprintf(stderr, "\n%d/%d sum_alternating checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum_alternating tests passed!\n", checks);
    return 0;
}
