/*
 * test_sum_rational.c -- Sum`Rational (Stage 5): infinite rational-function
 * summation via Hurwitz-Zeta / digamma / Coth closed forms.
 *
 * check() compares the InputForm of the evaluated input against an expected
 * string; same() asserts Simplify[a - b] == 0 (robust to printed form);
 * nclose() asserts two numeric values agree to a tolerance, computed entirely
 * inside the language (no parsing of printed reals).
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
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        failures++;
    }
    free(got);
}

/* Symbolic equality: Simplify[(a) - (b)] must be 0. */
static void same(const char* a, const char* b) {
    checks++;
    size_t n = strlen(a) + strlen(b) + 32;
    char* buf = malloc(n);
    snprintf(buf, n, "Simplify[(%s) - (%s)]", a, b);
    char* s = eval_str(buf);
    if (strcmp(s, "0") != 0) {
        fprintf(stderr, "FAIL (oracle): %s  vs  %s\n  Simplify[diff] = %s\n",
                a, b, s);
        failures++;
    }
    free(s); free(buf);
}

/* Numeric closeness inside the language: Abs[N[a-b, 20]] < 10^-tolexp must be
 * True.  Used to cross-check a closed form against a long partial sum. */
static void nclose(const char* a, const char* b, int tolexp) {
    checks++;
    size_t n = strlen(a) + strlen(b) + 64;
    char* buf = malloc(n);
    snprintf(buf, n, "TrueQ[Abs[N[(%s) - (%s), 25]] < 10^(-%d)]", a, b, tolexp);
    char* s = eval_str(buf);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL (numeric): %s  vs  %s  (tol 1e-%d) -> %s\n",
                a, b, tolexp, s);
        failures++;
    }
    free(s); free(buf);
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- Phase A: rational poles -> Zeta ---- */
    check("Sum[1/i^2, {i, 1, Infinity}]", "1/6 Pi^2");
    check("Sum[1/i^3, {i, 1, Infinity}]", "Zeta[3]");
    check("Sum[1/i^4, {i, 1, Infinity}]", "1/90 Pi^4");
    check("Sum[1/(i (i + 2)), {i, 1, Infinity}]", "3/4");
    /* shifted lower bound: Zeta[2, 1] and Zeta[2, 2] = Pi^2/6 - 1 */
    check("Sum[1/(i - 1)^2, {i, 2, Infinity}]", "1/6 Pi^2");
    check("Sum[1/i^2, {i, 2, Infinity}]", "-1 + 1/6 Pi^2");

    /* divergent / non-applicable inputs stay held */
    check("Sum[i, {i, 1, Infinity}]", "Sum[i, {i, 1, Infinity}]");
    check("Sum[1/i, {i, 1, Infinity}]", "Sum[1/i, {i, 1, Infinity}]");

    /* ---- Phase B: complex / radical poles -> PolyGamma ---- */
    check("Sum[1/(i (i^2 + 1)), {i, 1, Infinity}]",
          "1/2 (2 EulerGamma + PolyGamma[0, 1 - I] + PolyGamma[0, 1 + I])");
    /* four-PolyGamma form: assert numeric value against a long partial sum */
    nclose("Sum[1/((i^2 + 3 i + 1) (i^2 + 1)), {i, 1, Infinity}]",
           "Sum[1/((i^2 + 3 i + 1) (i^2 + 1)), {i, 1, 3000}]", 4);

    /* ---- Phase C: conjugate quadratic -> Coth ---- */
    check("Sum[1/(i^2 (i^2 + 1)), {i, 1, Infinity}]",
          "1/6 (3 + Pi^2 - 3 Pi Coth[Pi])");
    check("Sum[1/(i^2 + 1), {i, 1, Infinity}]", "1/2 (-1 + Pi Coth[Pi])");
    /* beta = 2 (Coth[2 Pi]) */
    same("Sum[1/(i^2 (i^2 + 4)), {i, 1, Infinity}]",
         "1/96 (3 + 4 Pi^2 - 6 Pi Coth[2 Pi])");
    /* shifted t = 2 uses the finite correction; exact value = full - 1/2 */
    same("Sum[1/(i^2 + 1), {i, 2, Infinity}]", "1/2 (-2 + Pi Coth[Pi])");

    /* ---- numeric cross-checks of the closed forms vs partial sums ---- */
    nclose("Sum[1/i^2, {i, 1, Infinity}]", "Sum[1/i^2, {i, 1, 3000}]", 3);
    nclose("Sum[1/(i^2 (i^2 + 1)), {i, 1, Infinity}]",
           "Sum[1/(i^2 (i^2 + 1)), {i, 1, 3000}]", 8);
    nclose("Sum[1/(i (i^2 + 1)), {i, 1, Infinity}]",
           "Sum[1/(i (i^2 + 1)), {i, 1, 3000}]", 4);

    if (failures) {
        fprintf(stderr, "\n%d/%d sum_rational checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum_rational tests passed!\n", checks);
    return 0;
}
