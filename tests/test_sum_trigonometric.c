/*
 * test_sum_trigonometric.c -- Sum`Trigonometric: infinite Fourier-type sums.
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

    /* ---- Headline elementary closed forms (s = 1) ---- */
    check("Sum[Sin[k]/k, {k, 1, Infinity}]", "1/2 (-1 + Pi)");
    same("Sum[Sin[k]/k, {k, 1, Infinity}]", "(Pi - 1)/2");
    same("Sum[Cos[k]/k, {k, 1, Infinity}]", "-Log[2 Sin[1/2]]");
    check("Sum[Sin[2 k]/k, {k, 1, Infinity}]", "1/2 (-2 + Pi)");
    same("Sum[Sin[3 k]/k, {k, 1, Infinity}]", "(Pi - 3)/2");

    /* ---- Linearisation of trig polynomials ---- */
    /* Sin[k]^2/k^2 = 1/2 (1 - Cos[2k])/k^2 -> Zeta[2]/2 - (1/2) Re PolyLog[2,e^{2i}]. */
    check("Sum[Sin[k]^2/k^2, {k, 1, Infinity}]",
          "1/12 (Pi^2 - 6 Re[PolyLog[2, E^(2*I)]])");

    /* ---- s >= 2 stays as Im/Re PolyLog (Wolfram parity) ---- */
    check("Sum[Sin[k]/k^2, {k, 1, Infinity}]", "Im[PolyLog[2, E^I]]");
    check("Sum[Cos[2 k]/k^2, {k, 1, Infinity}]", "Re[PolyLog[2, E^(2*I)]]");

    /* ---- lower limit != 1 : subtract the head ---- */
    same("Sum[Sin[k]/k, {k, 2, Infinity}]", "(Pi - 1)/2 - Sin[1]");

    /* ---- numeric cross-checks (conditionally convergent, lenient tol) ---- */
    check("Chop[N[Sum[Sin[k]/k, {k, 1, Infinity}] "
          "- Sum[Sin[k]/k, {k, 1, 6000}], 15], 10^-2]", "0");
    check("Chop[N[Sum[Cos[k]/k, {k, 1, Infinity}] "
          "- Sum[Cos[k]/k, {k, 1, 6000}], 15], 10^-2]", "0");

    /* ---- divergent / non-trig / symbolic inputs stay honest ---- */
    check("Sum[Sin[k]^2/k, {k, 1, Infinity}]",
          "Sum[Sin[k]^2/k, {k, 1, Infinity}]");        /* DC term at s=1 diverges */
    check("Sum[Sin[a k]/k, {k, 1, Infinity}]",
          "Im[-Log[1 - E^(I a)]]");                     /* symbolic coeff: no fabrication */
    check("Sum[1/k^2, {k, 1, Infinity}]", "1/6 Pi^2");  /* not intercepted */

    if (failures) {
        fprintf(stderr, "\n%d/%d sum_trigonometric checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum_trigonometric tests passed!\n", checks);
    return 0;
}
