/*
 * test_sum_euler.c -- Sum`Euler: infinite linear Euler sums.
 *
 * check() compares the InputForm string of the evaluated input against an
 * expected string; same() cross-checks a symbolic closed form via
 * Simplify[a - b] == 0; nclose() checks a numeric bound.
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

/* Simplify[(a) - (b)] must be 0 (robust to differing printed forms). */
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

    /* ---- Order 1 (Euler's formula), any outer power q >= 2 ---- */
    check("Sum[HarmonicNumber[k]/k^2, {k, 1, Infinity}]", "2 Zeta[3]");
    same("Sum[HarmonicNumber[k]/k^3, {k, 1, Infinity}]", "Pi^4/72");
    /* Sum H_k/k^4 = 3 Zeta[5] - Zeta[2] Zeta[3] */
    same("Sum[HarmonicNumber[k]/k^4, {k, 1, Infinity}]", "3 Zeta[5] - Zeta[2] Zeta[3]");
    /* constant coefficient is carried through */
    check("Sum[3 HarmonicNumber[k]/k^2, {k, 1, Infinity}]", "6 Zeta[3]");

    /* ---- Diagonal p = q ---- */
    /* Sum H_k^(2)/k^2 = (Zeta[2]^2 + Zeta[4])/2 = 7/4 Zeta[4] = 7 Pi^4/360 */
    same("Sum[HarmonicNumber[k, 2]/k^2, {k, 1, Infinity}]", "7/4 Zeta[4]");
    same("Sum[HarmonicNumber[k, 3]/k^3, {k, 1, Infinity}]",
         "(Zeta[3]^2 + Zeta[6])/2");

    /* ---- Non-diagonal linear, odd weight (double-zeta reduction) ---- */
    /* Sum H_k^(2)/k^3 = 3 Zeta[2] Zeta[3] - 9/2 Zeta[5]  (outer q = 3 odd) */
    same("Sum[HarmonicNumber[k, 2]/k^3, {k, 1, Infinity}]",
         "3 Zeta[2] Zeta[3] - 9/2 Zeta[5]");
    /* Sum H_k^(3)/k^2 = 11/2 Zeta[5] - 2 Zeta[2] Zeta[3]  (outer q = 2 even: reflection) */
    same("Sum[HarmonicNumber[k, 3]/k^2, {k, 1, Infinity}]",
         "11/2 Zeta[5] - 2 Zeta[2] Zeta[3]");
    /* Sum H_k^(4)/k^3 (weight 7 odd) numeric cross-check (converges ~1/k^3, so
     * 800 terms give ~1e-6; the O(N^2) exact-harmonic partial sum keeps N small). */
    check("Chop[N[Sum[HarmonicNumber[k, 4]/k^3, {k, 1, Infinity}] "
          "- Sum[HarmonicNumber[k, 4]/k^3, {k, 1, 800}], 12], 10^-4]", "0");

    /* ---- Quadratic sums H_k^2/k^q, q = 2..5 ---- */
    same("Sum[HarmonicNumber[k]^2/k^2, {k, 1, Infinity}]", "17/4 Zeta[4]");
    same("Sum[HarmonicNumber[k]^2/k^3, {k, 1, Infinity}]",
         "7/2 Zeta[5] - Zeta[2] Zeta[3]");
    same("Sum[HarmonicNumber[k]^2/k^4, {k, 1, Infinity}]",
         "97/24 Zeta[6] - 2 Zeta[3]^2");
    same("Sum[HarmonicNumber[k]^2/k^5, {k, 1, Infinity}]",
         "6 Zeta[7] - Zeta[2] Zeta[5] - 5/2 Zeta[3] Zeta[4]");

    /* ---- Divergent / out-of-scope inputs stay unevaluated ---- */
    /* q = 1 diverges */
    check("Sum[HarmonicNumber[k]/k, {k, 1, Infinity}]",
          "Sum[HarmonicNumber[k]/k, {k, 1, Infinity}]");
    /* even weight non-diagonal: no zeta reduction */
    check("Sum[HarmonicNumber[k, 2]/k^4, {k, 1, Infinity}]",
          "Sum[HarmonicNumber[k, 2]/k^4, {k, 1, Infinity}]");
    /* quadratic weight 8 (q = 6): irreducible MZV, stays held */
    check("Sum[HarmonicNumber[k]^2/k^6, {k, 1, Infinity}]",
          "Sum[HarmonicNumber[k]^2/k^6, {k, 1, Infinity}]");

    if (failures) {
        fprintf(stderr, "\n%d/%d sum_euler checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum_euler tests passed!\n", checks);
    return 0;
}
