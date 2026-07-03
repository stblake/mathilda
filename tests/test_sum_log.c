/*
 * test_sum_log.c -- extensive unit tests for the log-weighted infinite Sum
 * families:
 *
 *   Sum`LogZeta     Sum[c Log[k]/k^s] = -c Zeta'[s]  (elementary only for s==2,
 *                   via the Glaisher bridge);
 *   Sum`LogRational Sum[rational + Log[rational]] convergent combinations, in
 *                   EulerGamma / Log via matched digamma/LogGamma asymptotics.
 *
 * same()  Simplify[(a)-(b)] == 0.
 * held()  asserts the sum stayed unevaluated (head still Sum).
 * nclose() numeric cross-check against a long partial sum.
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
        fprintf(stderr, "FAIL (same): %s  vs  %s\n  Simplify[diff] = %s\n", a, b, s);
        failures++;
    }
    free(s); free(buf);
}

static void held(const char* input) {
    checks++;
    char* got = eval_str(input);
    if (strncmp(got, "Sum[", 4) != 0) {
        fprintf(stderr, "FAIL (held): %s\n  expected held, got: %s\n", input, got);
        failures++;
    }
    free(got);
}

static void nclose(const char* closed, const char* finite, const char* tol) {
    checks++;
    char buf[1400];
    snprintf(buf, sizeof buf, "Abs[N[(%s) - (%s)]] < %s", closed, finite, tol);
    char* s = eval_str(buf);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL (nclose): %s  vs  %s -> %s\n", closed, finite, s);
        failures++;
    }
    free(s);
}

/* The Glaisher closed form of -Zeta'[2]. */
#define NEGZ2 "(Pi^2/6)(12 Log[Glaisher] - EulerGamma - Log[2 Pi])"

int main(void) {
    symtab_init();
    core_init();

    /* ================= Sum`LogZeta (s == 2 only) ================= */
    same("Sum[Log[k]/k^2, {k, 1, Infinity}]", NEGZ2);
    same("Sum[5 Log[k]/k^2, {k, 1, Infinity}]", "5 (" NEGZ2 ")");
    same("Sum[Log[k]/(2 k^2), {k, 1, Infinity}]", "(1/2)(" NEGZ2 ")");
    nclose(NEGZ2, "Sum[Log[k]/k^2, {k, 1, 3000}]", "10^-2");
    /* s != 2 has no elementary Zeta' form -> held */
    held("Sum[Log[k]/k^3, {k, 1, Infinity}]");
    held("Sum[Log[k]/k^(3/2), {k, 1, Infinity}]");
    /* s <= 1 diverges -> held */
    held("Sum[Log[k]/k, {k, 1, Infinity}]");

    /* ================= Sum`LogRational ================= */
    check("Sum[1/k + Log[(k-1)/k], {k, 2, Infinity}]", "-1 + EulerGamma");
    check("Sum[1/k - Log[(k+1)/k], {k, 1, Infinity}]", "EulerGamma");
    same("Sum[2/k + 2 Log[(k-1)/k], {k, 2, Infinity}]", "2 (EulerGamma - 1)");
    /* pure convergent telescoping (harmonic pair cancels, no gamma constant) */
    same("Sum[1/k - 1/(k+1), {k, 1, Infinity}]", "1");
    nclose("EulerGamma - 1", "Sum[1/k + Log[(k-1)/k], {k, 2, 4000}]", "10^-3");
    nclose("EulerGamma", "Sum[1/k - Log[(k+1)/k], {k, 1, 4000}]", "10^-3");
    /* divergent: harmonic without a cancelling log, or a lone divergent log */
    held("Sum[1/k + Log[k], {k, 1, Infinity}]");
    held("Sum[Log[(k-1)/k], {k, 2, Infinity}]");
    held("Sum[1/k + Log[(k-1)/k^2], {k, 2, Infinity}]"); /* counts unequal -> diverge */

    if (failures == 0)
        printf("All %d sum-log tests passed!\n", checks);
    else
        fprintf(stderr, "%d/%d sum-log checks FAILED\n", failures, checks);
    return failures ? 1 : 0;
}
