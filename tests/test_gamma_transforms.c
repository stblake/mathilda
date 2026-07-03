/*
 * test_gamma_transforms.c -- unit tests for the Gamma-family transformation
 * rules in src/internal/simp/transforms/gamma.m, exercised through FullSimplify.
 *
 * Covers the recurrence / raising rules and the reflection + conjugate-pair
 * (Re == 1 -> Sinh, Re == 1/2 -> Cosh) reductions added for the infinite
 * rational-product work, plus the guards that keep the reflection from emitting
 * a spurious ComplexInfinity at integer arguments.
 *
 * check() exact InputForm; same() Simplify[(a)-(b)]==0; nsame() numeric equality.
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

/* Numeric equality |N[a-b]| < 10^-9. */
static void nsame(const char* a, const char* b) {
    checks++;
    char buf[512];
    snprintf(buf, sizeof buf, "Abs[N[(%s) - (%s)]] < 10^-9", a, b);
    char* s = eval_str(buf);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL (nsame): %s  vs  %s -> %s\n", a, b, s);
        failures++;
    }
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    /* -------- recurrence / raising (pre-existing, must still hold) -------- */
    check("FullSimplify[Gamma[x+1]/Gamma[x]]", "x");
    check("FullSimplify[Gamma[x+1] - x Gamma[x]]", "0");

    /* -------- reflection pair: Gamma[z] Gamma[1-z] = Pi/Sin[Pi z] -------- */
    same("FullSimplify[Gamma[1/4] Gamma[3/4]]", "Pi Sqrt[2]");
    same("FullSimplify[Gamma[1/3] Gamma[2/3]]", "2 Pi/Sqrt[3]");
    same("FullSimplify[Gamma[1/6] Gamma[5/6]]", "2 Pi");
    same("FullSimplify[Gamma[1/2]^2]", "Pi");
    /* symbolic order: reflection fires (n assumed non-integer) */
    same("FullSimplify[Gamma[n] Gamma[1-n]]", "Pi Csc[Pi n]");
    /* reflection numeric sanity */
    nsame("FullSimplify[Gamma[1/5] Gamma[4/5]]", "Pi/Sin[Pi/5]");

    /* -------- conjugate pair, Re == 1: Gamma[1+I b]Gamma[1-I b] = Pi b/Sinh[Pi b] -- */
    same("FullSimplify[Gamma[1+I] Gamma[1-I]]", "Pi/Sinh[Pi]");
    same("FullSimplify[Gamma[1+2 I] Gamma[1-2 I]]", "2 Pi/Sinh[2 Pi]");
    nsame("FullSimplify[Gamma[1+I/2] Gamma[1-I/2]]", "(Pi/2)/Sinh[Pi/2]");

    /* -------- conjugate pair, Re == 1/2: Gamma[1/2+I b]Gamma[1/2-I b] = Pi/Cosh[Pi b] -- */
    same("FullSimplify[Gamma[1/2+I] Gamma[1/2-I]]", "Pi/Cosh[Pi]");
    same("FullSimplify[Gamma[1/2+2 I] Gamma[1/2-2 I]]", "Pi/Cosh[2 Pi]");
    nsame("FullSimplify[Gamma[1/2+3 I] Gamma[1/2-3 I]]", "Pi/Cosh[3 Pi]");

    /* -------- integer-argument guard: no spurious ComplexInfinity -------- */
    check("FullSimplify[Gamma[2] Gamma[3]]", "2");     /* 1 * 2, reflection must NOT fire */
    check("FullSimplify[Gamma[5]]", "24");
    check("FullSimplify[Gamma[1] Gamma[1]]", "1");

    /* -------- non-pairs must be left alone (no false reduction) -------- */
    check("FullSimplify[Gamma[1/4]^2]", "Gamma[1/4]^2");   /* z+w != 1 */
    /* a lone gamma is not touched by the pair rules */
    check("FullSimplify[Gamma[1/4]]", "Gamma[1/4]");

    if (failures == 0)
        printf("All %d gamma-transform tests passed!\n", checks);
    else
        fprintf(stderr, "%d/%d gamma-transform checks FAILED\n", failures, checks);
    return failures ? 1 : 0;
}
