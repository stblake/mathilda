/*
 * test_product_special.c -- prerequisite special functions for Product
 * Stages 5-6: Hyperfactorial, BarnesG, QPochhammer.  Exact integer values,
 * symbolic behaviour, and machine + MPFR N (including an N[x, 35] case each).
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

static int failures = 0, checks = 0;

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
        fprintf(stderr, "FAIL (oracle): %s vs %s => %s\n", a, b, s);
        failures++;
    }
    free(s); free(buf);
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- Hyperfactorial[n] = prod k^k ---- */
    check("Hyperfactorial[0]", "1");
    check("Hyperfactorial[1]", "1");
    check("Hyperfactorial[3]", "108");                 /* 1*4*27 */
    check("Hyperfactorial[4]", "27648");
    same("Hyperfactorial[5]", "Product[k^k, {k, 1, 5}]");
    check("Hyperfactorial[n]", "Hyperfactorial[n]");   /* symbolic stays */
    check("Hyperfactorial[5/2]", "Hyperfactorial[5/2]");
    check("Chop[N[Hyperfactorial[4], 35] - 27648]", "0");

    /* ---- BarnesG[z] ---- */
    check("BarnesG[1]", "1");
    check("BarnesG[2]", "1");
    check("BarnesG[3]", "1");
    check("BarnesG[4]", "2");
    check("BarnesG[5]", "12");
    check("BarnesG[6]", "288");
    check("BarnesG[0]", "0");
    check("BarnesG[-5]", "0");
    same("BarnesG[6]", "Product[Gamma[k], {k, 1, 5}]");  /* G(n+1)=prod_{1}^{n-1} k! */
    check("BarnesG[z]", "BarnesG[z]");
    check("Chop[N[BarnesG[6], 35] - 288]", "0");

    /* ---- QPochhammer[a, q, n] ---- */
    check("QPochhammer[a, q, 0]", "1");
    same("QPochhammer[a, q, 3]", "(1 - a) (1 - a q) (1 - a q^2)");
    check("QPochhammer[1, q, 4]", "0");                  /* k=0 factor (1-1)=0 */
    check("QPochhammer[a, q, n]", "QPochhammer[a, q, n]"); /* symbolic n stays */
    same("QPochhammer[q, q, 4]", "Product[1 - q^k, {k, 1, 4}]"); /* (q;q)_4 */
    /* numeric: 3-arg exact at rational args, MPFR N to 35 digits */
    same("QPochhammer[1/2, 1/3, 4]", "(1/2) (5/6) (17/18) (53/54)");
    check("Chop[N[QPochhammer[1/2, 1/3, 4], 35] - (1/2)(5/6)(17/18)(53/54)]", "0");
    /* infinite (a;q)_inf for machine-real args, |q|<1 */
    check("Chop[N[QPochhammer[0.5, 0.3] - (1 - 0.5) QPochhammer[0.5 0.3, 0.3]]]", "0");

    if (failures) {
        fprintf(stderr, "\n%d/%d special-function checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d product-special tests passed!\n", checks);
    return 0;
}
