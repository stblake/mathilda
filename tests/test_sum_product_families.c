/*
 * test_sum_product_families.c -- new infinite Sum/Product closed-form families:
 *   Sum`Dirichlet, Sum`ZetaSeries, Sum`LogZeta (general s),
 *   Sum`EulerNonlinear, Product`QProduct (infinite), Product`BesselZero.
 *
 * check()  compares the InputForm string of the result to an expected string.
 * same()   cross-checks mathematical equality via Simplify[a - b] == 0 (robust
 *          to equivalent-but-differently-printed forms, e.g. Zeta[6] vs Pi^6).
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

    /* ---- Sum`Dirichlet: arithmetic-function Dirichlet series ---- */
    check("Sum[MoebiusMu[k]/k^3, {k, 1, Infinity}]", "1/Zeta[3]");
    check("Sum[MoebiusMu[k]/k^s, {k, 1, Infinity}]", "1/Zeta[s]");
    check("Sum[LiouvilleLambda[k]/k^2, {k, 1, Infinity}]", "1/15 Pi^2");
    check("Sum[DivisorSigma[0, k]/k^s, {k, 1, Infinity}]", "Zeta[s]^2");
    check("Sum[EulerPhi[k]/k^s, {k, 1, Infinity}]", "Zeta[-1 + s]/Zeta[s]");
    check("Sum[MoebiusMu[k]^2/k^s, {k, 1, Infinity}]", "Zeta[s]/Zeta[2 s]");
    same("Sum[DivisorSigma[0, k]^2/k^3, {k, 1, Infinity}]", "945 Zeta[3]^4/Pi^6");
    /* not a Dirichlet series -> falls through, stays a plain zeta/rational sum */
    check("Sum[1/k^s, {k, 1, Infinity}]", "Sum[1/k^s, {k, 1, Infinity}]");

    /* ---- Sum`ZetaSeries: order-swap of Zeta[a k + b] series ---- */
    check("Sum[Zeta[k] - 1, {k, 2, Infinity}]", "1");
    check("Sum[(-1)^k (Zeta[k] - 1), {k, 2, Infinity}]", "1/2");
    check("Sum[Zeta[2 k] - 1, {k, 1, Infinity}]", "3/4");
    check("Product[Exp[Zeta[k] - 1], {k, 2, Infinity}]", "E");

    /* ---- Sum`LogZeta: general integer s -> -c Zeta'[s] ---- */
    check("Sum[Log[k]/k^3, {k, 1, Infinity}]", "-Derivative[1][Zeta][3]");
    check("Sum[2 Log[k]/k^4, {k, 1, Infinity}]", "-2 Derivative[1][Zeta][4]");
    check("Product[k^(1/k^3), {k, 1, Infinity}]", "E^(-Derivative[1][Zeta][3])");

    /* ---- Sum`EulerNonlinear: nonlinear Euler sums (weight <= 7) ---- */
    check("Sum[PolyGamma[1, k]^2, {k, 1, Infinity}]", "3 Zeta[3]");
    same("Sum[HarmonicNumber[k] HarmonicNumber[k, 2]/k^3, {k, 1, Infinity}]",
         "-(101/48) Zeta[6] + (5/2) Zeta[3]^2");
    same("Sum[HarmonicNumber[k]^3/k^4, {k, 1, Infinity}]",
         "(231/16) Zeta[7] + 2 Zeta[2] Zeta[5] - (51/4) Zeta[3] Zeta[4]");
    /* regression: the narrow Sum`Euler cases are untouched */
    same("Sum[HarmonicNumber[k]^2/k^2, {k, 1, Infinity}]", "17/4 Zeta[4]");
    same("Sum[HarmonicNumber[k]/k^2, {k, 1, Infinity}]", "2 Zeta[3]");
    /* weight 8 is irreducible -> falls through (never a wrong value) */
    check("Sum[HarmonicNumber[k]^2/k^6, {k, 1, Infinity}]",
          "Sum[HarmonicNumber[k]^2/k^6, {k, 1, Infinity}]");

    /* ---- Product`QProduct: infinite q-products -> QPochhammer ---- */
    check("Product[1 - q^k, {k, 1, Infinity}]", "QPochhammer[q, q]");
    check("Product[1/(1 - x^k), {k, 1, Infinity}]", "1/QPochhammer[x, x]");
    check("Product[1 - q^(2 k), {k, 1, Infinity}]", "QPochhammer[q^2, q^2]");
    check("Product[1 - a q^k, {k, 0, Infinity}]", "QPochhammer[a, q]");
    /* finite form regression */
    check("Product[1 - q^k, {k, 1, n}]", "QPochhammer[q, q, n]");

    /* ---- Product`BesselZero: Bessel Hadamard product ---- */
    check("Product[1 - x^2/BesselJZero[0, k]^2, {k, 1, Infinity}]", "BesselJ[0, x]");
    same("Product[1 - x^2/BesselJZero[1, k]^2, {k, 1, Infinity}]", "2 BesselJ[1, x]/x");
    check("BesselJZero[0, 3]", "BesselJZero[0, 3]");

    if (failures) {
        fprintf(stderr, "\n%d/%d checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum/product family checks passed.\n", checks);
    return 0;
}
