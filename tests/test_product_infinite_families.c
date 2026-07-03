/*
 * test_product_infinite_families.c -- extensive unit tests for the infinite
 * Product[] families added on top of Stages 0-3:
 *
 *   Product`Geometric (infinite base^summable-exponent),
 *   Product`Viete       (cosine double-angle telescoping),
 *   Product`Cantor      (double-exponential telescoping),
 *   Product`EulerPrime  (Euler products / Dirichlet chi_4),
 *   Product`LogSum      (Exp/log-sum bridge),
 *   Product`RationalInfinite (Gamma canonical form, complex-conjugate roots).
 *
 * check()  exact InputForm string (clean rational / symbolic outputs).
 * same()   Simplify[(a)-(b)] == 0 (robust to Csch vs 1/Sinh, term order, ...).
 * held()   asserts the product stayed unevaluated (head still Product) --
 *          the correct behaviour for divergent / unrecognised inputs.
 * nclose() numeric cross-check of a closed form against a long finite product.
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

/* Assert the input stays unevaluated (head still Product). */
static void held(const char* input) {
    checks++;
    char* got = eval_str(input);
    if (strncmp(got, "Product[", 8) != 0) {
        fprintf(stderr, "FAIL (held): %s\n  expected held, got: %s\n", input, got);
        failures++;
    }
    free(got);
}

/* Numeric cross-check: |N[closed - finite]| < tol.  finite is a long partial
 * product; tol is loose because product tails converge slowly (~1/N or ~log). */
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

int main(void) {
    symtab_init();
    core_init();

    /* ================= Product`Geometric (infinite) ================= */
    check("Product[2^(k/2^k), {k, 1, Infinity}]", "4");        /* 2^Sum[k/2^k]=2^2 */
    check("Product[3^(k/2^k), {k, 1, Infinity}]", "9");        /* 3^2 */
    check("Product[2^(1/2^k), {k, 1, Infinity}]", "2");        /* 2^Sum[1/2^k]=2^1 */
    same("Product[4^(k/3^k), {k, 1, Infinity}]", "4^(3/4)");   /* Sum[k/3^k]=3/4 */
    same("Product[E^(k/2^k), {k, 1, Infinity}]", "E^2");
    /* divergent exponent series -> held */
    held("Product[2^k, {k, 1, Infinity}]");
    held("Product[2^(1/k), {k, 1, Infinity}]");                /* Sum[1/k] diverges */

    /* ================= Product`Viete ================= */
    check("Product[Cos[Pi/2^(k+1)], {k, 1, Infinity}]", "2/Pi");
    same("Product[Cos[x/2^k], {k, 1, Infinity}]", "Sin[x]/x");
    same("Product[Cos[Pi/2^k], {k, 2, Infinity}]", "2/Pi");
    same("Product[Cos[1/2^k], {k, 1, Infinity}]", "Sin[1]");   /* Sin[1]/1 */
    nclose("2/Pi", "Product[Cos[Pi/2^(k+1)], {k, 1, 25}]", "10^-6");
    /* not halving / not a cosine -> held */
    held("Product[Cos[k], {k, 1, Infinity}]");
    held("Product[Cos[Pi/3^k], {k, 1, Infinity}]");
    held("Product[Cosh[Pi/2^k], {k, 1, Infinity}]");

    /* ================= Product`Cantor ================= */
    check("Product[1 + (1/3)^(2^k), {k, 0, Infinity}]", "3/2");
    check("Product[1 + (1/2)^(2^k), {k, 0, Infinity}]", "2");
    check("Product[1 + (1/4)^(2^k), {k, 0, Infinity}]", "4/3");
    check("Product[1 + (1/3)^(2^k), {k, 1, Infinity}]", "9/8"); /* M=2 -> 1/(1-1/9) */
    nclose("3/2", "Product[1 + (1/3)^(2^k), {k, 0, 12}]", "10^-6");
    /* |base| >= 1 -> held (would diverge / not the Cantor identity) */
    held("Product[1 + 2^(2^k), {k, 0, Infinity}]");
    held("Product[1 + (1/3)^(3^k), {k, 0, Infinity}]");        /* tripling, not doubling */

    /* ================= Product`EulerPrime ================= */
    same("Product[1/(1 - 1/Prime[k]^2), {k, 1, Infinity}]", "Pi^2/6");
    same("Product[1/(1 - 1/Prime[k]^4), {k, 1, Infinity}]", "Pi^4/90");
    same("Product[1/(1 - 1/Prime[k]^6), {k, 1, Infinity}]", "Pi^6/945");
    check("Product[1/(1 - 1/Prime[k]^3), {k, 1, Infinity}]", "Zeta[3]");
    same("Product[1/(1 - (-1)^((Prime[k]-1)/2)/Prime[k]), {k, 2, Infinity}]", "Pi/4");
    nclose("Pi^2/6", "Product[1/(1 - 1/Prime[k]^2), {k, 1, 500}]", "10^-2");
    nclose("Pi/4", "Product[1/(1 - (-1)^((Prime[k]-1)/2)/Prime[k]), {k, 2, 2000}]", "10^-2");
    /* divergent (s <= 1) -> held */
    held("Product[1/(1 - 1/Prime[k]), {k, 1, Infinity}]");
    /* a non-prime rational body is not EulerPrime's -- but still evaluates via the
     * ordinary rational/limit route (reciprocal of prod(1-1/k^2) = 1/2). */
    check("Product[1/(1 - 1/k^2), {k, 2, Infinity}]", "2");

    /* ================= Product`LogSum ================= */
    check("Product[Exp[(-1)^k/k], {k, 1, Infinity}]", "1/2");
    same("Product[Exp[1/k^2], {k, 1, Infinity}]", "E^(Pi^2/6)");
    same("Product[Exp[1/k^3], {k, 1, Infinity}]", "E^Zeta[3]");
    same("Product[Exp[(-1)^k/k^2], {k, 1, Infinity}]", "E^(-Pi^2/12)"); /* Sum[(-1)^k/k^2]=-Pi^2/12 */
    /* via Sum`LogZeta: k^(1/k^2) -> Exp[-Zeta'[2]] in the explicit Glaisher form
     * (the system has no symbolic Zeta', so we compare against that form). */
    same("Product[k^(1/k^2), {k, 1, Infinity}]",
         "Exp[(Pi^2/6)(12 Log[Glaisher] - EulerGamma - Log[2 Pi])]");
    nclose("Exp[(Pi^2/6)(12 Log[Glaisher] - EulerGamma - Log[2 Pi])]",
           "Product[k^(1/k^2), {k, 1, 4000}]", "10^-2");
    /* via Sum`LogRational */
    same("Product[E^(1/k) (1 - 1/k), {k, 2, Infinity}]", "E^(EulerGamma - 1)");
    nclose("E^(EulerGamma - 1)", "Product[E^(1/k) (1 - 1/k), {k, 2, 4000}]", "10^-3");
    /* divergent log-sum -> held */
    held("Product[Exp[1/k], {k, 1, Infinity}]");

    /* ================= Product`RationalInfinite (Gamma) ================= */
    same("Product[(k^2-1)/(k^2+1), {k, 2, Infinity}]", "Pi/Sinh[Pi]");
    check("Product[(k^3-1)/(k^3+1), {k, 2, Infinity}]", "2/3");
    same("Product[1 + 1/(2k-1)^2, {k, 1, Infinity}]", "Cosh[Pi/2]");
    same("Product[(k^2-4)/(k^2+4), {k, 3, Infinity}]", "(10 Pi)/(3 Sinh[2 Pi])");
    nclose("Pi/Sinh[Pi]", "Product[(k^2-1)/(k^2+1), {k, 2, 5000}]", "10^-2");
    nclose("2/3", "Product[(k^3-1)/(k^3+1), {k, 2, 5000}]", "10^-2");
    nclose("Cosh[Pi/2]", "Product[1 + 1/(2k-1)^2, {k, 1, 5000}]", "10^-2");
    /* irreducible quartic (roots have irrational real part 1/Sqrt[2]) has no
     * elementary sinh/cosh collapse -> RationalInfinite bails and it stays held. */
    held("Product[(k^4-1)/(k^4+1), {k, 2, Infinity}]");
    /* real-root convergent products keep their existing (Weierstrass/limit) route */
    same("Product[1 + 1/k^2, {k, 1, Infinity}]", "Sinh[Pi]/Pi");
    check("Product[1 - 1/k^2, {k, 2, Infinity}]", "1/2");
    /* divergent rational product -> held */
    held("Product[(k+1)/k, {k, 1, Infinity}]");

    if (failures == 0)
        printf("All %d product-infinite-family tests passed!\n", checks);
    else
        fprintf(stderr, "%d/%d product-infinite-family checks FAILED\n", failures, checks);
    return failures ? 1 : 0;
}
