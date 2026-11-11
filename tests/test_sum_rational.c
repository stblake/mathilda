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

    /* ---- Finite ranges: polygamma closed form, cost independent of width ----
     * A range of more than 10^6 terms tries the closed form before enumerating;
     * Sum[1/(k + 3/7), {k, 12, 3*10^6}] used to enumerate 3*10^6 big rationals
     * and never finish.  Forms are the reference's. */
    check("Sum[1/(k + 3/7), {k, 12, 3000000}]",
          "-PolyGamma[0, 87/7] + PolyGamma[0, 21000010/7]");
    check("Sum[1/k, {k, 1, 10^7}]", "EulerGamma + PolyGamma[0, 10000001]");
    check("Sum[1/k^2, {k, 1, 10^8}]", "1/6 Pi^2 - PolyGamma[1, 100000001]");
    check("Sum[1/k^3, {k, 1, 10^7}]",
          "1/2 (-PolyGamma[2, 1] + PolyGamma[2, 10000001])");
    check("Sum[1/(2 k + 1), {k, 0, 10^7}]",
          "1/2 (-PolyGamma[0, 1/2] + PolyGamma[0, 20000003/2])");
    check("Sum[1/(k + a), {k, 1, 2*10^6}]",
          "-PolyGamma[0, 1 + a] + PolyGamma[0, 2000001 + a]");
    /* telescoping bodies still go to Gosper first */
    check("Sum[1/(k (k + 1)), {k, 1, 10^7}]", "10000000/10000001");
    /* irreducible cubic: RootSum over the Root poles, residue reduced mod q */
    check("Sum[1/(k^3 + k + 1), {k, 1, 10^7}]",
          "-RootSum[1 + #1 + #1^3 &, (4/31 - 9/31 #1 + 6/31 #1^2) PolyGamma[0, 1 - #1] &]"
          " + RootSum[1 + #1 + #1^3 &, (4/31 - 9/31 #1 + 6/31 #1^2) PolyGamma[0, 10000001 - #1] &]");
    /* the 10^7 value against the infinite sum: the tail is 1/(2 10^14) + O(10^-21) */
    nclose("Sum[1/(k^3 + k + 1), {k, 1, 10^7}] + 1/(2 10^14)",
           "Sum[1/(k^3 + k + 1), {k, 1, Infinity}]", 20);

    /* The stage itself on short ranges, against term-by-term enumeration:
     * exact where the polygammas evaluate (integer arguments), else to 10^-20. */
    check("Sum`Rational[1/k^2, k, 1, 300] == HarmonicNumber[300, 2]", "True");
    check("Sum`Rational[(k^2 + 1)/(k^2 + 3 k), k, 1, 300] == "
          "Sum[(j^2 + 1)/(j^2 + 3 j), {j, 1, 300}]", "True");
    check("Sum`Rational[1/(2 k + 1)^2, k, 0, 300] == Sum[1/(2 j + 1)^2, {j, 0, 300}]", "True");
    nclose("Sum`Rational[1/(k + 3/7), k, 12, 300]", "Sum[1/(j + 3/7), {j, 12, 300}]", 20);
    nclose("Sum`Rational[1/(k^2 + 1), k, 1, 300]", "Sum[1/(j^2 + 1), {j, 1, 300}]", 20);
    nclose("Sum`Rational[1/(k^2 + 1)^2, k, 1, 300]", "Sum[1/(j^2 + 1)^2, {j, 1, 300}]", 20);
    nclose("Sum`Rational[k/(k^2 + 1), k, 1, 300]", "Sum[j/(j^2 + 1), {j, 1, 300}]", 20);
    nclose("Sum`Rational[1/(k^2 - 2), k, 2, 300]", "Sum[1/(j^2 - 2), {j, 2, 300}]", 20);
    nclose("Sum`Rational[1/(k^3 + k + 1), k, 5, 300]", "Sum[1/(j^3 + j + 1), {j, 5, 300}]", 20);
    nclose("Sum`Rational[1/(k^5 + k + 1), k, 1, 200]", "Sum[1/(j^5 + j + 1), {j, 1, 200}]", 20);
    nclose("Sum`Rational[1/((k^3 + k + 1) (k^2 + 2)), k, 1, 200]",
           "Sum[1/((j^3 + j + 1) (j^2 + 2)), {j, 1, 200}]", 20);
    nclose("Sum`Rational[(3 k^4 - k + 2)/((k + 1/2)^2 (k^2 + k + 1)), k, -5, 300]",
           "Sum[(3 j^4 - j + 2)/((j + 1/2)^2 (j^2 + j + 1)), {j, -5, 300}]", 18);
    nclose("Sum`Rational[k^3/(k + 1/3)^3, k, -20, 200]", "Sum[j^3/(j + 1/3)^3, {j, -20, 200}]", 18);
    nclose("Sum`Rational[1/(k + 1/2), k, -300, -2]", "Sum[1/(j + 1/2), {j, -300, -2}]", 20);
    /* a pole on an index of the range has no finite value: the stage declines */
    check("Sum`Rational[1/(k - 5), k, 1, 100]", "Sum`Rational[1/(-5 + k), k, 1, 100]");
    /* symbolic finite bounds are not handled by this stage */
    check("Sum`Rational[1/(k + 3/7), k, 12, n]", "Sum`Rational[1/(3/7 + k), k, 12, n]");

    /* ---- Non-unit step: reindexed to a unit-step range first ---- */
    check("Sum[1/(k + 3/7), {k, 12, 3000000, 2}]",
          "1/2 (-PolyGamma[0, 87/14] + PolyGamma[0, 21000017/14])");
    check("Sum[1/(k + 3/7), {k, 12, 3000000, 3/2}]",
          "2/3 (-PolyGamma[0, 58/7] + PolyGamma[0, 14000009/7])");
    check("Sum[1/k, {k, 1, 10^7, 3}]",
          "1/3 (-PolyGamma[0, 1/3] + PolyGamma[0, 10000003/3])");
    check("Sum[k^2, {k, 1, 10^9, 7}]", "47619047404761904500000000");
    check("Sum[k^2, {k, 10^9, 1, -7}] == Sum[k^2, {k, 6, 10^9, 7}]", "True");

    /* ---- Inexact terms keep the numeric path (no exact closed form) ---- */
    check("Head[Sum[1.0/k, {k, 1, 10^7}]]", "Real");
    check("Head[Sum[N[1/k], {k, 1, 2*10^6}]]", "Real");
    nclose("Sum[1.0/k, {k, 1, 10^7}]", "EulerGamma + PolyGamma[0, 10000001]", 9);
    check("Head[Sum[1.0 k^2, {k, 1, 10^9}]]", "Real");

    /* ---- Short ranges still enumerate to exact rationals ---- */
    check("Sum[1/k, {k, 1, 10}]", "7381/2520");
    check("Sum[1/(k + 3/7), {k, 12, 3000}] == Sum[1/(j + 3/7), {j, 12, 3000}]", "True");
    check("Head[Sum[1/(k^3 + k + 1), {k, 1, 1000}]]", "Rational");
    check("Sum[x^i, {i, 0, 5}]", "1 + x + x^2 + x^3 + x^4 + x^5");

    if (failures) {
        fprintf(stderr, "\n%d/%d sum_rational checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum_rational tests passed!\n", checks);
    return 0;
}
