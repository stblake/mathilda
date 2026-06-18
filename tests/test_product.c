/*
 * test_product.c -- Product (Stages 0-3): finite expansion, Telescoping,
 * Rational (Pochhammer/Gamma), Geometric, the Method polyalgorithm, and the
 * multiplicative edge cases (empty/reversed range -> 1).
 *
 * check()  compares the InputForm string against an expected canonical string
 *          (reserved for cases whose output form is canonical).
 * same()   cross-checks mathematical equality via Simplify[(a)-(b)] == 0
 *          (robust to Pochhammer/Gamma vs rational output forms).
 * oracle() substitutes the symbolic bound n -> {0..8} into a closed form and
 *          compares it to the direct finite expansion -- the multiplicative
 *          ground truth (n=0 exercises the empty-product == 1 convention).
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

/* Substitute n -> {0..8} into a symbolic closed form `closed_n` and compare to
 * the direct finite expansion produced by `finite_fmt` (one %d for the bound
 * that `n` controls).  Both sides evaluate to exact numbers. */
/* `n_start` is the smallest n the closed form covers: for a product whose
 * range starts at lo, the natural domain is n >= lo-1 (below that the range is
 * already empty and the closed form's continuation may diverge). */
static void oracle_from(const char* closed_n, const char* finite_fmt, int n_start) {
    for (int N = n_start; N <= 8; N++) {
        char a[768], b[768], buf[1700];
        snprintf(a, sizeof a, "(%s) /. n -> %d", closed_n, N);
        snprintf(b, sizeof b, finite_fmt, N);
        snprintf(buf, sizeof buf, "Simplify[(%s) - (%s)]", a, b);
        char* s = eval_str(buf);
        if (strcmp(s, "0") != 0) {
            fprintf(stderr, "FAIL (oracle n=%d): %s\n  vs %s\n  diff = %s\n",
                    N, closed_n, b, s);
            failures++;
        }
        checks++;
        free(s);
    }
}

/* Default oracle: lower bound 1, so n >= 0 (n=0 is the empty product == 1). */
static void oracle(const char* closed_n, const char* finite_fmt) {
    oracle_from(closed_n, finite_fmt, 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- Stage 0: finite explicit expansion ---- */
    check("Product[k, {k, 1, 5}]", "120");
    check("Product[i^2, {i, 1, 6}]", "518400");
    check("Product[2, {k, 1, 3}]", "8");
    check("Product[3, {k, 1, 4}]", "81");
    check("Product[1/i, {i, 1, 4}]", "1/24");
    check("Product[i, {i, 2, 10, 2}]", "3840");        /* 2*4*6*8*10 */
    check("Product[f[i], {i, 1, 4}]", "f[1] f[2] f[3] f[4]");
    check("Product[f[i], {i, 1, 4, 2}]", "f[1] f[3]"); /* step 2 */
    check("Product[f[i], {i, {a, b, c}}]", "f[a] f[b] f[c]");
    check("Product[k, {k, {2, 3, 5}}]", "30");
    /* empty / reversed range -> multiplicative identity 1 (NOT 0) */
    check("Product[k, {k, 1, 0}]", "1");
    check("Product[k, {k, 5, 1}]", "1");
    check("Product[i^2, {i, 3, 2}]", "1");
    /* multiple iterators: prod_i prod_j (i j) = (1*2)(2*4) = 16 */
    check("Product[i j, {i, 1, 2}, {j, 1, 2}]", "16");
    check("Product[f[i, j], {i, 1, 2}, {j, 1, 2}]",
          "f[1, 1] f[1, 2] f[2, 1] f[2, 2]");
    /* outer bound depends on inner variable */
    check("Product[f[i, j], {i, 1, 3}, {j, 1, i}]",
          "f[1, 1] f[2, 1] f[2, 2] f[3, 1] f[3, 2] f[3, 3]");
    /* attributes + degenerate calls stay unevaluated */
    check("MemberQ[Attributes[Product], HoldAll]", "True");
    check("MemberQ[Attributes[Product], Protected]", "True");
    check("Product[f]", "Product[f]");
    check("Product[k, {k, 1, n}, Method -> \"Bogus\"]",
          "Product[k, {k, 1, n}, Method -> \"Bogus\"]");
    /* closed-form-first on a wide unit-step numeric range (no per-term blowup) */
    check("Product[k, {k, 1, 50}]",
          "30414093201713378043612608166064768844377641568960512000000000000");

    /* ---- Stage 1: Telescoping (Gamma-free rational anti-quotient) ---- */
    check("Product[1 + 1/k, {k, 1, n}]", "1 + n");
    check("Product[1/k + 1, {k, 1, n}]", "1 + n");
    same("Product[k/(k + 1), {k, 1, n}]", "1/(1 + n)");
    oracle("Product[1 + 1/k, {k, 1, n}]", "Product[1 + 1/k, {k, 1, %d}]");
    oracle("Product[k/(k + 1), {k, 1, n}]", "Product[k/(k + 1), {k, 1, %d}]");
    oracle_from("Product[1 - 1/k^2, {k, 2, n}]", "Product[1 - 1/k^2, {k, 2, %d}]", 1);
    oracle("Product[(k + 1)/(k + 2), {k, 1, n}]", "Product[(k + 1)/(k + 2), {k, 1, %d}]");
    /* Method -> "Telescoping" forces the stage and yields Gamma-free output */
    check("FreeQ[Product[1 - 1/k^2, {k, 2, n}, Method -> \"Telescoping\"], Gamma]",
          "True");
    check("FreeQ[Product[1 - 1/k^2, {k, 2, n}, Method -> \"Telescoping\"], Pochhammer]",
          "True");

    /* ---- Stage 2: Rational (Pochhammer / Gamma / Factorial) ---- */
    check("Product[k, {k, 1, n}]", "Factorial[n]");
    check("Product[k^2, {k, 1, n}]", "Factorial[n]^2");
    check("Product[k, k]", "Factorial[-1 + k]");        /* indefinite -> (k-1)! */
    oracle("Product[k, {k, 1, n}]", "Product[k, {k, 1, %d}]");
    oracle("Product[k^2, {k, 1, n}]", "Product[k^2, {k, 1, %d}]");
    oracle("Product[k + a, {k, 1, n}]", "Product[k + a, {k, 1, %d}]");
    oracle("Product[(2 k - 1)/(2 k), {k, 1, n}]", "Product[(2 k - 1)/(2 k), {k, 1, %d}]");
    oracle("Product[k + 1, {k, 1, n}]", "Product[k + 1, {k, 1, %d}]");
    oracle("Product[(k + 1) (k + 2), {k, 1, n}]", "Product[(k + 1) (k + 2), {k, 1, %d}]");
    /* shifted lower bound */
    same("Product[x + i, {i, 0, n - 1}] /. {x -> 3, n -> 5}", "Product[3 + i, {i, 0, 4}]");
    same("Product[x - i, {i, 0, n - 1}] /. {x -> 7, n -> 4}", "Product[7 - i, {i, 0, 3}]");
    /* Method -> "Rational" forces the stage */
    same("Product[k, {k, 1, n}, Method -> \"Rational\"] /. n -> 6", "720");

    /* ---- Stage 3: Geometric (base^k via Sum exponent) ---- */
    check("Product[2^k, {k, 1, n}]", "2^(1/2 n (1 + n))");
    oracle("Product[2^k, {k, 1, n}]", "Product[2^k, {k, 1, %d}]");
    oracle("Product[a^k, {k, 0, n}]", "Product[a^k, {k, 0, %d}]");
    oracle("Product[k 2^k, {k, 1, n}]", "Product[k 2^k, {k, 1, %d}]");
    oracle("Product[k^2 a^k, {k, 0, n}]", "Product[k^2 a^k, {k, 0, %d}]");
    same("Product[a^i, i] /. {a -> 2, i -> 6}", "2^15");   /* indefinite geometric */
    /* HANG REGRESSION: a symbolic exponent must NOT reach Together/Factor;
     * this must return promptly (a^(n(n+1)/2)), never hang. */
    check("FreeQ[Product[a^k, {k, 1, n}], Product]", "True");
    same("Product[a^k, {k, 1, n}] /. {a -> 3, n -> 5}", "Product[3^k, {k, 1, 5}]");

    /* ---- Method polyalgorithm: per-method == Automatic in-class, falls
     *      through (unevaluated) out-of-class. ---- */
    same("Product[1 + 1/k, {k, 1, n}, Method -> \"Telescoping\"]",
         "Product[1 + 1/k, {k, 1, n}]");
    same("Product[k, {k, 1, n}, Method -> \"Rational\"]", "Product[k, {k, 1, n}]");
    same("Product[2^k, {k, 1, n}, Method -> \"Geometric\"]", "Product[2^k, {k, 1, n}]");
    /* out-of-class: Rational cannot do a geometric body -> stays held */
    check("Product[2^k, {k, 1, n}, Method -> \"Rational\"]",
          "Product[2^k, {k, 1, n}, Method -> \"Rational\"]");
    /* out-of-class: Geometric on a pure rational body -> stays held */
    check("Product[k, {k, 1, n}, Method -> \"Geometric\"]",
          "Product[k, {k, 1, n}, Method -> \"Geometric\"]");

    /* ---- Stage 4: infinite products + convergence (Product`Infinite) ---- */
    /* Weierstrass sine/cosine family (recognised in C). */
    check("Product[1 + 1/k^2, {k, 1, Infinity}]", "Sinh[Pi]/Pi");
    same("Product[1 + 5/i^2, {i, 1, Infinity}]", "Sinh[Sqrt[5] Pi]/(Sqrt[5] Pi)");
    /* numeric oracle for the z-form (Sqrt[-z^2/Pi^2] does not auto-reduce) */
    check("Chop[N[(Product[1 - z^2/(Pi^2 i^2), {i, 1, Infinity}] - Sin[z]/z) /. z -> 2]]",
          "0");
    /* rational-telescoping infinite product via the limit route */
    check("Product[1 - 1/k^2, {k, 2, Infinity}]", "1/2");
    /* divergent products stay held (Product::div is printed to stderr) */
    check("Product[1 + 1/k, {k, 1, Infinity}]", "Product[1 + 1/k, {k, 1, Infinity}]");
    check("Product[k, {k, 1, Infinity}]", "Product[k, {k, 1, Infinity}]");
    /* VerifyConvergence -> False skips the gate (still no closed form here, but
     * must not error and must stay well-defined) */
    check("Head[Product[k, {k, 1, Infinity}, VerifyConvergence -> False]]", "Product");

    /* ---- Stage 5: q-products (Product`QProduct -> QPochhammer) ---- */
    check("Product[1 - q^k, {k, 1, n}]", "QPochhammer[q, q, n]");
    check("Product[1 - a q^k, {k, 0, n - 1}]", "QPochhammer[a, q, n]");
    check("Product[1 - a q^k, {k, 1, n}]", "QPochhammer[a q, q, n]");
    check("Product[(1 - a q^k)/(1 - b q^k), {k, 0, n - 1}]",
          "QPochhammer[a, q, n]/QPochhammer[b, q, n]");
    oracle("Product[1 - q^k, {k, 1, n}]", "Product[1 - q^k, {k, 1, %d}]");
    same("Product[1 - a q^k, {k, 0, n - 1}] /. {a -> 3, q -> 1/2, n -> 4}",
         "Product[1 - 3 (1/2)^k, {k, 0, 3}]");
    same("Product[1 - q^k, {k, 1, n}, Method -> \"QProduct\"]", "QPochhammer[q, q, n]");

    /* ---- Stage 6: named special-function products ---- */
    check("Product[k^k, {k, 1, n}]", "Hyperfactorial[n]");
    check("Product[Gamma[i], {i, 1, n - 1}]", "BarnesG[n]");
    oracle("Product[k^k, {k, 1, n}]", "Product[k^k, {k, 1, %d}]");
    same("Product[Gamma[i], {i, 1, n - 1}] /. n -> 6", "BarnesG[6]");

    /* ---- Unsummable / non-rational body stays held ---- */
    check("Product[f[i], {i, 1, n}]", "Product[f[i], {i, 1, n}]");

    if (failures) {
        fprintf(stderr, "\n%d/%d product checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d product tests passed!\n", checks);
    return 0;
}
