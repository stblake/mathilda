/*
 * test_sum.c -- Sum (Stages 0-3) and DifferenceDelta.
 *
 * Uses an explicit aborting checker rather than assert_eval_eq because libc
 * assert() compiles out under -DNDEBUG (Release).  check() compares the
 * InputForm string of the evaluated input against an expected string;
 * oracle() cross-checks a symbolic closed form against direct finite
 * expansion at a concrete bound.
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

/* Cross-check mathematical equality: Simplify[(a) - (b)] must be 0.  Robust to
 * the two sides being printed in different (but equivalent) forms. */
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

int main(void) {
    symtab_init();
    core_init();

    /* ---- Stage 0: finite explicit expansion ---- */
    check("Sum[i, {i, 1, 10}]", "55");
    check("Sum[i^2, {i, 1, 100}]", "338350");
    check("Sum[2^i, {i, 0, 10}]", "2047");
    check("Sum[1/i, {i, 1, 5}]", "137/60");
    check("Sum[i, {i, 2, 10, 2}]", "30");
    check("Sum[i, {i, {1, 3, 5}}]", "9");
    check("Sum[i, {i, {2, 3, 5, 7}}]", "17");
    check("Sum[i j, {i, 1, 3}, {j, 1, 3}]", "36");
    /* outer bound depends on inner variable */
    check("Sum[f[i, j], {i, 1, 3}, {j, 1, i}]",
          "f[1, 1] + f[2, 1] + f[2, 2] + f[3, 1] + f[3, 2] + f[3, 3]");
    check("Sum[f[i], {i, 1, 5}]", "f[1] + f[2] + f[3] + f[4] + f[5]");
    /* symbolic body over a symbolic list */
    check("Sum[f[i], {i, {a, b, c}}]", "f[a] + f[b] + f[c]");
    /* unsummable symbolic-bound case stays held */
    check("Sum[f[i], {i, 1, n}]", "Sum[f[i], {i, 1, n}]");

    /* ---- Stage 1: polynomial ---- */
    check("Sum[i, {i, 1, n}]", "1/2 n (1 + n)");
    check("Sum[i^2, {i, 1, n}]", "1/6 n (1 + n) (1 + 2 n)");
    check("Sum[i^3, {i, 1, n}]", "1/4 n^2 (1 + n)^2");
    check("Sum[i^2, i]", "1/6 i (-1 + i) (-1 + 2 i)");
    check("Sum[i, i]", "1/2 i (-1 + i)");
    /* oracle: closed form at n=7 equals direct finite sum */
    same("Sum[i^2, {i, 1, n}] /. n -> 7", "Sum[i^2, {i, 1, 7}]");
    same("Sum[i^4, {i, 1, n}] /. n -> 9", "Sum[i^4, {i, 1, 9}]");
    same("Sum[i^3 + (i + 3)^5, {i, 1, n}] /. n -> 6",
         "Sum[i^3 + (i + 3)^5, {i, 1, 6}]");

    /* ---- Stage 2: geometric / polynomial-exponential ---- */
    check("Sum[a^i, i]", "a^i/(-1 + a)");
    check("Sum[2^i, i]", "2^i");
    check("Sum[q1^i q2^i, i]", "(q1 q2)^i/(-1 + q1 q2)");
    same("Sum[a^i, {i, 0, n}] /. {a -> 3, n -> 5}", "Sum[3^i, {i, 0, 5}]");
    same("Sum[i 2^i, {i, 1, n}] /. n -> 6", "Sum[i 2^i, {i, 1, 6}]");
    same("Sum[i^2 a^i, {i, 0, n}] /. {a -> 2, n -> 4}",
         "Sum[i^2 2^i, {i, 0, 4}]");

    /* ---- Stage 3: Gosper + DifferenceDelta ---- */
    check("Sum[k k!, k]", "Factorial[k]");
    check("Sum[k k!, {k, 1, n}]", "-1 + Factorial[1 + n]");
    check("Sum[1/(i (i + 1)), i]", "-1/i");
    same("Sum[1/(i (i + 1)), {i, 1, n}] /. n -> 8", "Sum[1/(i (i + 1)), {i, 1, 8}]");
    same("Sum[1/(i (i + 1) (i + 2)), {i, 1, n}] /. n -> 7",
         "Sum[1/(i (i + 1) (i + 2)), {i, 1, 7}]");
    /* DifferenceDelta is the left inverse of indefinite Sum */
    same("DifferenceDelta[Sum[k k!, k], k]", "k k!");
    same("DifferenceDelta[Sum[1/(i (i + 1)), i], i]", "1/(i (i + 1))");
    check("DifferenceDelta[i^2, i]", "1 + 2 i");

    if (failures) {
        fprintf(stderr, "\n%d/%d sum checks FAILED\n", failures, checks);
        return 1;
    }
    printf("All %d sum tests passed!\n", checks);
    return 0;
}
