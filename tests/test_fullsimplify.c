/*
 * test_fullsimplify.c -- framework tests for FullSimplify.
 *
 * Exercises the FullSimplify driver and relevance engine: genuine gaps that
 * Simplify alone does not close, the >=-Simplify guarantee, option handling
 * (ComplexityFunction, TransformationFunctions, Assumptions, TimeConstraint),
 * and protection. The per-identity corpus lives in
 * test_fullsimplify_corpus.c / fullsimplify_corpus.m.
 *
 * Uses a hard-failing checker (exit(1)) rather than test_utils.h's
 * assert_eval_eq, whose libc assert() is compiled out under -DNDEBUG.
 */

#include "expr.h"
#include "symtab.h"
#include "parse.h"
#include "eval.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void core_init(void);

static int failures = 0;

/* Assert that evaluating `input` prints exactly `expected`. */
static void chk(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: parse: %s\n", input); failures++; return; }
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(res);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
        failures++;
    }
    free(s);
    expr_free(res);
}

/* Assert that evaluating `input` prints a string containing `needle`. */
static void chk_contains(const char* input, const char* needle) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: parse: %s\n", input); failures++; return; }
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(res);
    if (!strstr(s, needle)) {
        fprintf(stderr, "FAIL: %s\n  expected to contain: %s\n  actual: %s\n",
                input, needle, s);
        failures++;
    }
    free(s);
    expr_free(res);
}

int main(void) {
    symtab_init();
    core_init();

    /* --- The driver is available and protected. --- */
    chk("Attributes[FullSimplify]", "{Protected}");

    /* --- Genuine gaps: FullSimplify reduces where Simplify does not. --- */
    chk("Simplify[Gamma[x+1]/Gamma[x]]", "Gamma[1 + x]/Gamma[x]"); /* sanity: gap */
    chk("FullSimplify[Gamma[x+1]/Gamma[x]]", "x");
    chk("FullSimplify[Gamma[x+1] - x Gamma[x]]", "0");
    chk("FullSimplify[LogGamma[x+1] - LogGamma[x]]", "Log[x]");
    chk("FullSimplify[PolyGamma[0, x+1] - PolyGamma[0, x]]", "1/x");
    chk("FullSimplify[Erf[x] + Erfc[x]]", "1");
    chk("FullSimplify[PolyLog[2, z] + PolyLog[2, -z]]", "1/2 PolyLog[2, z^2]");
    chk("FullSimplify[Surd[x, 3]^3]", "x");
    /* Conversion rule (Pochhammer -> Gamma) enables a downstream cancellation. */
    chk("FullSimplify[Pochhammer[a, n]/Gamma[a+n]]", "1/Gamma[a]");

    /* --- >= Simplify: never worse on inputs Simplify already handles. --- */
    chk("FullSimplify[(x-1)(x+1)(x^2+1)+1]", "x^4");
    chk("FullSimplify[Sin[x]^2 + Cos[x]^2]", "1");

    /* --- Assumptions (positional). --- */
    chk("FullSimplify[Log[E^(x+y)], Element[{x, y}, Reals]]", "x + y");
    chk("FullSimplify[Sqrt[x^2], x > 0]", "x");

    /* --- TransformationFunctions: user transforms merge with the relevance
     * set (Automatic is always kept). --- */
    chk("FullSimplify[a + b, TransformationFunctions -> {(# /. a -> 0 &)}]", "b");

    /* --- ComplexityFunction is forwarded to Simplify. --- */
    chk_contains("FullSimplify[100 Log[2], ComplexityFunction -> LeafCount]", "Log");

    /* --- TimeConstraint -> {tLoc, tTot} | t | Infinity all preserve the
     * result for fast transforms (the budgets are plumbed, not tripped). --- */
    chk("FullSimplify[Gamma[x+1]/Gamma[x], TimeConstraint -> {1, 5}]", "x");
    chk("FullSimplify[Gamma[x+1]/Gamma[x], TimeConstraint -> 5]", "x");
    chk("FullSimplify[Gamma[x+1]/Gamma[x], TimeConstraint -> Infinity]", "x");

    /* --- Relevance: a plain expression with no special-function heads still
     * simplifies (and behaves exactly like Simplify). --- */
    chk("FullSimplify[a x + b x + c]", "c + (a + b) x");

    if (failures == 0) {
        printf("All FullSimplify framework tests passed.\n");
        return 0;
    }
    fprintf(stderr, "%d FullSimplify framework test(s) failed.\n", failures);
    return 1;
}
