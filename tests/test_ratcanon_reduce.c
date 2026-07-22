/* test_ratcanon_reduce.c — Phase 3 of RATCANON_REWRITE_PLAN.md.
 *
 * Validates rat_canon_reduce / rat_canon_normalize: one fmpz_mpoly_q combine +
 * reduce-mod-relation-ideal (flint_tower_reduce), accepting only a fully
 * radical-FREE result, with the §0.4 output convention.  Tests:
 *   - accepted forms: the exact canonical output on the covered regimes;
 *   - parity: where it does not decline, the result is math-equal to the
 *     classical Together/Cancel builtin (Together[rat - builtin] == 0);
 *   - declines: residual-radical / cube-root / forward-trig inputs return NULL
 *     (handed to the classical path in Phase 4);
 *   - the Sqrt[k] open-bug reproducer reduces.
 */
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "sym_names.h"
#include "ratcanon.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* rat_canon_normalize on the evaluated input (as the builtins receive it). */
static Expr* rcn(const char* input, RcMode mode) {
    Expr* e = parse_expression(input);
    Expr* v = evaluate(e);
    Expr* r = rat_canon_normalize(v, mode);
    expr_free(e); expr_free(v);
    return r;
}

static Expr* builtin_of(const char* input, RcMode mode) {
    Expr* e = parse_expression(input);
    Expr* call = expr_new_function(
        expr_new_symbol(mode == RCM_CANCEL ? SYM_Cancel : SYM_Together),
        (Expr*[]){ e }, 1);
    Expr* r = evaluate(call);   /* evaluate borrows */
    expr_free(call);
    return r;
}

static int is_zero(const Expr* e) {
    return e && ((e->type == EXPR_INTEGER && e->data.integer == 0) ||
                 (e->type == EXPR_REAL && e->data.real == 0.0));
}

/* Accepted: non-NULL result with the exact expected FullForm. */
static void accept_form(const char* input, RcMode mode, const char* expected) {
    Expr* r = rcn(input, mode);
    if (!r) { printf("FAIL accept: %s declined\n", input); ASSERT(r != NULL); return; }
    char* got = expr_to_string_fullform(r);
    if (strcmp(got, expected) != 0) {
        printf("FAIL form: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got); expr_free(r);
}

/* Declined: rat_canon_normalize returns NULL (classical fallback in Phase 4). */
static void declines(const char* input, RcMode mode) {
    Expr* r = rcn(input, mode);
    if (r) { char* s = expr_to_string(r); printf("FAIL decline: %s -> %s\n", input, s); free(s); }
    ASSERT(r == NULL);
    if (r) expr_free(r);
}

/* Parity: where not declined, math-equal to the classical builtin. */
static void parity(const char* input, RcMode mode) {
    Expr* r = rcn(input, mode);
    if (!r) return;                       /* declined — classical covers it */
    Expr* b = builtin_of(input, mode);
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_new_integer(-1), b }, 2);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ r, neg }, 2);
    Expr* tg = expr_new_function(expr_new_symbol(SYM_Together), (Expr*[]){ sum }, 1);
    Expr* diff = evaluate(tg);
    if (!is_zero(diff)) {
        char* s = expr_to_string(diff);
        printf("FAIL parity: %s -> nonzero diff %s\n", input, s);
        free(s);
    }
    ASSERT(is_zero(diff));
    expr_free(tg); expr_free(diff);
}

/* -------------------------------------------------------------------------- */

static void test_accept_forms(void) {
    accept_form("(1+x)/(1-x^2)", RCM_TOGETHER, "Times[-1, Power[Plus[-1, x], -1]]");
    accept_form("(x^2-1)/(x-1)", RCM_CANCEL, "Plus[1, x]");
    accept_form("1/x+1/y", RCM_TOGETHER,
                "Times[Power[x, -1], Power[y, -1], Plus[x, y]]");
    accept_form("1/(1+Log[x])+1/Log[x]", RCM_TOGETHER,
                "Times[Power[Plus[Log[x], Power[Log[x], 2]], -1], Plus[1, Times[2, Log[x]]]]");
    accept_form("(E^(2x)-1)/(E^x-1)", RCM_TOGETHER, "Plus[1, Power[E, x]]");
    accept_form("1/(x-I)+1/(x+I)", RCM_TOGETHER,
                "Times[2, x, Power[Plus[1, Power[x, 2]], -1]]");
    accept_form("1/(x-Sqrt[2])+1/(x+Sqrt[2])", RCM_TOGETHER,
                "Times[2, x, Power[Plus[-2, Power[x, 2]], -1]]");
    accept_form("1/(1+Tan[x])+1/(1-Tan[x])", RCM_TOGETHER,
                "Times[-2, Power[Plus[-1, Power[Tan[x], 2]], -1]]");
    /* Sqrt[k] open-bug reproducer: the conjugate sum eliminates the radical. */
    accept_form("1/(b-a Sqrt[k])+1/(b+a Sqrt[k])", RCM_TOGETHER,
                "Times[-2, b, Power[Plus[Times[-1, Power[b, 2]], Times[Power[a, 2], k]], -1]]");
    /* Cancel leaves a sum of fractions uncombined (per-term). */
    accept_form("a/b+c/d", RCM_CANCEL,
                "Plus[Times[a, Power[b, -1]], Times[c, Power[d, -1]]]");
    /* Phase 3b: VARIABLE-radicand pre-formed cancellations (radicand-var ordering
     * eliminates the radicand; the radical lands only in the numerator). */
    accept_form("(a^2-b)/(a-Sqrt[b])", RCM_CANCEL, "Plus[a, Power[b, Rational[1, 2]]]");
    accept_form("(y-1)/(y^(1/3)-1)", RCM_CANCEL,
                "Plus[1, Power[y, Rational[1, 3]], Power[y, Rational[2, 3]]]");
    accept_form("(x^3-y)/(x-y^(1/3))", RCM_CANCEL,
                "Plus[Power[x, 2], Times[x, Power[y, Rational[1, 3]]], Power[y, Rational[2, 3]]]");
    /* Phase 3c: CONSTANT-radicand pre-formed cancellations via number-field GCD. */
    accept_form("(x^2-2)/(x-Sqrt[2])", RCM_CANCEL, "Plus[Power[2, Rational[1, 2]], x]");
    accept_form("(x^2-3)/(x+Sqrt[3])", RCM_CANCEL,
                "Plus[Times[-1, Power[3, Rational[1, 2]]], x]");
}

static void test_declines(void) {
    /* Coprime / WL-kept radical in the denominator: classical gives the cleaner
     * form; decline. */
    declines("1/(x-Sqrt[2])", RCM_TOGETHER);             /* radical WL-kept in denom       */
    declines("1/(x-Sqrt[3])+1/(x-Sqrt[5])", RCM_TOGETHER); /* coprime multi-radical sum    */
    declines("1/(1+Cos[x])+1/Sin[x]", RCM_TOGETHER);     /* forward trig (dependent)       */
}

static void test_parity(void) {
    const char* tog[] = {
        "(1+x)/(1-x^2)", "1/x+1/y", "x/(x+1)+1/(x+1)", "1/(x^2-1)+1/(x+1)",
        "2/(x-1)-2/(x+1)", "1/(a+b)+1/(a-b)", "b/(x-a)+b/(x+a)", "a/(x-1)+b/(x-2)",
        "1/(1+Log[x])+1/Log[x]", "1/Log[x]-1/Log[y]", "(E^(2x)-1)/(E^x-1)",
        "1/(E^x-1)+1/(E^x+1)", "1/(1+Tan[x])+1/(1-Tan[x])",
        "1/(x-I)+1/(x+I)", "1/(1-I x)+1/(1+I x)",
        "1/(x-Sqrt[2])+1/(x+Sqrt[2])", "1/(b-a Sqrt[k])+1/(b+a Sqrt[k])",
        /* declined ones return NULL and are skipped: */
        "(x^2-2)/(x-Sqrt[2])", "(y-1)/(y^(1/3)-1)", "1/(x-Sqrt[2])",
        NULL };
    for (int i = 0; tog[i]; i++) parity(tog[i], RCM_TOGETHER);
    const char* can[] = {
        "(x^2-1)/(x-1)", "(x^3-1)/(x-1)", "(2 x^2+4x)/(x^2-4)", "a/b+c/d",
        "(x^2+2x+1)/(x+1)", "(a x^2 - a)/(x-1)",
        /* Phase 3b variable-radicand cancellations: */
        "(a^2-b)/(a-Sqrt[b])", "(y-1)/(y^(1/3)-1)", "(x^3-y)/(x-y^(1/3))",
        "(k x^2-1)/(x-1/Sqrt[k])",
        /* Phase 3c constant-radicand cancellations: */
        "(x^2-2)/(x-Sqrt[2])", "(x^2-3)/(x+Sqrt[3])",
        "(x-Sqrt[2])(x+Sqrt[2])/(x-Sqrt[2])", NULL };
    for (int i = 0; can[i]; i++) parity(can[i], RCM_CANCEL);
}

int main(void) {
    symtab_init();
    core_init();
    printf("Running ratcanon reduce (Phase 3) tests...\n");
    TEST(test_accept_forms);
    TEST(test_declines);
    TEST(test_parity);
    printf("All ratcanon reduce tests passed!\n");
    return 0;
}
