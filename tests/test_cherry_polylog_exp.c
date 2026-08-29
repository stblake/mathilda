/* test_cherry_polylog_exp.c — exponential-tower polylogarithm-ladder integration.
 *
 * x^n / Q(E^(c x)) -> polylogarithms up to weight n+1 (cherry_polylog_exp.c),
 * over rational AND algebraic roots of Q.  A plain Simplify diff-back cannot
 * reduce the exp-log branch relations, so correctness is a PowerExpand diff-back
 * that splits logs of quotients (Log[a] -> Log[Factor[Together[a]]]) plus a
 * numeric interior-point diff-back.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

static bool eval_is(const char* s, const char* expected) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got, expected);
    free(got); expr_free(r); expr_free(e);
    return ok;
}

/* Integrate[f] closes to a PolyLog form of the stated top weight; PowerExpand +
 * numeric diff-back verified. */
static void assert_polylog(const char* f, int top_weight) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x]}, Head[r] =!= Integrate "
        "&& !FreeQ[r, PolyLog[%d, _]]]", f, top_weight);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a PolyLog[%d,..] form", f, top_weight);
    snprintf(buf, sizeof(buf),
        "Simplify[PowerExpand[(D[Integrate[%s, x], x] - (%s)) "
        "/. Log[a_] :> Log[Factor[Together[a]]]]]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: PowerExpand diff-back nonzero", f);
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate[%s, x], x] - (%s)) /. x -> 13/10]] < 1/10^6", f, f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: numeric diff-back nonzero", f);
}

static void test_polylog_ladder(void) {
    /* higher weight over a rational root (weight 3, 5) */
    assert_polylog("x^2/(-1 + E^x)", 3);        /* -> ... - 2 PolyLog[3, E^-x] */
    assert_polylog("x^3/(-1 + E^x)", 4);
    assert_polylog("x^4/(-1 + E^(5 x))", 5);     /* rate c = 5 */
    /* algebraic roots (Q(Sqrt[5]) from theta^2 + theta - 1) */
    assert_polylog("x/(-1 + E^x + E^(2 x))", 2);
    assert_polylog("x^2/(-1 + E^x + E^(2 x))", 3);

    /* exact form of the x^4/(E^(5x)-1) ladder (numeric: E^-5x vs E^5x branch). */
    ASSERT_MSG(eval_is(
        "Abs[N[(Integrate[x^4/(-1 + E^(5 x)), x] - ("
        "1/5 x^4 Log[1 - E^(-5 x)] - 4/25 x^3 PolyLog[2, E^(-5 x)] "
        "- 12/125 x^2 PolyLog[3, E^(-5 x)] - 24/625 x PolyLog[4, E^(-5 x)] "
        "- 24/3125 PolyLog[5, E^(-5 x)])) /. x -> 7/10]] < 1/10^6", "True"),
        "x^4/(E^(5x)-1) exact ladder form");
}

/* Decline-safety: elementary exp integrands and non-P(x)/Q(theta) shapes must
 * not be turned into a spurious polylog. */
static void test_polylog_declines(void) {
    ASSERT_MSG(eval_is("FreeQ[Integrate[1/(-1 + E^x), x], PolyLog]", "True"),
        "1/(E^x-1) should stay elementary");
    ASSERT_MSG(eval_is("FreeQ[Integrate[E^x/(-1 + E^x), x], PolyLog]", "True"),
        "E^x/(E^x-1) should stay elementary");
    /* outer-log integrand is handled by the dilog engine, not this one. */
    ASSERT_MSG(eval_is(
        "Head[Integrate`Cherry`PolyLogExp[Log[1 + E^x], x]] "
        "=== Integrate`Cherry`PolyLogExp", "True"),
        "PolyLogExp should decline the outer-log form");
    /* direct debug surface fires on a genuine case */
    ASSERT_MSG(eval_is(
        "With[{r = Integrate`Cherry`PolyLogExp[x^2/(-1 + E^x), x]}, "
        "!FreeQ[r, PolyLog[3, _]]]", "True"),
        "Integrate`Cherry`PolyLogExp direct surface should fire");
}

/* No regression on the LOG-tower dilog engine (top-level depth gate). */
static void test_log_tower_unregressed(void) {
    ASSERT_MSG(eval_is("Simplify[Integrate[Log[x]/(1 - x), x] - PolyLog[2, 1 - x]]", "0"),
        "Log[x]/(1-x) must stay PolyLog[2,1-x]");
}

int main(void) {
    core_init();
    TEST(test_polylog_ladder);
    TEST(test_polylog_declines);
    TEST(test_log_tower_unregressed);
    printf("All Cherry exponential-tower polylogarithm tests passed.\n");
    return 0;
}
