/* test_cherry_dilog_exp.c — exponential-tower dilogarithm integration (Cherry).
 *
 * Rational-in-E^(cx) x-weighted forms (x/(E^x-1)) and outer-log forms
 * (Log[1+E^x]) -> PolyLog[2, exp] (cherry_dilog_exp.c).  A plain Simplify
 * diff-back cannot reduce the exp-log branch relations (Log[1-E^-x] =
 * Log[E^x-1] - x), so correctness is checked by a PowerExpand diff-back that
 * first splits logs of quotients (Log[a] -> Log[Factor[Together[a]]]) AND a
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

/* Integrate[f] closes to a PolyLog form; PowerExpand + numeric diff-back verified.
 * The Together inside the Log split is what turns E^-x = 1/theta so
 * Log[1-E^-x] = Log[(E^x-1)/E^x] reduces — a plain Log[Factor[a]] cannot. */
static void assert_dilog_exp(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x]}, Head[r] =!= Integrate && !FreeQ[r, PolyLog]]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a PolyLog form", f);
    snprintf(buf, sizeof(buf),
        "Simplify[PowerExpand[(D[Integrate[%s, x], x] - (%s)) "
        "/. Log[a_] :> Log[Factor[Together[a]]]]]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: PowerExpand diff-back nonzero", f);
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate[%s, x], x] - (%s)) /. x -> 13/10]] < 1/10^6", f, f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: numeric diff-back nonzero", f);
}

static void test_dilog_exp(void) {
    /* rational-in-E^(cx) x-weighted: the flagship family. */
    assert_dilog_exp("x/(-1 + E^x)");      /* = x Log[1-E^-x] - PolyLog[2, E^-x] */
    assert_dilog_exp("x/(1 + E^x)");
    assert_dilog_exp("x E^x/(-1 + E^x)");
    assert_dilog_exp("x/(-1 + E^(2 x))");  /* c = 2 rate */
    /* outer-log forms (weight-1 log = the outer Log[theta - rho]). */
    assert_dilog_exp("Log[1 + E^x]");      /* = -PolyLog[2, -E^x] */
    assert_dilog_exp("Log[1 + E^(-x)]");
    /* mixture: a linear combination of the two shapes. */
    assert_dilog_exp("2 x/(-1 + E^x) + 3 Log[1 + E^x]");

    /* exact flagship form (numeric: Log[1-E^-x] vs Log[E^x-1]-x confounds Simplify). */
    ASSERT_MSG(eval_is(
        "Abs[N[(Integrate[x/(-1 + E^x), x] "
        "- (x Log[1 - E^-x] - PolyLog[2, E^-x])) /. x -> 13/10]] < 1/10^6", "True"),
        "x/(E^x-1) exact form");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[Log[1 + E^x], x] - (-PolyLog[2, -E^x])]", "0"),
        "Log[1+E^x] exact form");
}

/* Decline-safety.  Weight-2 integrands (-> PolyLog[3]) and elementary exp
 * integrands must NOT be turned into a spurious dilog. */
static void test_dilog_exp_declines(void) {
    /* weight-2 OUTER-log x Log[1+E^x] -> trilog is out of scope for both the
     * dilog and the polylog engines (the latter is rational-in-theta only). */
    ASSERT_MSG(eval_is("Head[Integrate`RischTranscendental[x Log[1 + E^x], x]]"
                       " === Integrate`RischTranscendental", "True"),
        "x Log[1+E^x] should decline (weight-2 outer log)");
    /* x^2/(E^x-1) is NOT out of scope any more — the general polylog-ladder
     * engine (cherry_polylog_exp) closes it to a PolyLog[3] form. */
    ASSERT_MSG(eval_is("!FreeQ[Integrate[x^2/(-1 + E^x), x], PolyLog[3, _]]", "True"),
        "x^2/(E^x-1) now closes via the polylog ladder");
    /* elementary exp integrands stay elementary (handled before the Cherry stage). */
    ASSERT_MSG(eval_is("FreeQ[Integrate[1/(-1 + E^x), x], PolyLog]", "True"),
        "1/(E^x-1) should stay elementary");
    ASSERT_MSG(eval_is("FreeQ[Integrate[E^x/(-1 + E^x), x], PolyLog]", "True"),
        "E^x/(E^x-1) should stay elementary");
}

/* No regression on the LOG-tower dilog engine: a DerivativeDivides u=Log[x]
 * substitution must not route these through the exp-tower engine (which is gated
 * to the top level), so the cleaner PolyLog[2, 1-x] form survives. */
static void test_log_tower_unregressed(void) {
    ASSERT_MSG(eval_is("Simplify[Integrate[Log[x]/(1 - x), x] - PolyLog[2, 1 - x]]", "0"),
        "Log[x]/(1-x) must stay PolyLog[2,1-x]");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[Log[x]/(1 + x), x] - (Log[x] Log[1 + x] + PolyLog[2, -x])]", "0"),
        "Log[x]/(1+x) must stay clean");
}

/* The direct debug surfaces Integrate`Cherry`<Name>[f, x] apply each engine
 * bypassing the cascade — confirm the DilogExp surface is registered and fires. */
static void test_direct_surface(void) {
    ASSERT_MSG(eval_is(
        "With[{r = Integrate`Cherry`DilogExp[x/(-1 + E^x), x]}, "
        "Head[r] =!= Integrate`Cherry`DilogExp && !FreeQ[r, PolyLog]]", "True"),
        "Integrate`Cherry`DilogExp direct surface should fire");
    ASSERT_MSG(eval_is(
        "With[{r = Integrate`Cherry`Dilog[Log[x]/(1 + x), x]}, !FreeQ[r, PolyLog]]", "True"),
        "Integrate`Cherry`Dilog direct surface should fire");
    ASSERT_MSG(eval_is(
        "With[{r = Integrate`Cherry`Ei[E^(1/x), x]}, !FreeQ[r, ExpIntegralEi]]", "True"),
        "Integrate`Cherry`Ei direct surface should fire");
    ASSERT_MSG(eval_is(
        "With[{r = Integrate`Cherry`Li[x/Log[x]^2, x]}, !FreeQ[r, LogIntegral]]", "True"),
        "Integrate`Cherry`Li direct surface should fire");
}

int main(void) {
    core_init();
    TEST(test_dilog_exp);
    TEST(test_dilog_exp_declines);
    TEST(test_log_tower_unregressed);
    TEST(test_direct_surface);
    printf("All Cherry exponential-tower dilogarithm tests passed.\n");
    return 0;
}
