/* test_cherry_dilog.c — dilogarithm integration (Cherry degree-2 Sigma-decomp).
 *
 * R(x) Log[w] -> Log-Log products + PolyLog[2, g] (cherry_dilog.c).  Because a
 * plain Simplify diff-back cannot reduce the log-branch relations, correctness is
 * checked by an exact PowerExpand diff-back AND a numeric interior-point diff-back.
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

/* Integrate[f] closes to a PolyLog form; PowerExpand + numeric diff-back verified. */
static void assert_dilog(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x]}, Head[r] =!= Integrate && !FreeQ[r, PolyLog]]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a PolyLog form", f);
    snprintf(buf, sizeof(buf),
        "Simplify[PowerExpand[D[Integrate[%s, x], x] - (%s)]]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: PowerExpand diff-back nonzero", f);
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate[%s, x], x] - (%s)) /. x -> 13/10]] < 1/10^6", f, f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: numeric diff-back nonzero", f);
}

static void test_dilog(void) {
    /* generalisation of rt_try_dilog to the Log-Log + PolyLog[2] form. */
    assert_dilog("Log[x]/(1+x)");     /* = Log[x]Log[1+x] + PolyLog[2,-x] */
    assert_dilog("Log[x]/(1-x)");     /* = PolyLog[2, 1-x]                */
    assert_dilog("Log[x]/(x^2-1)");   /* factored denominator            */
    assert_dilog("Log[2 x + 1]/(x+1)");
    /* exact flagship form */
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[Log[x]/(1+x), x] - (Log[x] Log[1+x] + PolyLog[2, -x])]", "0"),
        "Log[x]/(1+x) exact form");
    /* rt_try_dilog fast-path forms still fire */
    assert_dilog("Log[1+x]/x");       /* -PolyLog[2,-x] */
    assert_dilog("Log[1-x]/x");
}

/* Decline-safety: transcendental-constant root spacings and products of logs
 * (degree > 1 in the log tower) are later increments; they decline cleanly. */
static void test_dilog_declines(void) {
    ASSERT_MSG(eval_is("Head[Integrate`RischTranscendental[Log[2+x]/x, x]]"
                       " === Integrate`RischTranscendental", "True"),
        "Log[2+x]/x should decline");
    ASSERT_MSG(eval_is("Head[Integrate`RischTranscendental[Log[x] Log[1+x]/x, x]]"
                       " === Integrate`RischTranscendental", "True"),
        "Log[x]Log[1+x]/x should decline");
}

int main(void) {
    core_init();
    TEST(test_dilog);
    TEST(test_dilog_declines);
    printf("All Cherry dilogarithm tests passed.\n");
    return 0;
}
