/* test_cherry_li.c — Cherry's logarithmic integral: LogIntegral (li) engine.
 *
 * Pins for cherry_li.c (CHERRY_PLAN.md §5): the single-log tower multi-li cases
 * from G. W. Cherry, "Integration in Finite Terms with Special Functions: The
 * Logarithmic Integral" (1986), plus the decline-safety cases.  Because a plain
 * Simplify diff-back cannot reduce Log[w^k] = k Log[w] (branch-cut safety),
 * correctness is checked by an exact PowerExpand diff-back AND, independently, a
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

static Expr* eval_str(const char* s) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}
static bool eval_is(const char* s, const char* expected) {
    Expr* r = eval_str(s);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got, expected);
    free(got);
    expr_free(r);
    return ok;
}

/* Integrate f via RischTranscendental; assert it closes to a LogIntegral form
 * that (a) PowerExpand-diff-backs to 0 exactly and (b) numeric-diff-backs to 0 at
 * an interior point. */
static void assert_li(const char* f) {
    char buf[1500];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x, Method -> \"RischTranscendental\"]},"
        " Head[r] =!= Integrate && !FreeQ[r, LogIntegral]]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a closed LogIntegral form", f);
    snprintf(buf, sizeof(buf),
        "Simplify[PowerExpand[D[Integrate[%s, x, Method -> \"RischTranscendental\"], x]"
        " - (%s)]]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: PowerExpand diff-back nonzero", f);
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate[%s, x, Method -> \"RischTranscendental\"], x] - (%s))"
        " /. x -> 17/10]] < 1/10^6", f, f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: numeric diff-back nonzero", f);
}

static void assert_declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Head[Integrate`RischTranscendental[%s, x]] === Integrate`RischTranscendental", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a clean decline", f);
}

/* Cherry 1986 paper pins. */
static void test_cherry_1986_li(void) {
    /* Ex 1.1 / d1: INT x/Log[x]^2 dx = 2 li(x^2) - x^2/Log[x]. */
    assert_li("x/Log[x]^2");
    ASSERT_MSG(eval_is(
        "Simplify[PowerExpand[Integrate[x/Log[x]^2, x, Method -> \"RischTranscendental\"]"
        " - (2 LogIntegral[x^2] - x^2/Log[x])]]", "0"), "d1 exact form");
    /* d3: INT x^2/Log[x+1] dx = li((x+1)^3) - 2 li((x+1)^2) + li(x+1). */
    assert_li("x^2/Log[x+1]");
    /* Ex 5.1: INT x^3/Log[x^2-1] dx = (1/2) li(x^2-1) + (1/2) li((x^2-1)^2). */
    assert_li("x^3/Log[x^2-1]");
}

/* Broader single-log multi-li coverage. */
static void test_cherry_li_stress(void) {
    assert_li("1/Log[x]");            /* li(x)                    */
    assert_li("x/Log[x]");            /* li(x^2)                  */
    assert_li("1/Log[2 x]");          /* affine argument          */
    assert_li("x^2/Log[x-1]");        /* shifted                  */
    assert_li("x^3/Log[x^2]");        /* w = x^2                  */
    assert_li("x^4/Log[x+2]");        /* multi-li                 */
    assert_li("x^5/Log[x+1]");
    assert_li("(x+1)/Log[x]^2");
    assert_li("x^2/Log[3 x + 1]");
    assert_li("x/Log[x^2+1]");
    assert_li("2/Log[x]^3");          /* m = 3                    */
}

/* Transcendental-constant rescaling: a constant root rho of the theta-denominator
 * gives a rescaled li(e^(-rho) w) term (Log[w] - rho = Log[e^(-rho) w]). */
static void test_cherry_li_rescale(void) {
    /* d2: INT 1/(Log[x]+3) dx = e^-3 li(e^3 x). */
    assert_li("1/(Log[x]+3)");
    ASSERT_MSG(eval_is(
        "Simplify[PowerExpand[Integrate[1/(Log[x]+3), x, Method -> \"RischTranscendental\"]"
        " - LogIntegral[E^3 x]/E^3]]", "0"), "d2 exact form");
    /* d4: INT (Log[x]^2+3)/(Log[x]^2+3 Log[x]+2) dx
     *       = x + 4 e^-1 li(e x) - 7 e^-2 li(e^2 x). */
    assert_li("(Log[x]^2+3)/(Log[x]^2+3 Log[x]+2)");
    assert_li("1/(Log[x]+1)");
    assert_li("1/((Log[x]+1)(Log[x]+2))");
}

/* Decline-safety (later increments): multi-log towers requiring a product
 * Sigma-decomposition, and a non-polynomial log argument — decline cleanly. */
static void test_cherry_li_declines(void) {
    assert_declines("x/Log[Sin[x]]");       /* non-polynomial w                */
    assert_declines("x^3/Log[x^3-x]");      /* reducible w -> product args     */
}

int main(void) {
    core_init();
    TEST(test_cherry_1986_li);
    TEST(test_cherry_li_stress);
    TEST(test_cherry_li_rescale);
    TEST(test_cherry_li_declines);
    printf("All Cherry LogIntegral tests passed.\n");
    return 0;
}
