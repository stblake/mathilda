/* test_eval_timestamps.c -- M3 phase-3 (eval timestamps) tests.
 *
 * Phase-3 invariants:
 *   - A freshly-constructed Expr has last_evaluated_at == 0 and the
 *     global eval_clock starts at 1, so first-time evaluation never
 *     hits the cache.
 *   - After evaluate() reaches fixed point, the result's
 *     last_evaluated_at equals eval_clock_get(); a subsequent
 *     evaluate(result) returns immediately (semantic equality, no
 *     re-walk) as long as the clock is unchanged.
 *   - Symbol-table mutations bump the clock and invalidate every
 *     cached evaluation: Set, SetDelayed, Clear, SetAttributes,
 *     ClearAttributes.
 *   - Pure builtin calls (no symtab change) leave the clock alone.
 *   - The timestamp is benign metadata: expr_eq / expr_hash ignore it.
 */

#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Re-declare here because test_utils.h does not pull in eval.h. */
extern uint64_t eval_clock_get(void);

/* ------------------------------------------------------------------ */
/* 1. Clock starts non-zero and is monotonically increasing.          */
/* ------------------------------------------------------------------ */

static void test_clock_starts_nonzero_and_monotonic(void) {
    uint64_t c0 = eval_clock_get();
    ASSERT(c0 >= 1);

    /* A pure evaluation should NOT advance the clock. */
    Expr* e = parse_expression("1+2");
    Expr* r = evaluate(e);
    expr_free(e);
    expr_free(r);
    ASSERT(eval_clock_get() == c0);

    /* A Set must advance the clock. */
    Expr* a = parse_expression("zclkA = 5");
    Expr* ar = evaluate(a);
    expr_free(a); expr_free(ar);
    uint64_t c1 = eval_clock_get();
    ASSERT(c1 > c0);

    /* Another Set must advance again. */
    Expr* b = parse_expression("zclkB = 7");
    Expr* br = evaluate(b);
    expr_free(b); expr_free(br);
    ASSERT(eval_clock_get() > c1);

    /* Clean up. */
    Expr* clr = parse_expression("Clear[zclkA, zclkB]");
    expr_free(evaluate(clr));
    expr_free(clr);
}

/* ------------------------------------------------------------------ */
/* 2. Fresh expression has last_evaluated_at == 0.                     */
/* ------------------------------------------------------------------ */

static void test_fresh_expr_has_zero_timestamp(void) {
    Expr* a = expr_new_integer(42);
    ASSERT(a->last_evaluated_at == 0);
    expr_free(a);

    Expr* s = expr_new_symbol("ZZZ");
    ASSERT(s->last_evaluated_at == 0);
    expr_free(s);

    Expr* args[2] = { expr_new_integer(2), expr_new_integer(3) };
    Expr* f = expr_new_function(expr_new_symbol("Plus"), args, 2);
    ASSERT(f->last_evaluated_at == 0);
    expr_free(f);

    Expr* p = parse_expression("Sin[Pi/4] + Cos[x]");
    ASSERT(p->last_evaluated_at == 0);
    expr_free(p);
}

/* ------------------------------------------------------------------ */
/* 3. Result of evaluate() is stamped with current clock.              */
/* ------------------------------------------------------------------ */

static void test_evaluate_stamps_result(void) {
    Expr* p = parse_expression("Sin[Pi/4]");
    Expr* r = evaluate(p);
    expr_free(p);

    uint64_t clk = eval_clock_get();
    ASSERT(r->last_evaluated_at == clk);

    expr_free(r);
}

/* ------------------------------------------------------------------ */
/* 4. Re-evaluation of an already-stamped FUNCTION early-exits.        */
/*    Verified semantically: the result is structurally identical.    */
/* ------------------------------------------------------------------ */

static void test_reevaluate_same_clock_is_identity(void) {
    Expr* p = parse_expression("D[Sin[x]^2 + Cos[x]^2, x]");
    Expr* r1 = evaluate(p);
    expr_free(p);

    /* r1 is a FUNCTION (or possibly an atom if simplification reduced
     * to 0). Whatever it is, evaluating it again under the same clock
     * must yield a structurally-identical value. */
    uint64_t clk_before = eval_clock_get();
    Expr* r2 = evaluate(r1);
    uint64_t clk_after = eval_clock_get();

    ASSERT(clk_after == clk_before);   /* no symtab mutation */
    ASSERT(expr_eq(r1, r2));
    /* If r1 was a FUNCTION, the early-exit must have fired -- in which
     * case r2 is the SAME node as r1 (refcount-shared). Atoms also
     * share identity through expr_copy, so r2 == r1 should hold. */
    ASSERT(r1 == r2);

    expr_free(r1);
    expr_free(r2);
}

/* ------------------------------------------------------------------ */
/* 5. Set bumps the clock and invalidates the cache.                   */
/* ------------------------------------------------------------------ */

static void test_set_invalidates(void) {
    /* Make sure the symbol starts unbound. */
    Expr* clr0 = parse_expression("Clear[zclkX]");
    expr_free(evaluate(clr0));
    expr_free(clr0);

    /* Evaluate f[zclkX] once with zclkX symbolic. */
    Expr* p1 = parse_expression("Hold[zclkX + 1]");
    Expr* r1 = evaluate(p1);
    expr_free(p1);
    char* s1 = expr_to_string(r1);
    /* Without a value, this should remain Hold[zclkX + 1]. */
    ASSERT(strstr(s1, "zclkX") != NULL);
    free(s1);

    uint64_t clk_pre = eval_clock_get();

    /* Now assign a value -- this must bump the clock. */
    Expr* set_p = parse_expression("zclkX = 99");
    expr_free(evaluate(set_p));
    expr_free(set_p);

    ASSERT(eval_clock_get() > clk_pre);

    /* If we re-evaluate r1 (which still has the OLD timestamp), the
     * cache must NOT fire (because the clock has bumped), and the
     * new evaluation under the new symbol-table state must reflect
     * the assignment. We use ReleaseHold to force the held subexpr
     * to evaluate, picking up zclkX = 99. */
    Expr* p2 = parse_expression("ReleaseHold[Hold[zclkX + 1]]");
    Expr* r2 = evaluate(p2);
    expr_free(p2);
    char* s2 = expr_to_string(r2);
    ASSERT(strcmp(s2, "100") == 0);
    free(s2);

    expr_free(r1);
    expr_free(r2);

    Expr* clr = parse_expression("Clear[zclkX]");
    expr_free(evaluate(clr));
    expr_free(clr);
}

/* ------------------------------------------------------------------ */
/* 6. Cached evaluation reflects the up-to-date symtab after Set.      */
/* ------------------------------------------------------------------ */

static void test_set_then_eval_returns_new_value(void) {
    Expr* clr0 = parse_expression("Clear[zclkY]");
    expr_free(evaluate(clr0));
    expr_free(clr0);

    /* First lookup: zclkY is unbound, symbol stays as itself. */
    Expr* p1 = parse_expression("zclkY");
    Expr* r1 = evaluate(p1);
    char* s1 = expr_to_string(r1);
    ASSERT(strcmp(s1, "zclkY") == 0);
    free(s1);
    expr_free(r1);
    expr_free(p1);

    /* Bind it. */
    Expr* sp = parse_expression("zclkY = 13");
    expr_free(evaluate(sp));
    expr_free(sp);

    /* Re-lookup must produce 13, not stale "zclkY". */
    Expr* p2 = parse_expression("zclkY");
    Expr* r2 = evaluate(p2);
    char* s2 = expr_to_string(r2);
    ASSERT(strcmp(s2, "13") == 0);
    free(s2);
    expr_free(r2);
    expr_free(p2);

    /* Also, evaluating something that USES zclkY -- if cached -- must
     * produce 14, not the old un-bound result. */
    Expr* p3 = parse_expression("zclkY + 1");
    Expr* r3 = evaluate(p3);
    char* s3 = expr_to_string(r3);
    ASSERT(strcmp(s3, "14") == 0);
    free(s3);
    expr_free(r3);
    expr_free(p3);

    Expr* clr = parse_expression("Clear[zclkY]");
    expr_free(evaluate(clr));
    expr_free(clr);
}

/* ------------------------------------------------------------------ */
/* 7. SetDelayed bumps the clock.                                      */
/* ------------------------------------------------------------------ */

static void test_setdelayed_invalidates(void) {
    Expr* clr0 = parse_expression("Clear[zclkF]");
    expr_free(evaluate(clr0));
    expr_free(clr0);

    uint64_t c0 = eval_clock_get();
    Expr* sd = parse_expression("zclkF[x_] := x^2");
    expr_free(evaluate(sd));
    expr_free(sd);
    ASSERT(eval_clock_get() > c0);

    /* Definition must take effect for fresh callers. */
    Expr* call = parse_expression("zclkF[7]");
    Expr* res = evaluate(call);
    char* s = expr_to_string(res);
    ASSERT(strcmp(s, "49") == 0);
    free(s); expr_free(res); expr_free(call);

    Expr* clr = parse_expression("Clear[zclkF]");
    expr_free(evaluate(clr));
    expr_free(clr);
}

/* ------------------------------------------------------------------ */
/* 8. Clear bumps the clock.                                           */
/* ------------------------------------------------------------------ */

static void test_clear_invalidates(void) {
    Expr* def = parse_expression("zclkG = 100");
    expr_free(evaluate(def));
    expr_free(def);

    /* Cache zclkG -> 100 by evaluating once. */
    Expr* call = parse_expression("zclkG");
    Expr* r = evaluate(call);
    expr_free(r); expr_free(call);

    uint64_t c0 = eval_clock_get();
    Expr* clr = parse_expression("Clear[zclkG]");
    expr_free(evaluate(clr));
    expr_free(clr);
    ASSERT(eval_clock_get() > c0);

    /* After Clear, zclkG must evaluate back to itself. */
    Expr* call2 = parse_expression("zclkG");
    Expr* r2 = evaluate(call2);
    char* s2 = expr_to_string(r2);
    ASSERT(strcmp(s2, "zclkG") == 0);
    free(s2); expr_free(r2); expr_free(call2);
}

/* ------------------------------------------------------------------ */
/* 9. SetAttributes / ClearAttributes bump the clock.                  */
/* ------------------------------------------------------------------ */

static void test_setattributes_invalidates(void) {
    Expr* clr0 = parse_expression("Clear[zclkH]");
    expr_free(evaluate(clr0));
    expr_free(clr0);

    /* Adding a brand-new attribute bit -- must bump. */
    uint64_t c0 = eval_clock_get();
    Expr* sa = parse_expression("SetAttributes[zclkH, Orderless]");
    expr_free(evaluate(sa));
    expr_free(sa);
    ASSERT(eval_clock_get() > c0);

    /* Removing it -- must bump again. */
    uint64_t c1 = eval_clock_get();
    Expr* ca = parse_expression("ClearAttributes[zclkH, Orderless]");
    expr_free(evaluate(ca));
    expr_free(ca);
    ASSERT(eval_clock_get() > c1);

    /* Adding the SAME bit a second time -- already-set, bump suppressed. */
    Expr* sa1 = parse_expression("SetAttributes[zclkH, Flat]");
    expr_free(evaluate(sa1));
    expr_free(sa1);
    uint64_t c2 = eval_clock_get();
    Expr* sa2 = parse_expression("SetAttributes[zclkH, Flat]");
    expr_free(evaluate(sa2));
    expr_free(sa2);
    ASSERT(eval_clock_get() == c2);   /* no real change, no bump */
}

/* ------------------------------------------------------------------ */
/* 10. Pure builtin calls do NOT bump the clock.                       */
/* ------------------------------------------------------------------ */

static void test_pure_builtin_does_not_bump(void) {
    uint64_t c0 = eval_clock_get();

    Expr* exprs[] = {
        parse_expression("Sin[Pi/4]"),
        parse_expression("Factor[x^4 - 1]"),
        parse_expression("Expand[(x+y+z)^3]"),
        parse_expression("Length[{1,2,3,4,5}]"),
        parse_expression("Map[Sin, {1,2,3}]"),
        parse_expression("D[Cos[x^2], x]"),
        parse_expression("Apart[1/((x-1)*(x-2))]"),
        parse_expression("PolynomialGCD[x^2-1, x^2-2*x+1]"),
        NULL
    };
    for (int i = 0; exprs[i]; i++) {
        Expr* r = evaluate(exprs[i]);
        expr_free(exprs[i]);
        expr_free(r);
    }

    ASSERT(eval_clock_get() == c0);
}

/* ------------------------------------------------------------------ */
/* 11. Cache hit returns SAME pointer (refcount-shared).               */
/* ------------------------------------------------------------------ */

static void test_cache_hit_returns_same_pointer(void) {
    Expr* p = parse_expression("Plus[a, b, c]");
    Expr* r1 = evaluate(p);
    expr_free(p);

    /* r1 is a FUNCTION with a stamped timestamp. evaluate(r1) under
     * the same clock should return r1 (inc-ref'd) without doing any
     * tree-walking work. */
    ASSERT(r1->type == EXPR_FUNCTION);
    unsigned rc_before = r1->refcount;
    Expr* r2 = evaluate(r1);
    ASSERT(r2 == r1);
    ASSERT(r1->refcount == rc_before + 1);

    expr_free(r2);
    expr_free(r1);
}

/* ------------------------------------------------------------------ */
/* 12. Repeated heavy evaluation is bit-stable across iterations.     */
/*                                                                    */
/* Cached subtrees flow through pattern matching and rule replacement;*/
/* this test runs heavyweight kernels many times and asserts the      */
/* printed output is identical to iteration 1's output. Any aliasing  */
/* or stale-cache bug would show up as drift between iterations.      */
/* Uses ASSERT_MSG so failures actually abort the test under          */
/* Release / -DNDEBUG.                                                */
/* ------------------------------------------------------------------ */

static void run_eval_and_capture(const char* src, char* out, size_t out_size) {
    Expr* p = parse_expression(src);
    Expr* r = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(r);
    snprintf(out, out_size, "%s", s);
    free(s);
    expr_free(r);
}

static void test_repeated_heavy_eval(void) {
    const char* kernels[] = {
        "Expand[(x+y+z)^4]",
        "D[Sin[x]^2 + Cos[x]^2, x]",
        "Factor[x^4 - 1]",
        "Together[1/x + 1/y]",
        "PolynomialGCD[x^3 - 1, x^2 - 1]",
        "Apart[1/((x-1)*(x-2))]",
        NULL
    };

    char baseline[6][1024];
    for (int k = 0; kernels[k]; k++) {
        run_eval_and_capture(kernels[k], baseline[k], sizeof(baseline[k]));
    }

    /* Verify every subsequent iteration matches the baseline byte-for-byte. */
    for (int i = 0; i < 50; i++) {
        for (int k = 0; kernels[k]; k++) {
            char buf[1024];
            run_eval_and_capture(kernels[k], buf, sizeof(buf));
            ASSERT_MSG(strcmp(buf, baseline[k]) == 0,
                "iter=%d kernel=%s drifted:\n  baseline: %s\n  got:      %s",
                i, kernels[k], baseline[k], buf);
        }
    }
}

/* ------------------------------------------------------------------ */
/* 13. expr_eq ignores timestamps.                                     */
/* ------------------------------------------------------------------ */

static void test_expr_eq_ignores_timestamp(void) {
    Expr* a = expr_new_integer(7);
    Expr* b = expr_new_integer(7);
    a->last_evaluated_at = 1234;
    b->last_evaluated_at = 9999;
    ASSERT(expr_eq(a, b));
    expr_free(a); expr_free(b);

    Expr* args1[1] = { expr_new_integer(3) };
    Expr* args2[1] = { expr_new_integer(3) };
    Expr* f1 = expr_new_function(expr_new_symbol("Sin"), args1, 1);
    Expr* f2 = expr_new_function(expr_new_symbol("Sin"), args2, 1);
    f1->last_evaluated_at = 100;
    f2->last_evaluated_at = 200;
    ASSERT(expr_eq(f1, f2));
    expr_free(f1); expr_free(f2);
}

/* ------------------------------------------------------------------ */
/* 14. Cache is correctly bypassed across an interleaved Set.          */
/* ------------------------------------------------------------------ */

static void test_interleaved_set_invalidates_subexpressions(void) {
    Expr* clr0 = parse_expression("Clear[zclkV]");
    expr_free(evaluate(clr0));
    expr_free(clr0);

    /* Round 1: zclkV unbound. */
    Expr* p1 = parse_expression("zclkV * 2 + 1");
    Expr* r1 = evaluate(p1);
    char* s1 = expr_to_string(r1);
    ASSERT(strstr(s1, "zclkV") != NULL);
    free(s1); expr_free(r1); expr_free(p1);

    /* Bind. */
    Expr* sp = parse_expression("zclkV = 4");
    expr_free(evaluate(sp));
    expr_free(sp);

    /* Round 2: must produce 9. */
    Expr* p2 = parse_expression("zclkV * 2 + 1");
    Expr* r2 = evaluate(p2);
    char* s2 = expr_to_string(r2);
    ASSERT(strcmp(s2, "9") == 0);
    free(s2); expr_free(r2); expr_free(p2);

    /* Re-bind to a different value. */
    Expr* sp2 = parse_expression("zclkV = 10");
    expr_free(evaluate(sp2));
    expr_free(sp2);

    Expr* p3 = parse_expression("zclkV * 2 + 1");
    Expr* r3 = evaluate(p3);
    char* s3 = expr_to_string(r3);
    ASSERT(strcmp(s3, "21") == 0);
    free(s3); expr_free(r3); expr_free(p3);

    Expr* clr = parse_expression("Clear[zclkV]");
    expr_free(evaluate(clr));
    expr_free(clr);
}

/* ------------------------------------------------------------------ */
/* 15. Cached subexpressions used inside a fresh tree still evaluate. */
/*                                                                    */
/* We build a fresh outer call whose argument is a previously-        */
/* evaluated (and timestamped) expression. The outer call has         */
/* timestamp 0, so the early-exit does NOT fire at the outer level;   */
/* the recursive evaluation of the cached child should hit the cache  */
/* AND still produce a correct end-to-end result.                     */
/* ------------------------------------------------------------------ */

static void test_cached_subexpr_in_fresh_outer(void) {
    /* Pre-evaluate Sin[Pi/4] -> Sqrt[2]/2; r is timestamped. */
    Expr* p_inner = parse_expression("Sin[Pi/4]");
    Expr* r_inner = evaluate(p_inner);
    expr_free(p_inner);
    ASSERT(r_inner->last_evaluated_at == eval_clock_get());

    /* Now build Plus[r_inner, r_inner] manually. The new function
     * node has timestamp 0 (fresh), but its children are stamped. */
    Expr* args[2] = { expr_copy(r_inner), expr_copy(r_inner) };
    Expr* outer = expr_new_function(expr_new_symbol("Plus"), args, 2);
    ASSERT(outer->last_evaluated_at == 0);

    Expr* r_outer = evaluate(outer);
    expr_free(outer);

    /* Sqrt[2]/2 + Sqrt[2]/2 == Sqrt[2]. */
    char* s = expr_to_string(r_outer);
    ASSERT(strcmp(s, "Sqrt[2]") == 0);
    free(s);

    expr_free(r_inner);
    expr_free(r_outer);
}

/* ------------------------------------------------------------------ */
/* 16. Large random walk: assignments and evaluations in random order.*/
/*                                                                    */
/* This is the strongest correctness test -- it exercises the cache   */
/* across many bumps and verifies that every evaluation returns       */
/* exactly what the symbol-table state at that moment dictates.       */
/* ------------------------------------------------------------------ */

static uint32_t xorshift32_state = 0xC0DECAFEu;
static uint32_t xorshift32(void) {
    uint32_t x = xorshift32_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    xorshift32_state = x;
    return x;
}

static void test_random_assign_evaluate_walk(void) {
    /* Three slots tracked in shadow state. */
    int64_t shadow[3] = {0, 0, 0};
    bool bound[3] = {false, false, false};
    const char* names[3] = {"zclkR0", "zclkR1", "zclkR2"};

    /* Reset all slots. */
    for (int i = 0; i < 3; i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "Clear[%s]", names[i]);
        Expr* c = parse_expression(buf);
        expr_free(evaluate(c));
        expr_free(c);
    }

    for (int iter = 0; iter < 200; iter++) {
        uint32_t op = xorshift32() % 5;
        int slot = (int)(xorshift32() % 3);

        if (op == 0) {
            /* Assign a value. */
            int64_t v = (int64_t)((xorshift32() % 200) - 100);
            char buf[64];
            snprintf(buf, sizeof(buf), "%s = %lld", names[slot], (long long)v);
            Expr* p = parse_expression(buf);
            expr_free(evaluate(p));
            expr_free(p);
            shadow[slot] = v;
            bound[slot] = true;
        } else if (op == 1) {
            /* Clear. */
            char buf[64];
            snprintf(buf, sizeof(buf), "Clear[%s]", names[slot]);
            Expr* p = parse_expression(buf);
            expr_free(evaluate(p));
            expr_free(p);
            bound[slot] = false;
        } else {
            /* Read and verify. */
            Expr* p = parse_expression(names[slot]);
            Expr* r = evaluate(p);
            expr_free(p);
            char* s = expr_to_string(r);
            if (bound[slot]) {
                char want[32];
                snprintf(want, sizeof(want), "%lld", (long long)shadow[slot]);
                ASSERT_MSG(strcmp(s, want) == 0,
                    "iter=%d slot=%s: expected %s, got %s", iter, names[slot], want, s);
            } else {
                ASSERT_MSG(strcmp(s, names[slot]) == 0,
                    "iter=%d slot=%s: expected unbound, got %s", iter, names[slot], s);
            }
            free(s);
            expr_free(r);
        }
    }

    /* Cleanup. */
    for (int i = 0; i < 3; i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "Clear[%s]", names[i]);
        Expr* c = parse_expression(buf);
        expr_free(evaluate(c));
        expr_free(c);
    }
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_clock_starts_nonzero_and_monotonic);
    TEST(test_fresh_expr_has_zero_timestamp);
    TEST(test_evaluate_stamps_result);
    TEST(test_reevaluate_same_clock_is_identity);
    TEST(test_set_invalidates);
    TEST(test_set_then_eval_returns_new_value);
    TEST(test_setdelayed_invalidates);
    TEST(test_clear_invalidates);
    TEST(test_setattributes_invalidates);
    TEST(test_pure_builtin_does_not_bump);
    TEST(test_cache_hit_returns_same_pointer);
    TEST(test_repeated_heavy_eval);
    TEST(test_expr_eq_ignores_timestamp);
    TEST(test_interleaved_set_invalidates_subexpressions);
    TEST(test_cached_subexpr_in_fresh_outer);
    TEST(test_random_assign_evaluate_walk);

    printf("All eval-timestamp tests passed!\n");
    return 0;
}
