/* test_eval_eager_exit.c -- M3 phase-4 (§3.4) tests.
 *
 * Phase-4 invariants:
 *   - evaluate_step now takes a `bool* changed` out-parameter and sets
 *     it true iff a real rewrite fired during the step. The outer
 *     fixed-point loop uses this signal to skip the O(tree) expr_eq
 *     compare on the common case where nothing fires.
 *
 *   - For atoms (Integer, Real, String, BigInt), evaluate_step never
 *     reports a change.
 *
 *   - For symbols with no OwnValue, evaluate_step never reports a
 *     change.
 *
 *   - For function nodes whose head/args are already at a fixed point
 *     and to which no DownValue / built-in / special primitive applies,
 *     evaluate_step never reports a change.
 *
 *   - When a built-in (e.g. Plus, Times) unconditionally rebuilds its
 *     output even though no terms combined, the change flag is a
 *     false-positive but the expr_eq fallback still recognises the
 *     fixed point and the loop terminates in one iteration. This is
 *     the correctness-preserving guarantee from §3.4: false positives
 *     cost at most one extra expr_eq, never an infinite loop.
 *
 *   - Combining §3.3 (timestamps) and §3.4: a re-evaluation of a
 *     stamped expression hits the early-exit at the top of evaluate()
 *     and never enters evaluate_step at all.
 *
 *   - Heavy mixed workloads still produce identical, byte-stable
 *     output across many repeated evaluations.
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
#include <stdbool.h>

/* Re-declare here because test_utils.h does not pull in eval.h. */
extern uint64_t eval_clock_get(void);
extern Expr* evaluate_step(Expr* e, bool* changed);

/* ------------------------------------------------------------------ */
/* 1. Atoms never report a change.                                    */
/* ------------------------------------------------------------------ */

static void test_evaluate_step_atom_no_change(void) {
    /* Build a few atoms and verify evaluate_step leaves *changed=false. */
    Expr* atoms[] = {
        parse_expression("42"),
        parse_expression("3.14"),
        parse_expression("\"hello\""),
        parse_expression("123456789012345678901234567890"),  /* BigInt */
    };
    for (size_t i = 0; i < sizeof(atoms) / sizeof(atoms[0]); i++) {
        bool changed = true;  /* poison so we observe a write */
        Expr* result = evaluate_step(atoms[i], &changed);
        ASSERT_MSG(!changed, "atom %zu reported a change spuriously", i);
        expr_free(result);
        expr_free(atoms[i]);
    }
}

/* ------------------------------------------------------------------ */
/* 2. Bare symbol with no OwnValue never reports a change.            */
/* ------------------------------------------------------------------ */

static void test_evaluate_step_bare_symbol_no_change(void) {
    Expr* x = parse_expression("eagerSym1");
    bool changed = true;
    Expr* result = evaluate_step(x, &changed);
    ASSERT_MSG(!changed, "bare symbol reported a change spuriously");
    expr_free(result);
    expr_free(x);
}

/* ------------------------------------------------------------------ */
/* 3. Bare symbol WITH an OwnValue reports a change.                  */
/* ------------------------------------------------------------------ */

static void test_evaluate_step_bound_symbol_changes(void) {
    Expr* set = parse_expression("eagerSym2 = 99");
    expr_free(evaluate(set));
    expr_free(set);

    Expr* x = parse_expression("eagerSym2");
    bool changed = false;
    Expr* result = evaluate_step(x, &changed);
    ASSERT_MSG(changed, "bound symbol failed to report a change");
    expr_free(result);
    expr_free(x);

    expr_free(evaluate(parse_expression("Clear[eagerSym2]")));
}

/* ------------------------------------------------------------------ */
/* 4. Built-in that fires reports a change.                           */
/* ------------------------------------------------------------------ */

static void test_evaluate_step_builtin_change(void) {
    /* Sin[Pi/4] has a known closed form; Sin[] should fire. */
    Expr* e = parse_expression("Sin[Pi/4]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "Sin[Pi/4] failed to report a change");
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 5. NULL out-param is allowed.                                      */
/* ------------------------------------------------------------------ */

static void test_evaluate_step_null_changed_param(void) {
    Expr* e = parse_expression("1+2");
    Expr* result = evaluate_step(e, NULL);
    ASSERT(result != NULL);
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 6. Plus-of-symbols converges (regression test for the false-positive */
/*    found during §3.4 bring-up: builtin_plus rebuilds the tree even   */
/*    when no terms combined; without the expr_eq fallback the outer    */
/*    loop would never terminate on `a + b + c`).                       */
/* ------------------------------------------------------------------ */

static void test_plus_of_symbols_converges(void) {
    /* The whole point of the test is that this returns within
     * MAX_ITERATIONS without printing $IterationLimit. */
    assert_eval_eq("eagerA + eagerB + eagerC", "eagerA + eagerB + eagerC", 0);
    assert_eval_eq("eagerA * eagerB * eagerC", "eagerA eagerB eagerC", 0);
    assert_eval_eq("eagerA - eagerB + eagerA", "2 eagerA - eagerB", 0);
}

/* ------------------------------------------------------------------ */
/* 7. Re-evaluation of an already-evaluated expression hits the §3.3   */
/*    early-exit at the top of evaluate() and short-circuits without   */
/*    invoking evaluate_step at all (so the change-flag instrumentation*/
/*    is irrelevant on the cache hit).                                 */
/* ------------------------------------------------------------------ */

static void test_repeated_evaluate_uses_timestamp_path(void) {
    Expr* e = parse_expression("Sin[x] + Cos[x] + Tan[x]");
    Expr* first = evaluate(e);
    expr_free(e);

    /* Re-evaluating the result must (a) be a no-op semantically,
     * (b) leave the timestamp untouched and equal to the live clock. */
    uint64_t clock_before = eval_clock_get();
    uint64_t stamp = first->last_evaluated_at;
    ASSERT_MSG(stamp == clock_before, "first evaluation did not stamp");

    Expr* second = evaluate(first);
    ASSERT_MSG(second->last_evaluated_at == clock_before,
               "timestamp moved on a no-op re-evaluation");

    /* Output must be byte-stable. */
    char* a = expr_to_string(first);
    char* b = expr_to_string(second);
    ASSERT_STR_EQ(a, b);
    free(a); free(b);

    expr_free(first); expr_free(second);
}

/* ------------------------------------------------------------------ */
/* 8. Listable threading reports a change.                            */
/* ------------------------------------------------------------------ */

static void test_listable_change(void) {
    Expr* e = parse_expression("Sin[{x, y}]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "Listable threading failed to report a change");
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 9. Flat flattening reports a change.                               */
/* ------------------------------------------------------------------ */

static void test_flat_change(void) {
    /* Plus[a, Plus[b, c]] -- evaluate sub-Plus first, but for evaluate_step
     * the explicit nested Plus must be flattened by the FLAT pass. */
    Expr* e = parse_expression("Plus[eagerP, Plus[eagerQ, eagerR]]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "Flat flattening failed to report a change");
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 10. Sequence flattening reports a change.                          */
/* ------------------------------------------------------------------ */

static void test_sequence_flatten_change(void) {
    Expr* e = parse_expression("List[1, Sequence[2, 3], 4]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "Sequence flatten failed to report a change");
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 11. Orderless sorting reports a change when args are out of order. */
/* ------------------------------------------------------------------ */

static void test_orderless_unsorted_change(void) {
    /* Plus[c, b, a] -- argument evaluation is a no-op (each is a bare
     * symbol with no OwnValue), so the change flag should come from
     * either Orderless re-sort or the builtin_plus rebuild. Both of
     * those count as "fired"; we just need the flag to be true on
     * the first step. */
    Expr* e = parse_expression("Plus[eagerZ, eagerY, eagerX]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "out-of-order Orderless args failed to fire");
    expr_free(result);
    expr_free(e);
}

/* ------------------------------------------------------------------ */
/* 12. DownValue match reports a change.                              */
/* ------------------------------------------------------------------ */

static void test_downvalue_change(void) {
    Expr* def = parse_expression("eagerF[x_] := x^2");
    expr_free(evaluate(def));
    expr_free(def);

    Expr* e = parse_expression("eagerF[7]");
    bool changed = false;
    Expr* result = evaluate_step(e, &changed);
    ASSERT_MSG(changed, "DownValue match failed to report a change");
    expr_free(result);
    expr_free(e);

    expr_free(evaluate(parse_expression("Clear[eagerF]")));
}

/* ------------------------------------------------------------------ */
/* 13. Heavy idempotence sweep: many random expressions must produce  */
/*    byte-stable output across 50 repeated evaluations.              */
/* ------------------------------------------------------------------ */

static void test_repeated_heavy_eval(void) {
    const char* kernels[] = {
        "Expand[(x + 1)^7]",
        "Factor[x^4 - 1]",
        "Sin[Pi/3] + Cos[Pi/4]",
        "Integrate[D[Sin[x]^2, x], x]",
        "Simplify[Sin[x]^2 + Cos[x]^2]",
        "Plus @@ Range[20]",
        "Map[#^2 &, {1, 2, 3, 4, 5}]",
        "Table[i j, {i, 1, 4}, {j, 1, 4}]",
    };
    const size_t nk = sizeof(kernels) / sizeof(kernels[0]);

    char* baselines[8];
    for (size_t k = 0; k < nk; k++) {
        Expr* e = parse_expression(kernels[k]);
        Expr* r = evaluate(e);
        expr_free(e);
        baselines[k] = expr_to_string(r);
        expr_free(r);
    }

    for (int rep = 0; rep < 50; rep++) {
        for (size_t k = 0; k < nk; k++) {
            Expr* e = parse_expression(kernels[k]);
            Expr* r = evaluate(e);
            expr_free(e);
            char* out = expr_to_string(r);
            ASSERT_MSG(strcmp(out, baselines[k]) == 0,
                       "Kernel %zu drifted on rep %d:\n  baseline: %s\n  got:      %s",
                       k, rep, baselines[k], out);
            free(out);
            expr_free(r);
        }
    }

    for (size_t k = 0; k < nk; k++) free(baselines[k]);
}

/* ------------------------------------------------------------------ */
/* 14. Deep nested no-op: a fully-evaluated expression should not   */
/*    drift across many re-evaluations.                              */
/* ------------------------------------------------------------------ */

static void test_deep_nested_noop_stability(void) {
    Expr* e = parse_expression(
        "{Plus[a, b, c], Times[d, e, f], "
        " Plus[a, Times[b, c], Plus[d, e]], "
        " Sin[x] + Cos[y] + Tan[z]}");
    Expr* r = evaluate(e);
    expr_free(e);

    char* baseline = expr_to_string(r);

    for (int i = 0; i < 30; i++) {
        Expr* re = evaluate(r);
        char* s = expr_to_string(re);
        ASSERT_MSG(strcmp(s, baseline) == 0,
                   "Deep nested expression drifted at iter %d:\n  baseline: %s\n  got:      %s",
                   i, baseline, s);
        free(s);
        expr_free(re);
    }

    free(baseline);
    expr_free(r);
}

/* ------------------------------------------------------------------ */
/* 15. Set+evaluate cycle produces correct results across redefinitions.*/
/*    Validates that the change flag interacts correctly with the     */
/*    eval-clock invalidation from §3.3.                              */
/* ------------------------------------------------------------------ */

static void test_set_then_reeval(void) {
    /* Define f[x_] := x + 1, evaluate f[10] -> 11. Redefine f[x_] := x*2,
     * evaluate f[10] -> 20. The clock bump from the second Set must
     * invalidate any cached evaluation. */
    expr_free(evaluate(parse_expression("Clear[eagerG]")));
    expr_free(evaluate(parse_expression("eagerG[x_] := x + 1")));
    assert_eval_eq("eagerG[10]", "11", 0);

    expr_free(evaluate(parse_expression("Clear[eagerG]")));
    expr_free(evaluate(parse_expression("eagerG[x_] := x * 2")));
    assert_eval_eq("eagerG[10]", "20", 0);

    expr_free(evaluate(parse_expression("Clear[eagerG]")));
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_evaluate_step_atom_no_change);
    TEST(test_evaluate_step_bare_symbol_no_change);
    TEST(test_evaluate_step_bound_symbol_changes);
    TEST(test_evaluate_step_builtin_change);
    TEST(test_evaluate_step_null_changed_param);
    TEST(test_plus_of_symbols_converges);
    TEST(test_repeated_evaluate_uses_timestamp_path);
    TEST(test_listable_change);
    TEST(test_flat_change);
    TEST(test_sequence_flatten_change);
    TEST(test_orderless_unsorted_change);
    TEST(test_downvalue_change);
    TEST(test_repeated_heavy_eval);
    TEST(test_deep_nested_noop_stability);
    TEST(test_set_then_reeval);

    printf("\nAll §3.4 eager-exit tests passed!\n");
    return 0;
}
