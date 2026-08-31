/* test_graph_slow.c - RG-2 acceptance rows that are too slow for the default suite.
 *
 * WHY A SEPARATE BINARY. These rows exercise FindVertexColoring's guards against
 * the unbounded case, and the only way to do that honestly is to actually spend
 * the budget -- on the order of a minute or two. Two reasons that cannot live in
 * test_graph.c:
 *
 *   1. test_utils.h arms alarm(60) in a constructor, so a >60s test is killed
 *      mid-run and its buffered output is lost -- which looks like a silent pass
 *      or an inexplicable hang, not a failure. main() below re-arms the alarm.
 *   2. A two-minute test in the default suite gets deleted or commented out by
 *      the third person who trips over it.
 *
 * So this target is EXCLUDE_FROM_ALL and has no add_test(): `make` does not
 * build it and ctest does not run it. Run it deliberately:
 *
 *     cd tests/build && make graph_slow_tests && ./graph_slow_tests
 *
 * The tradeoff is real and worth stating: rows verified here are NOT verified on
 * every build, so a regression in the budget logic will not be caught by CI. The
 * alternative was leaving these as prose in the plan ("verified, numbers
 * documented"), which is unreproducible by the next person -- strictly worse.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "test_utils.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/* Run the exact search over a generator expression, reporting nodes and wall
 * time. Returns fvc_search's verdict: chi, or 0 meaning "gave up". */
static int timed_search(const char* src, long* steps_out, double* secs_out, int* n_out) {
    Expr* sd = evaluate(parse_expression("SeedRandom[1]"));
    expr_free(sd);
    Expr* g = evaluate(parse_expression(src));
    GraphAdj* a = graph_build_adj(g);
    ASSERT(a != NULL);
    int* colour = calloc((size_t)(a->n > 0 ? a->n : 1), sizeof(int));
    long steps = 0;
    clock_t t0 = clock();
    int chi = fvc_search(a, colour, &steps);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    /* When it DID answer, the answer must be a proper colouring -- a budget
     * change must not be able to turn a correct result into a plausible one. */
    if (chi > 0) {
        for (int u = 0; u < a->n; u++) {
            ASSERT_MSG(colour[u] >= 1 && colour[u] <= chi, "colour outside 1..chi");
            for (int j = 0; j < a->outdeg[u]; j++)
                ASSERT_MSG(colour[u] != colour[a->out[u][j]], "adjacent share a colour");
        }
    }
    printf("    %-26s n=%d ret=%d nodes=%ld %.2fs\n", src, a->n, chi, steps, secs);
    fflush(stdout);
    if (steps_out) *steps_out = steps;
    if (secs_out)  *secs_out  = secs;
    if (n_out)     *n_out     = a->n;
    free(colour); graph_adj_free(a); expr_free(g);
    return chi;
}

/* AC-10d: on exhausting FVC_MAX_STEPS the search reports failure (0) rather
 * than handing back its incumbent. The incumbent is a VALID colouring, so the
 * bug this guards against is not a crash -- it is a plausible answer that
 * quietly is not minimal, which is the whole failure mode the exact search
 * exists to prevent. */
static void test_budget_refuses_rather_than_guessing(void) {
    long steps = 0; double secs = 0;
    int chi = timed_search("RandomGraph[{128, 2000}]", &steps, &secs, NULL);
    ASSERT_MSG(chi == 0, "dense n=128 must exhaust the budget and report failure");
    ASSERT_MSG(steps > 8000000L, "must have actually spent the budget, not bailed early");
    /* RNG DEPENDENCE, stated rather than hidden: this row needs an instance hard
     * enough to exhaust 8M nodes, and it identifies one only as
     * `SeedRandom[1]; RandomGraph[{128, 2000}]`. If the RNG stream ever changes,
     * this instance may become solvable and the row will fail even though the
     * budget logic is fine. That failure is loud and self-explaining -- pick a
     * harder instance -- which is why the row is left seed-based rather than
     * carrying 2000 serialised edges. */
}

/* The companion row to AC-10d, and the reason the budget was raised to 8M: a
 * dense n=100 instance is genuinely solvable (3.9M nodes, ~30s) and must
 * ANSWER. At the earlier 2M budget it refused after 14s -- a correct result
 * converted into a refusal by a guard that was sized too tight. */
static void test_budget_does_not_refuse_a_solvable_instance(void) {
    long steps = 0;
    int chi = timed_search("RandomGraph[{100, 1200}]", &steps, NULL, NULL);
    /* Assert the PROPERTY, not the value. chi==8 was what this RNG produced on
     * this host; it has no external ground truth, so pinning it would turn any
     * change to SeedRandom/RandomGraph -- or a platform RNG difference -- into a
     * colouring-test failure with no colouring bug present. The property under
     * test is "the budget did not convert a correct answer into a refusal",
     * which is exactly `answered` plus `inside the budget`. The observed value
     * is printed by timed_search above for anyone comparing runs. */
    ASSERT_MSG(chi > 0, "dense n=100 must ANSWER, not be refused by the budget");
    ASSERT_MSG(steps < 8000000L, "must finish strictly inside the budget");
}

/* Phase 1 manual verification asked for a NUMBER, not "promptly": sparse graphs
 * at the vertex cap. Two shapes, and the measured result is that BOTH meet their
 * bounds (lb == ub) and search nothing -- which is the real finding here, since
 * it shows cost tracks DENSITY rather than n. A sparse graph at n=128 is free;
 * a dense one at the same n is unbounded (see the budget rows below). That is
 * precisely why a vertex cap alone was never a sufficient guard. */
static void test_sparse_at_the_cap_is_fast(void) {
    long steps = 0; double secs = 0;

    int chi = timed_search("CycleGraph[128]", &steps, &secs, NULL);
    ASSERT_MSG(chi == 2, "chi(C128) should be 2");
    ASSERT_MSG(steps == 0, "an even cycle should meet its bounds, searching nothing");
    ASSERT_MSG(secs < 1.0, "sparse at the cap must be well under a second");

    chi = timed_search("RandomGraph[{128, 200}]", &steps, &secs, NULL);
    ASSERT_MSG(chi > 0, "a sparse random graph at the cap must ANSWER, not refuse");
    /* Asserted, not merely observed: the plan ticks a manual-verification box on
     * this instance searching zero nodes, and an unasserted measurement in a
     * commit message is not a criterion anyone can re-check. */
    ASSERT_MSG(steps == 0, "sparse random at the cap should also meet its bounds");
    ASSERT_MSG(secs < 1.0, "sparse at the cap must be well under a second");
}

/* AC-10e: TimeConstrained is the abort channel the docstring points users at,
 * so it has to actually work on this head -- and the head is exactly the shape
 * that breaks it, a single builtin call that never returns to the evaluator's
 * own deadline poll. Two things asserted: the result is $Aborted (not a
 * valid-but-unproven colouring, and not a hang), and the session survives to
 * evaluate afterwards.
 *
 * Lives here rather than in the default suite because it needs an instance slow
 * enough that a 2-second deadline is genuinely hit -- the same dense graph the
 * budget rows use. It returns after ~2s rather than the full ~100s, so it is by
 * far the cheapest row in this file.
 *
 * Not a leak test: aborting via siglongjmp unwinds past the search's frees, a
 * few KB, which is the tree's standing behaviour for every abortable builtin
 * (see ## Risks and Rollback in the plan). This row pins termination and
 * survival only. */
static void test_timeconstrained_aborts_a_slow_instance(void) {
    Expr* sd = evaluate(parse_expression("SeedRandom[1]"));
    expr_free(sd);

    clock_t t0 = clock();
    Expr* r = evaluate(parse_expression(
        "TimeConstrained[FindVertexColoring[RandomGraph[{128, 2000}]], 2]"));
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    char* s = expr_to_string(r);
    printf("    TimeConstrained[...,2] -> %s  %.2fs\n", s, secs);

    ASSERT_MSG(strcmp(s, "$Aborted") == 0,
               "a deadline must yield $Aborted, never an unproven colouring");
    /* The whole point of the deadline: it fires far short of FVC_MAX_STEPS
     * (~100s). Generous margin so a loaded host does not fail spuriously. */
    ASSERT_MSG(secs < 30.0, "the deadline must fire long before the node budget");
    free(s);
    expr_free(r);

    /* The session survives the abort and still evaluates correctly. */
    assert_eval_eq("Max[FindVertexColoring[CycleGraph[5]]]", "3", 0);
}

int main(void) {
    /* test_utils.h's constructor armed alarm(60); these rows are designed to
     * outlast that. Re-arm generously -- 0 would remove the safety net. */
    alarm(600);

    symtab_init();
    core_init();

    TEST(test_sparse_at_the_cap_is_fast);
    TEST(test_timeconstrained_aborts_a_slow_instance);
    TEST(test_budget_does_not_refuse_a_solvable_instance);
    TEST(test_budget_refuses_rather_than_guessing);

    printf("All slow graph tests passed!\n");
    return 0;
}
