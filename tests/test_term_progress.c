/*
 * test_term_progress.c -- the loop-progress stack the REPL spinner reads
 * (src/term/term_progress.c) and the iterator hooks that write it: tokens,
 * the slot cap, range counts, the TimeConstrained unwind, nesting, and the
 * worker-thread counter.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include "term/term_progress.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static void eval_discard(const char* src) {
    Expr* p = parse_expression(src);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* A builtin that records what the stack looks like from inside a loop body. */
static int    g_seen_depth = -1;
static MthProgressView g_seen_outer, g_seen_inner;
static Expr* builtin_probe(Expr* res) {
    (void)res;
    g_seen_depth = term_progress_sample(&g_seen_outer, &g_seen_inner);
    return expr_new_symbol(SYM_Null);
}

static void test_tokens(void) {
    ASSERT(g_mth_progress_depth == 0);
    int a = mth_progress_push("Outer", 10);
    int b = mth_progress_push("Inner", 0);
    ASSERT(a == 0 && b == 1 && g_mth_progress_depth == 2);
    mth_progress_step(a, 3);
    mth_progress_step(b, 7);
    MthProgressView o, i;
    ASSERT(term_progress_sample(&o, &i) == 2);
    ASSERT_STR_EQ(o.head, "Outer");
    ASSERT(o.done == 3 && o.total == 10);
    ASSERT_STR_EQ(i.head, "Inner");
    ASSERT(i.done == 7 && i.total == 0);
    ASSERT(i.gen == o.gen + 1);
    /* a missed inner pop is repaired by the outer one */
    mth_progress_pop(a);
    ASSERT(g_mth_progress_depth == 0);
    ASSERT(term_progress_sample(&o, &i) == 0);
    /* negative totals are "unknown" */
    int c = mth_progress_push("X", -5);
    ASSERT(g_mth_progress[c].total == 0);
    mth_progress_pop(c);
}

static void test_cap(void) {
    int tok[MTH_PROGRESS_SLOTS + 4];
    for (int k = 0; k < MTH_PROGRESS_SLOTS + 4; k++) tok[k] = mth_progress_push("L", k + 1);
    ASSERT(g_mth_progress_depth == MTH_PROGRESS_SLOTS + 4);
    mth_progress_step(tok[MTH_PROGRESS_SLOTS + 2], 99);     /* beyond the cap: ignored */
    MthProgressView o, i;
    ASSERT(term_progress_sample(&o, &i) == MTH_PROGRESS_SLOTS);
    ASSERT(o.total == 1 && i.total == 2);
    mth_progress_pop(tok[0]);
    ASSERT(g_mth_progress_depth == 0);
}

static void test_range_count(void) {
    ASSERT(mth_progress_range_count(1, 10, 1) == 10);
    ASSERT(mth_progress_range_count(1, 10, 2) == 5);
    ASSERT(mth_progress_range_count(0, 1, 0.1) == 11);
    ASSERT(mth_progress_range_count(10, 1, -1) == 10);
    ASSERT(mth_progress_range_count(10, 1, 1) == 0);    /* empty */
    ASSERT(mth_progress_range_count(1, 10, 0) == 0);
    ASSERT(mth_progress_range_count(1, HUGE_VAL, 1) == 0);
    ASSERT(mth_progress_range_count(1, 1e30, 1) == 0);  /* too long for a long */
}

static void test_iterators(void) {
    symtab_add_builtin("ProgressProbe", builtin_probe);
    struct { const char* src; const char* head; long done, total; } cases[] = {
        { "Do[If[i == 3, ProgressProbe[]], {i, 5}]",          "Do",       2, 5 },
        { "Do[If[i == 5, ProgressProbe[]], {i, 1, 7, 2}]",    "Do",       2, 4 },
        { "Do[If[i === c, ProgressProbe[]], {i, {a, b, c}}]", "Do",       2, 3 },
        { "k = 0; Do[k++; If[k == 2, ProgressProbe[]], 6]",   "Do",       1, 6 },
        { "Table[If[i == 5, ProgressProbe[]]; i, {i, 8}]",    "Table",    4, 8 },
        { "Table[If[i === b, ProgressProbe[]]; i, {i, {a, b}}]", "Table", 1, 2 },
        { "k = 0; Table[k++; If[k == 3, ProgressProbe[]], 4]", "Table",   2, 4 },
        { "Sum[If[i == 2, ProgressProbe[]]; f[i], {i, 3}]",   "Sum",      1, 3 },
        { "Sum[If[i === y, ProgressProbe[]]; f[i], {i, {x, y}}]", "Sum",  1, 2 },
        { "Product[If[i == 3, ProgressProbe[]]; g[i], {i, 4}]", "Product", 2, 4 },
        { "Nest[(If[# === h[h[x]], ProgressProbe[]]; h[#]) &, x, 5]", "Nest", 2, 5 },
        { "NestList[(If[# === h[x], ProgressProbe[]]; h[#]) &, x, 3]", "NestList", 1, 3 },
        { "FoldList[(If[#2 == 3, ProgressProbe[]]; #1 + #2) &, 0, {1, 2, 3, 4}]", "FoldList", 2, 4 },
    };
    for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        g_seen_depth = -1;
        memset(&g_seen_outer, 0, sizeof g_seen_outer);
        eval_discard(cases[k].src);
        ASSERT_MSG(g_seen_depth == 1, "%s: depth %d", cases[k].src, g_seen_depth);
        ASSERT_MSG(g_seen_outer.head && strcmp(g_seen_outer.head, cases[k].head) == 0,
                   "%s: head %s", cases[k].src, g_seen_outer.head ? g_seen_outer.head : "(null)");
        ASSERT_MSG(g_seen_outer.done == cases[k].done && g_seen_outer.total == cases[k].total,
                   "%s: %ld/%ld, want %ld/%ld", cases[k].src, g_seen_outer.done,
                   g_seen_outer.total, cases[k].done, cases[k].total);
        ASSERT_MSG(g_mth_progress_depth == 0, "%s left depth %d", cases[k].src, g_mth_progress_depth);
    }
    /* nesting: the multi-spec Table is an outer and an inner loop */
    g_seen_depth = -1;
    eval_discard("Table[If[i == 2 && j == 3, ProgressProbe[]]; i j, {i, 3}, {j, 4}]");
    ASSERT(g_seen_depth == 2);
    ASSERT_STR_EQ(g_seen_outer.head, "Table");
    ASSERT(g_seen_outer.done == 1 && g_seen_outer.total == 3);
    ASSERT_STR_EQ(g_seen_inner.head, "Table");
    ASSERT(g_seen_inner.done == 2 && g_seen_inner.total == 4);
    ASSERT(g_mth_progress_depth == 0);
    /* Break, Return, Throw and Abort all pop */
    eval_discard("Do[If[i == 2, Break[]], {i, 10}]");
    eval_discard("f[] := Do[If[i == 2, Return[i]], {i, 10}]; f[]");
    eval_discard("Catch[Table[If[i == 2, Throw[i]]; i, {i, 10}]]");
    eval_discard("CheckAbort[Sum[If[i == 2, Abort[]]; i, {i, 10}], x]");
    ASSERT(g_mth_progress_depth == 0);
}

static void test_timeconstrained(void) {
    term_progress_init();
    eval_discard("TimeConstrained[Table[Do[j, {j, 10^9}], {i, 3}], 0.3, timeout]");
    ASSERT_MSG(g_mth_progress_depth == 0, "a timeout left depth %d", g_mth_progress_depth);
    /* and the stack still works afterwards */
    g_seen_depth = -1;
    eval_discard("Do[If[i == 2, ProgressProbe[]], {i, 2}]");
    ASSERT(g_seen_depth == 1 && g_seen_outer.done == 1);
}

static void test_workers_and_limit(void) {
    ASSERT(mth_progress_workers() == 0);
    mth_progress_workers_add(4);
    mth_progress_workers_add(2);
    ASSERT(mth_progress_workers() == 6);
    mth_progress_workers_add(-6);
    ASSERT(mth_progress_workers() == 0);
    mth_progress_workers_add(-1);           /* never negative */
    ASSERT(mth_progress_workers() == 0);
    /* a parallel driver leaves the count balanced */
    eval_discard("FactorInteger[NextPrime[10^25] NextPrime[10^26]]");
    ASSERT(mth_progress_workers() == 0);
    ASSERT(mth_progress_memory_limit() == 0);
    mth_progress_set_memory_limit(1u << 30);
    ASSERT(mth_progress_memory_limit() == (1u << 30));
    mth_progress_set_memory_limit(0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_tokens);
    TEST(test_cap);
    TEST(test_range_count);
    TEST(test_iterators);
    TEST(test_timeconstrained);
    TEST(test_workers_and_limit);
    printf("All term_progress tests passed.\n");
    return 0;
}
