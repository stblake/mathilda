/* PVec: the persistent vector core.  A randomised model check against a plain
 * array (set / push / pop at sizes straddling every trie boundary), persistence
 * (an old version never changes after an update), structural sharing (an update
 * copies one path, not the vector), typed -> boxed conversion, and the
 * complexity gate: per-update cost must not grow with n. */
#include "expr.h"
#include "core.h"
#include "symtab.h"
#include "pvec.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static uint64_t rng = 88172645463325252ULL;
static uint64_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static int64_t elem_i(PVec* v, int64_t i) { int64_t x; ASSERT_MSG(pv_get_i64(v, i, &x), "typed get"); return x; }

static void check_vs_model(PVec* v, const int64_t* m, int64_t n) {
    ASSERT_MSG(pv_length(v) == n, "length");
    for (int64_t i = 0; i < n; i++) ASSERT_MSG(elem_i(v, i) == m[i], "element mismatch");
}

static void test_model(void) {
    const int64_t CAP = 40000;
    int64_t* m = malloc(sizeof(int64_t) * CAP);
    int64_t n = 0;
    PVec* v = pv_new();
    /* grow past 32, 1056, 32800 (trie levels 1, 2, 3), with sets and pops mixed in */
    for (int step = 0; step < 120000; step++) {
        uint64_t r = rnd() % 100;
        if (r < 60 && n < CAP) { m[n] = (int64_t)(rnd() % 1000); v = pv_push_i64(v, m[n]); n++; }
        else if (r < 80 && n > 0) {
            Expr* out = NULL; v = pv_pop(v, &out);
            ASSERT_MSG(out && out->data.integer == m[n - 1], "pop value"); expr_free(out); n--;
        } else if (n > 0) {
            int64_t i = (int64_t)(rnd() % (uint64_t)n); int64_t x = (int64_t)(rnd() % 1000);
            v = pv_set(v, i, expr_new_integer(x)); m[i] = x;
        }
        if (step % 5000 == 0) check_vs_model(v, m, n);
    }
    check_vs_model(v, m, n);
    /* drain completely */
    while (n > 0) { Expr* out = NULL; v = pv_pop(v, &out); ASSERT_MSG(out->data.integer == m[n - 1], "drain"); expr_free(out); n--; }
    ASSERT_MSG(pv_length(v) == 0, "empty after drain");
    v = pv_push_i64(v, 7);
    ASSERT_MSG(elem_i(v, 0) == 7, "reusable after drain");
    pv_release(v);
    free(m);
}

static void test_boundaries(void) {
    int64_t sizes[] = {0, 1, 31, 32, 33, 63, 64, 65, 1055, 1056, 1057, 1088, 1089, 32799, 32800, 32801, 33824};
    for (size_t s = 0; s < sizeof sizes / sizeof *sizes; s++) {
        int64_t n = sizes[s];
        PVec* v = pv_new();
        for (int64_t i = 0; i < n; i++) v = pv_push_i64(v, i * 3);
        ASSERT_MSG(pv_length(v) == n, "len at boundary");
        for (int64_t i = 0; i < n; i++) ASSERT_MSG(elem_i(v, i) == i * 3, "get at boundary");
        for (int64_t i = n; i > 0; i--) {
            Expr* out = NULL; v = pv_pop(v, &out);
            ASSERT_MSG(out->data.integer == (i - 1) * 3, "pop at boundary"); expr_free(out);
        }
        ASSERT_MSG(pv_length(v) == 0, "drained");
        pv_release(v);
    }
}

static void test_persistence(void) {
    PVec* v = pv_new();
    for (int64_t i = 0; i < 5000; i++) v = pv_push_i64(v, i);
    PVec* snap = pv_retain(v);                      /* keep the old version */
    PVec* w = pv_set(pv_retain(v), 2500, expr_new_integer(-1));
    PVec* x = pv_push_i64(pv_retain(v), 99);
    Expr* out = NULL;
    PVec* y = pv_pop(pv_retain(v), &out); expr_free(out);
    ASSERT_MSG(elem_i(w, 2500) == -1 && elem_i(snap, 2500) == 2500 && elem_i(v, 2500) == 2500, "set is persistent");
    ASSERT_MSG(pv_length(x) == 5001 && pv_length(snap) == 5000 && elem_i(x, 5000) == 99, "push is persistent");
    ASSERT_MSG(pv_length(y) == 4999 && pv_length(snap) == 5000 && elem_i(snap, 4999) == 4999, "pop is persistent");
    for (int64_t i = 0; i < 5000; i++) ASSERT_MSG(elem_i(snap, i) == i, "snapshot intact");
    /* structural sharing: one set copies a path, not the vector */
    int64_t nodes, shared;
    pv_node_stats(w, &nodes, &shared);
    printf("  nodes=%lld shared=%lld\n", (long long)nodes, (long long)shared);
    ASSERT_MSG(shared > nodes / 2, "an update shares most nodes with the original");
    pv_release(v); pv_release(snap); pv_release(w); pv_release(x); pv_release(y);
}

static void test_boxing(void) {
    PVec* v = pv_new();
    for (int64_t i = 0; i < 100; i++) v = pv_push_i64(v, i);
    ASSERT_MSG(pv_kind(v) == PV_I64, "ints stay typed");
    PVec* keep = pv_retain(v);
    v = pv_push(v, expr_new_string("s"));           /* a non-integer converts to boxed once */
    ASSERT_MSG(pv_kind(v) == PV_BOXED && pv_kind(keep) == PV_I64, "boxing leaves the snapshot typed");
    ASSERT_MSG(pv_length(v) == 101 && elem_i(keep, 99) == 99, "snapshot untouched by boxing");
    Expr* e = pv_get(v, 100);
    ASSERT_MSG(e->type == EXPR_STRING, "boxed element");
    expr_free(e);
    e = pv_get(v, 50);
    ASSERT_MSG(e->type == EXPR_INTEGER && e->data.integer == 50, "old elements survive boxing");
    expr_free(e);
    PVec* r = pv_new();
    r = pv_push_real(r, 1.5); r = pv_push_real(r, 2.5);
    r = pv_set(r, 0, expr_new_integer(3));          /* kind change on set */
    ASSERT_MSG(pv_kind(r) == PV_BOXED, "set of another kind boxes");
    pv_release(v); pv_release(keep); pv_release(r);
}

static void test_to_list(void) {
    PVec* v = pv_new();
    for (int64_t i = 0; i < 70; i++) v = pv_push_i64(v, i);
    Expr* l = pv_to_list(v);
    ASSERT_MSG(l->type == EXPR_FUNCTION && l->data.function.arg_count == 70, "to_list length");
    ASSERT_MSG(l->data.function.args[69]->data.integer == 69, "to_list last");
    expr_free(l);
    pv_release(v);
}

/* time n pushes, then n set-at-random-index updates (ns per op) */
static double time_ops(int64_t n, int set_phase) {
    PVec* v = pv_new();
    clock_t t0 = clock();
    for (int64_t i = 0; i < n; i++) v = pv_push_i64(v, i);
    if (set_phase) {
        t0 = clock();
        for (int64_t i = 0; i < n; i++) v = pv_set(v, (int64_t)(rnd() % (uint64_t)n), expr_new_integer(i));
    }
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC * 1e9 / (double)n;
    pv_release(v);
    return dt;
}
/* the persistent form: every update keeps the old version alive (path copy) */
static double time_persistent_sets(int64_t n) {
    PVec* v = pv_new();
    for (int64_t i = 0; i < n; i++) v = pv_push_i64(v, i);
    clock_t t0 = clock();
    PVec* cur = pv_retain(v);
    for (int64_t i = 0; i < n; i++) {
        PVec* nx = pv_set(pv_retain(cur), (int64_t)(rnd() % (uint64_t)n), expr_new_integer(i));
        pv_release(cur); cur = nx;
    }
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC * 1e9 / (double)n;
    pv_release(cur); pv_release(v);
    return dt;
}

static void test_complexity(void) {
    /* ns/op at n and 4n: O(log32 n) means the ratio stays near 1 (one more level
     * costs a few percent); a linear update would be ~4. Generous bound for noise. */
    double push1 = time_ops(200000, 0), push4 = time_ops(800000, 0);
    double set1 = time_ops(200000, 1), set4 = time_ops(800000, 1);
    double ps1 = time_persistent_sets(200000), ps4 = time_persistent_sets(800000);
    printf("  push        n=200k %7.1f ns/op   n=800k %7.1f ns/op   ratio %.2f\n", push1, push4, push4 / push1);
    printf("  set (rc=1)  n=200k %7.1f ns/op   n=800k %7.1f ns/op   ratio %.2f\n", set1, set4, set4 / set1);
    printf("  set (kept)  n=200k %7.1f ns/op   n=800k %7.1f ns/op   ratio %.2f\n", ps1, ps4, ps4 / ps1);
    ASSERT_MSG(push4 / push1 < 2.0, "push is amortised O(1)");
    ASSERT_MSG(set4 / set1 < 2.0, "in-place set is O(log32 n)");
    ASSERT_MSG(ps4 / ps1 < 2.0, "persistent set is O(log32 n), not O(n)");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_model);
    TEST(test_boundaries);
    TEST(test_persistence);
    TEST(test_boxing);
    TEST(test_to_list);
    TEST(test_complexity);
    printf("All pvec tests passed!\n");
    return 0;
}
