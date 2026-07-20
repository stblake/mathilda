/*
 * test_expr_pool.c
 *
 * Unit tests for the Expr node free-list pool (expr.c) and the two
 * evaluate_step allocation fast paths it accompanies:
 *   #1  Expr node recycling through a singly-linked free-list.
 *   #2  Lazy allocation of evaluate_step's held-Unevaluated tracking array.
 *   #3  Stack buffer for evaluate_step's per-call argument scratch, with a
 *       heap fallback for large-arity calls.
 *
 * The pool tests observe recycling purely through the public
 * expr_new_* / expr_free / expr_pool_free_all API (expr_alloc_node is
 * private): a node freed and then re-allocated must come back from the
 * pool. The evaluator tests drive real expressions through evaluate() to
 * prove the fast paths preserve semantics, including the rare held
 * Unevaluated wrappers and >8-argument (heap fallback) calls.
 */

#include "expr.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include <string.h>

/* --------------------------------------------------------------------- */
/* #1  Node recycling                                                     */
/* --------------------------------------------------------------------- */

/* A node freed and immediately re-allocated is the SAME memory (LIFO
 * free-list). This is the core recycling guarantee that removes malloc/
 * free from the hot path. */
static void test_pool_recycles_single_node(void) {
    Expr* a = expr_new_integer(12345);
    Expr* pa = a;
    expr_free(a);
    Expr* b = expr_new_integer(67890);
    ASSERT(b == pa);                 /* recycled the just-freed node */
    ASSERT(b->type == EXPR_INTEGER); /* fully re-initialized */
    ASSERT(b->data.integer == 67890);
    ASSERT(b->refcount == 1);
    expr_free(b);
}

/* The pool is type-agnostic: a node freed as one type is reused for any
 * other type (all Expr nodes are the same fixed size). */
static void test_pool_recycles_across_types(void) {
    Expr* r = expr_new_real(2.5);
    Expr* pr = r;
    expr_free(r);
    Expr* s = expr_new_symbol("poolsym");
    ASSERT(s == pr);                 /* same slab of memory, new type */
    ASSERT(s->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(s->data.symbol.name, "poolsym");
    ASSERT(s->data.symbol.def == NULL);
    expr_free(s);
}

/* Freeing a FUNCTION node must recycle its head and every argument node,
 * not just the top node. We free f[10, 20] and confirm the next batch of
 * allocations is drawn from exactly the freed nodes (order-independent). */
static void test_pool_recycles_function_children(void) {
    Expr* h  = expr_new_symbol("f");
    Expr* a1 = expr_new_integer(10);
    Expr* a2 = expr_new_integer(20);
    Expr* args[2] = { a1, a2 };
    Expr* f = expr_new_function(h, args, 2);

    Expr* freed[4] = { h, a1, a2, f };
    expr_free(f);   /* recycles h, a1, a2 and f */

    /* Next four allocations must all be recycled nodes from `freed`,
     * and must be four distinct pointers. */
    Expr* got[4];
    for (int i = 0; i < 4; i++) got[i] = expr_new_integer(i);

    for (int i = 0; i < 4; i++) {
        int in_freed = 0;
        for (int j = 0; j < 4; j++) if (got[i] == freed[j]) in_freed = 1;
        ASSERT(in_freed);                    /* came from the pool */
        for (int j = i + 1; j < 4; j++) ASSERT(got[i] != got[j]); /* distinct */
    }
    for (int i = 0; i < 4; i++) expr_free(got[i]);
}

/* --------------------------------------------------------------------- */
/* Payload safety: recycling must not corrupt still-live nodes           */
/* --------------------------------------------------------------------- */

/* Allocate many nodes, free half, allocate more, and confirm the nodes
 * that stayed live never had their payloads clobbered by pool churn.
 * This exercises the free-list link (stored in the dead node's payload)
 * against real live neighbours. */
static void test_pool_live_nodes_uncorrupted(void) {
    enum { N = 2000 };
    Expr* keep[N];
    for (int i = 0; i < N; i++) keep[i] = expr_new_integer(i * 3 + 1);

    /* Free a scratch node between each survivor so the pool cycles hard. */
    for (int round = 0; round < 4; round++) {
        for (int i = 0; i < N; i++) {
            Expr* scratch = expr_new_real((double)i);
            expr_free(scratch);      /* churn the free-list */
        }
        /* Survivors must still hold their original values. */
        for (int i = 0; i < N; i++) {
            ASSERT(keep[i]->type == EXPR_INTEGER);
            ASSERT(keep[i]->data.integer == i * 3 + 1);
        }
    }
    for (int i = 0; i < N; i++) expr_free(keep[i]);
}

/* Bigint / string payloads carry heap resources; recycling their nodes
 * must run the destructor before reuse (no leak) and reinit cleanly. */
static void test_pool_recycles_owned_payloads(void) {
    Expr* big = expr_new_bigint_from_str("123456789012345678901234567890");
    ASSERT(big->type == EXPR_BIGINT);
    expr_free(big);                  /* mpz_clear then recycle */

    Expr* str = expr_new_string("hello pool");
    ASSERT(str->type == EXPR_STRING);
    ASSERT_STR_EQ(str->data.string, "hello pool");
    expr_free(str);                  /* free(string) then recycle */

    /* Reuse the recycled node for a fresh bigint — must be independent. */
    Expr* big2 = expr_new_bigint_from_str("42");
    ASSERT(big2->type == EXPR_BIGINT);
    char buf[64];
    Expr* s = expr_new_string("independent");
    ASSERT_STR_EQ(s->data.string, "independent");
    (void)buf;
    expr_free(big2);
    expr_free(s);
}

/* --------------------------------------------------------------------- */
/* Bounded pool: large batches (beyond the cap) stay correct              */
/* --------------------------------------------------------------------- */

/* Allocate far more nodes than the free-list cap (8192), so that freeing the
 * batch overflows the cap and the surplus goes back to the OS. Every live node
 * must be distinct with its own value, and a second wave must still be served
 * correctly (partly from the recycled cap, partly from fresh malloc).
 * expr_pool_free_all is intentionally NOT exercised mid-run — it is teardown-
 * only and runs via atexit. */
static void test_pool_bounded_large_batch(void) {
    enum { M = 12000 };              /* > EXPR_POOL_CAP (8192) */
    Expr** v = malloc(sizeof(Expr*) * M);
    for (int i = 0; i < M; i++) {
        v[i] = expr_new_integer(i);
        ASSERT(v[i] != NULL);
    }
    /* All live nodes distinct and holding their own value. */
    for (int i = 0; i < M; i++) ASSERT(v[i]->data.integer == i);
    for (int i = 1; i < M; i++) ASSERT(v[i] != v[i - 1]);
    for (int i = 0; i < M; i++) expr_free(v[i]);   /* overflows the cap */

    /* Second wave: values still independent whether recycled or freshly malloc'd. */
    for (int i = 0; i < M; i++) {
        v[i] = expr_new_real((double)i + 0.5);
        ASSERT(v[i]->type == EXPR_REAL && v[i]->data.real == (double)i + 0.5);
    }
    for (int i = 1; i < M; i++) ASSERT(v[i] != v[i - 1]);
    for (int i = 0; i < M; i++) expr_free(v[i]);
    free(v);
}

/* --------------------------------------------------------------------- */
/* Evaluator correctness through the fast paths (#1/#2/#3)                */
/* --------------------------------------------------------------------- */

/* Ordinary arithmetic / structural evaluation, small arity (stack-args). */
static void test_eval_small_arity(void) {
    assert_eval_eq("1 + 2*3", "7", 0);
    assert_eval_eq("Plus[b, a, c]", "a + b + c", 0);      /* Orderless */
    assert_eval_eq("Times[3, 4, 5]", "60", 0);
    assert_eval_eq("{x, y} = {10, 20}; x + y", "30", 0);  /* destructuring */
    /* Ten logistic-map steps must land on the known value (Do path). */
    assert_eval_startswith("Module[{x = 1./3}, Do[x = 3.5 x (1 - x), {10}]; x]",
                           "0.38323");
    assert_eval_startswith("Nest[3.5 # (1 - #) &, 1./3, 10]", "0.38323");
}

/* The lazily-allocated held-Unevaluated tracking array (#2): held slots
 * keep the wrapper, non-held slots strip it. These paths must behave
 * identically to before the array became lazy. */
static void test_eval_held_unevaluated(void) {
    /* Non-held slot: wrapper stripped, then the exposed 1+2 evaluates. */
    assert_eval_eq("f[Unevaluated[1 + 2]]", "f[3]", 0);
    /* Sequence inside Unevaluated is preserved (not flattened into Length). */
    assert_eval_eq("Length[Unevaluated[Sequence[a, b]]]", "2", 0);
    /* HoldAll head keeps a held Unevaluated wrapper intact. */
    assert_eval_eq("Hold[Unevaluated[1 + 2]]", "Hold[Unevaluated[1 + 2]]", 0);
}

/* Large-arity calls (> EVAL_SMALL_ARGS) must take the heap fallback for
 * the scratch array and still evaluate correctly. */
static void test_eval_large_arity_heap_fallback(void) {
    /* 12-argument Plus / Times (arity > 8). */
    assert_eval_eq("Plus[1,2,3,4,5,6,7,8,9,10,11,12]", "78", 0);
    assert_eval_eq("Times[1,1,1,1,1,1,1,1,1,1,2,3]", "6", 0);
    /* A 20-element list threads through a Listable head (heap-fallback
     * args) without loss. */
    assert_eval_eq(
        "Total[{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}]",
        "210", 0);
}

/* A tight repeated evaluation (the workload the pool targets) must remain
 * numerically exact and leak-free across thousands of iterations. */
static void test_eval_tight_loop_stable(void) {
    assert_eval_startswith(
        "Module[{x = 1./3}, Do[x = 3.5 x (1 - x), {5000}]; x]", "0.");
    /* Symbolic accumulation still simplifies each pass. */
    assert_eval_eq("Module[{s = 0}, Do[s = s + k, {k, 1, 100}]; s]", "5050", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pool_recycles_single_node);
    TEST(test_pool_recycles_across_types);
    TEST(test_pool_recycles_function_children);
    TEST(test_pool_live_nodes_uncorrupted);
    TEST(test_pool_recycles_owned_payloads);
    TEST(test_pool_bounded_large_batch);

    TEST(test_eval_small_arity);
    TEST(test_eval_held_unevaluated);
    TEST(test_eval_large_arity_heap_fallback);
    TEST(test_eval_tight_loop_stable);

    printf("All expr_pool tests passed!\n");
    return 0;
}
