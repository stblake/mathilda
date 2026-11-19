/* Randomised differential stress test of the LinkedList storage core
 * (src/linkedlist.c) against a trivially-correct reference model: a plain
 * growable array of owned Expr*.
 *
 * Every operation is applied to both; after each one the lengths must agree,
 * a random index is compared, and periodically the whole list (ll_to_list), a
 * shallow copy (ll_copy) and the payload's memory bound are checked. The
 * element generator moves through regimes (machine integers, reals, complex,
 * then mixed and symbolic), so typed storage, the one-way typed->boxed
 * conversion, retyping of an emptied list, front/back runs, chunk growth,
 * tail trimming and the queue rebase are all driven through many states.
 *
 * Deterministic (fixed LCG seed per round). Usage: linkedlist_stress_tests [ops]
 * (default 2,000,000 operations over 8 rounds). Labelled nightly by name.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "linkedlist.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------- model */

typedef struct { Expr** v; int64_t n, cap; } Model;

static void m_reserve(Model* m, int64_t need) {
    if (need <= m->cap) return;
    int64_t c = m->cap ? m->cap : 16;
    while (c < need) c *= 2;
    m->v = realloc(m->v, (size_t)c * sizeof(Expr*));
    if (!m->v) abort();
    m->cap = c;
}
static void m_push_back(Model* m, Expr* e) { m_reserve(m, m->n + 1); m->v[m->n++] = e; }
static void m_push_front(Model* m, Expr* e) {
    m_reserve(m, m->n + 1);
    memmove(m->v + 1, m->v, (size_t)m->n * sizeof(Expr*));
    m->v[0] = e;
    m->n++;
}
static Expr* m_pop_back(Model* m)  { return m->v[--m->n]; }
static Expr* m_pop_front(Model* m) {
    Expr* e = m->v[0];
    memmove(m->v, m->v + 1, (size_t)(m->n - 1) * sizeof(Expr*));
    m->n--;
    return e;
}
static void m_clear(Model* m) {
    for (int64_t i = 0; i < m->n; i++) expr_free(m->v[i]);
    m->n = 0;
}

/* ------------------------------------------------------------- rng */

static uint64_t g_rng;
static uint64_t rnd(void) {
    g_rng = g_rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
    return g_rng >> 17;
}
static int64_t rnd_below(int64_t n) { return (int64_t)(rnd() % (uint64_t)n); }

/* Element regimes. */
enum { R_INT, R_REAL, R_CPLX, R_MIXED, R_NREG };

static Expr* gen_elem(int regime) {
    int pick = regime == R_MIXED ? (int)rnd_below(6) : regime;
    switch (pick) {
        case R_INT:  return expr_new_integer((int64_t)rnd() - (int64_t)(UINT64_C(1) << 45));
        case R_REAL: return expr_new_real((double)(int64_t)rnd() / 1024.0);
        case R_CPLX: {
            Expr* a[2] = { expr_new_real((double)rnd_below(1000) + 0.5),
                           expr_new_real(-(double)rnd_below(1000) - 0.25) };
            return expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
        }
        case 3: return expr_new_symbol(intern_symbol(rnd_below(2) ? "x" : "y"));
        case 4: {
            Expr* a[2] = { expr_new_integer(rnd_below(9)), expr_new_symbol(intern_symbol("z")) };
            return expr_new_function(expr_new_symbol(SYM_List), a, 2);
        }
        default: return expr_new_string("s");
    }
}

/* ------------------------------------------------------------- checks */

static void fail(const char* what, int64_t step) {
    fprintf(stderr, "FAIL at step %lld: %s\n", (long long)step, what);
    abort();
}

static void check_index(const LinkedList* ll, const Model* m, int64_t i, int64_t step) {
    Expr* g = ll_get(ll, i);
    if (!expr_eq(g, m->v[i])) {
        char* a = expr_to_string(g);
        char* b = expr_to_string(m->v[i]);
        fprintf(stderr, "index %lld: list %s, model %s\n", (long long)i, a, b);
        free(a); free(b);
        fail("element mismatch", step);
    }
    expr_free(g);
}

static void check_all(const LinkedList* ll, const Model* m, int64_t step) {
    if (ll_length(ll) != m->n) fail("length mismatch", step);
    /* Through the bulk path (packed when typed) ... */
    Expr* l = ll_to_list(ll);
    Expr* ml = expr_new_function(expr_new_symbol(SYM_List), m->v, (size_t)m->n);
    /* expr_new_function adopts the pointers: copy them first. */
    for (int64_t i = 0; i < m->n; i++) expr_copy(m->v[i]);
    if (!expr_eq(l, ml)) fail("ll_to_list mismatch", step);
    expr_free(ml);
    expr_free(l);
    /* ... and through a shallow copy, which must be independent. */
    LinkedList* c = ll_copy(ll);
    if (ll_length(c) != m->n || ll_kind(c) != ll_kind(ll)) {
        fprintf(stderr, "copy: length %lld vs %lld, kind %d vs %d\n",
                (long long)ll_length(c), (long long)m->n, (int)ll_kind(c), (int)ll_kind(ll));
        fail("copy shape", step);
    }
    for (int64_t k = 0; k < 8 && m->n; k++) check_index(c, m, rnd_below(m->n), step);
    ll_release(c);
    /* Memory stays O(live): allocated chunk slots are bounded by a small
     * multiple of the live count plus a few base chunks (rebase + trim). */
    int64_t bound = (int64_t)sizeof(void*) * 0 + 64 * 1024 + 48 * (m->n + 64);
    if (ll_payload_bytes(ll) > bound) {
        fprintf(stderr, "payload %lld bytes for %lld live\n",
                (long long)ll_payload_bytes(ll), (long long)m->n);
        fail("memory not O(live)", step);
    }
}

/* ------------------------------------------------------------- driver */

static void round_run(uint64_t seed, int64_t ops) {
    g_rng = seed;
    LinkedList* ll = ll_new();
    Model m = { NULL, 0, 0 };
    int regime = R_INT;
    /* Bias toward growth, then shrinkage, then queue traffic, cyclically. */
    for (int64_t step = 0; step < ops; step++) {
        if (step % 50000 == 0) regime = (int)rnd_below(R_NREG);
        int phase = (int)((step / 20000) % 3);        /* 0 grow, 1 shrink, 2 queue */
        int64_t r = rnd_below(100);
        if (phase == 2) {
            /* queue: push back, pop front */
            if (r < 50 || m.n == 0) {
                Expr* e = gen_elem(regime);
                m_push_back(&m, expr_copy(e));
                ll_push_back(ll, e);
            } else {
                Expr* a = ll_pop_front(ll);
                Expr* b = m_pop_front(&m);
                if (!expr_eq(a, b)) fail("queue pop_front mismatch", step);
                expr_free(a); expr_free(b);
            }
        } else if (r < (phase == 0 ? 35 : 15)) {
            Expr* e = gen_elem(regime);
            m_push_back(&m, expr_copy(e));
            ll_push_back(ll, e);
        } else if (r < (phase == 0 ? 60 : 25)) {
            Expr* e = gen_elem(regime);
            m_push_front(&m, expr_copy(e));
            ll_push_front(ll, e);
        } else if (r < 65) {
            /* typed fast pushes */
            int64_t v = (int64_t)rnd_below(1000);
            if (rnd_below(2)) { ll_push_back_i64(ll, v); m_push_back(&m, expr_new_integer(v)); }
            else { ll_push_front_real(ll, (double)v); m_push_front(&m, expr_new_real((double)v)); }
        } else if (r < 80) {
            if (m.n) {
                bool back = rnd_below(2);
                Expr* a = back ? ll_pop_back(ll) : ll_pop_front(ll);
                Expr* b = back ? m_pop_back(&m) : m_pop_front(&m);
                if (!expr_eq(a, b)) fail("pop mismatch", step);
                expr_free(a); expr_free(b);
            } else if (ll_pop_back(ll) != NULL || ll_pop_front(ll) != NULL) {
                fail("pop of empty not NULL", step);
            }
        } else if (r < 90) {
            if (m.n) {
                int64_t i = rnd_below(m.n);
                Expr* e = gen_elem(regime);
                expr_free(m.v[i]);
                m.v[i] = expr_copy(e);
                ll_set(ll, i, e);
            }
        } else if (r < 92) {
            /* bulk append of a small List of the regime's elements */
            int64_t k = rnd_below(40);
            Expr** a = malloc((size_t)(k ? k : 1) * sizeof(Expr*));
            for (int64_t t = 0; t < k; t++) { a[t] = gen_elem(regime); m_push_back(&m, expr_copy(a[t])); }
            Expr* lst = expr_new_function(expr_new_symbol(SYM_List), a, (size_t)k);
            free(a);
            if (!ll_append_elements(ll, lst)) fail("append_elements declined a List", step);
            expr_free(lst);
        } else if (r < 93 && phase == 1) {
            /* drain completely: exercises retyping of an emptied list */
            while (m.n) {
                Expr* a = ll_pop_back(ll);
                Expr* b = m_pop_back(&m);
                if (!expr_eq(a, b)) fail("drain mismatch", step);
                expr_free(a); expr_free(b);
            }
        } else if (m.n) {
            check_index(ll, &m, rnd_below(m.n), step);
        }
        if (ll_length(ll) != m.n) fail("length mismatch", step);
        if ((step & 4095) == 0) check_all(ll, &m, step);
    }
    check_all(ll, &m, ops);
    ll_release(ll);
    m_clear(&m);
    free(m.v);
}

/* A 10^6-deep chain of nested lists: the cycle check, ByteCount-style
 * traversal and release must all run in constant C stack. */
static int64_t elem_bytes(Expr* e) { (void)e; return 0; }

static void deep_chain(int64_t depth) {
    LinkedList* root = ll_new();
    LinkedList* cur = root;
    for (int64_t i = 0; i < depth; i++) {
        LinkedList* nx = ll_new();
        ll_push_back(cur, expr_new_linkedlist(ll_retain(nx)));
        ll_release(cur == root ? NULL : cur);    /* keep only the chain's refs */
        cur = nx;
    }
    Expr* back = expr_new_linkedlist(ll_retain(root));
    if (!ll_would_cycle(cur, back)) fail("deep cycle not detected", depth);
    expr_free(back);
    if (ll_byte_count(root, elem_bytes) <= depth) fail("deep byte count", depth);
    if (cur != root) ll_release(cur);
    ll_release(root);                            /* frees the whole chain */
}

int main(int argc, char** argv) {
    int64_t ops = argc > 1 ? atoll(argv[1]) : 2000000;
    symtab_init();
    core_init();
    const int rounds = 8;
    for (int r = 0; r < rounds; r++) {
        printf("round %d: %lld ops\n", r, (long long)(ops / rounds));
        fflush(stdout);
        round_run(UINT64_C(0x9E3779B97F4A7C15) * (uint64_t)(r + 1), ops / rounds);
    }
    printf("deep chain: 10^6\n");
    deep_chain(1000000);
    printf("All LinkedList stress tests passed!\n");
    symtab_clear();
    return 0;
}
