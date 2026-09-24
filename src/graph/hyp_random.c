/* hyp_random.c - RandomHypergraph.
 *
 *   RandomHypergraph[{n, m}, k]
 *       a k-uniform hypergraph on the vertices 1..n with m hyperedges, each an
 *       independent, uniformly random k-subset of 1..n, listed in increasing
 *       order. Hyperedges are drawn independently, so a hyperedge may repeat
 *       (a multi-hypergraph; for n^k >> m repeats are rare). 0 <= k <= n.
 *   RandomHypergraph[{n, m}, k, c]
 *       a List of c such hypergraphs.
 *   RandomHypergraph[{n, {e, a}}]
 *   RandomHypergraph[{n, {{e1, a1}, {e2, a2}, ...}}]
 *       the Function Repository form: e_i hyperedges of arity a_i whose
 *       vertices are drawn uniformly WITH replacement from 1..n (so a vertex
 *       may repeat inside a hyperedge, as in Wolfram-model states). The vertex
 *       list is the vertices that occur, in first-appearance order ("at most n
 *       nodes"). Where the FR function returns the bare edge List, this returns
 *       a Hypergraph; EdgeList recovers the List.
 *
 * Randomness comes from the user-visible stream (random_uniform_01), so every
 * form is reproducible under SeedRandom. Subsets use Floyd's algorithm: k draws
 * per hyperedge, no rejection, O(k) time with a stamp array for membership.
 *
 * Memory (SPEC section 4): returns a fresh Hypergraph (or List); res untouched.
 */

#include "graph_hyper.h"
#include "graph.h"
#include "expr.h"
#include "random.h"
#include "sym_names.h"
#include <stdint.h>
#include <stdlib.h>

static long count_arg(const Expr* e) {
    if (!e || e->type != EXPR_INTEGER || e->data.integer < 0) return -1;
    return (long)e->data.integer;
}

/* Uniform integer in [0, bound). */
static long draw(long bound) {
    long r = (long)(random_uniform_01() * (double)bound);
    return r < bound ? r : bound - 1;
}

static int cmp_long(const void* a, const void* b) {
    long x = *(const long*)a, y = *(const long*)b;
    return (x > y) - (x < y);
}

static Expr* uniform_one(long n, long m, long k, Expr** ints, Expr* lh, long* stamp, long* buf) {
    Expr** es = malloc((size_t)(m > 0 ? m : 1) * sizeof(Expr*));
    Expr** row = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr*));
    Expr** vs = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    if (!es || !row || !vs) { free(es); free(row); free(vs); return NULL; }
    for (long j = 0; j < m; j++) {
        /* Floyd: for t = n-k .. n-1, pick r in [0, t]; take r unless already
         * taken, in which case take t (which cannot have been taken yet). */
        long c = 0;
        for (long t = n - k; t < n; t++) {
            long r = draw(t + 1);
            long pick = (stamp[r] == j) ? t : r;
            stamp[pick] = j;
            buf[c++] = pick;
        }
        if (k <= 16) {
            for (long a = 1; a < k; a++) {             /* insertion sort */
                long x = buf[a], b = a - 1;
                while (b >= 0 && buf[b] > x) { buf[b + 1] = buf[b]; b--; }
                buf[b + 1] = x;
            }
        } else {
            qsort(buf, (size_t)k, sizeof(long), cmp_long);
        }
        for (long a = 0; a < k; a++) row[a] = expr_copy(ints[buf[a]]);
        es[j] = expr_new_function(expr_copy(lh), row, (size_t)k);
    }
    /* Reset the stamps for the next hypergraph (they are compared to j). */
    for (long i = 0; i < n; i++) stamp[i] = -1;
    for (long i = 0; i < n; i++) vs[i] = expr_copy(ints[i]);
    Expr* args[2] = { expr_new_function(expr_copy(lh), vs, (size_t)n),
                      expr_new_function(expr_copy(lh), es, (size_t)m) };
    free(es); free(row); free(vs);
    return expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), args, 2);
}

/* FR form: spec is {e, a} or a List of such pairs. */
static Expr* fr_form(long n, const Expr* spec) {
    if (!graph_is_list(spec)) return NULL;
    size_t np = 1;
    int nested = spec->data.function.arg_count > 0 && graph_is_list(spec->data.function.args[0]);
    if (nested) np = spec->data.function.arg_count;
    long total = 0;
    for (size_t p = 0; p < np; p++) {
        const Expr* pr = nested ? spec->data.function.args[p] : spec;
        if (!graph_is_list(pr) || pr->data.function.arg_count != 2) return NULL;
        long e = count_arg(pr->data.function.args[0]);
        long a = count_arg(pr->data.function.args[1]);
        if (e < 0 || a < 0 || (a > 0 && n < 1)) return NULL;
        if (e > INT32_MAX - total) return NULL;
        total += e;
    }
    Expr** ints = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    Expr** es = malloc((size_t)(total > 0 ? total : 1) * sizeof(Expr*));
    if (!ints || !es) { free(ints); free(es); return NULL; }
    for (long i = 0; i < n; i++) ints[i] = expr_new_integer(i + 1);
    Expr* lh = expr_new_symbol(SYM_List);
    long c = 0;
    for (size_t p = 0; p < np; p++) {
        const Expr* pr = nested ? spec->data.function.args[p] : spec;
        long e = (long)pr->data.function.args[0]->data.integer;
        long a = (long)pr->data.function.args[1]->data.integer;
        Expr** row = malloc((size_t)(a > 0 ? a : 1) * sizeof(Expr*));
        if (!row) break;
        for (long j = 0; j < e; j++) {
            for (long t = 0; t < a; t++) row[t] = expr_copy(ints[draw(n)]);
            es[c++] = expr_new_function(expr_copy(lh), row, (size_t)a);
        }
        free(row);
    }
    for (long i = 0; i < n; i++) expr_free(ints[i]);
    free(ints);
    Expr* el = expr_new_function(lh, es, (size_t)c);
    free(es);
    /* Hypergraph[edges]: the constructor derives the occurring vertices. */
    return expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), &el, 1);
}

Expr* builtin_random_hypergraph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;
    const Expr* spec = res->data.function.args[0];
    if (!graph_is_list(spec) || spec->data.function.arg_count != 2) return NULL;
    long n = count_arg(spec->data.function.args[0]);
    if (n < 0 || n > INT32_MAX) return NULL;
    if (argc == 1) return fr_form(n, spec->data.function.args[1]);

    long m = count_arg(spec->data.function.args[1]);
    long k = count_arg(res->data.function.args[1]);
    if (m < 0 || m > INT32_MAX || k < 0 || k > n) return NULL;
    long reps = 1;
    if (argc == 3) {
        reps = count_arg(res->data.function.args[2]);
        if (reps < 0 || reps > INT32_MAX) return NULL;
    }
    Expr** ints = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    long* stamp = malloc((size_t)(n > 0 ? n : 1) * sizeof(long));
    long* buf = malloc((size_t)(k > 0 ? k : 1) * sizeof(long));
    Expr** outs = malloc((size_t)(reps > 0 ? reps : 1) * sizeof(Expr*));
    if (!ints || !stamp || !buf || !outs) {
        free(ints); free(stamp); free(buf); free(outs); return NULL;
    }
    for (long i = 0; i < n; i++) { ints[i] = expr_new_integer(i + 1); stamp[i] = -1; }
    Expr* lh = expr_new_symbol(SYM_List);
    long made = 0;
    for (; made < reps; made++) {
        outs[made] = uniform_one(n, m, k, ints, lh, stamp, buf);
        if (!outs[made]) break;
    }
    for (long i = 0; i < n; i++) expr_free(ints[i]);
    free(ints); free(stamp); free(buf);
    Expr* out = NULL;
    if (made == reps) {
        if (argc == 2) { out = outs[0]; }
        else out = expr_new_function(expr_copy(lh), outs, (size_t)reps);
    } else {
        for (long i = 0; i < made; i++) expr_free(outs[i]);
    }
    expr_free(lh);
    free(outs);
    return out;
}
