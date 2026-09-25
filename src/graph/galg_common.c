/* galg_common.c - shared plumbing for the algos stream (see graph_algos.h).
 *
 * Everything here reads a graph through the validated-graph memo
 * (graph_edge_indices), so a builder costs O(V + E) integer work and no
 * expression hashing. The simple undirected CSR (GalgUG) is the working form of
 * every algorithm that ignores edge direction -- matchings, covers, cliques,
 * planarity, undirected cuts -- and is built by a counting sort, so neighbour
 * lists come out sorted ascending with duplicates removed in the same pass.
 *
 * Memory (SPEC section 4): builders return owned structures (free with
 * galg_ug_free); result builders return fresh Exprs that share vertex/edge
 * subexpressions with the source graph by reference (expr_copy).
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include "core.h"
#include "pack.h"
#include "ndarray.h"
#include <stdlib.h>
#include <string.h>

void galg_poll(void) { tc_check_deadline(); }

int galg_nv(const Expr* g) {
    return (int)g->data.function.args[0]->data.function.arg_count;
}

long galg_ne(const Expr* g) {
    return (long)g->data.function.args[1]->data.function.arg_count;
}

void galg_ug_free(GalgUG* u) {
    if (!u) return;
    free(u->off); free(u->adj); free(u);
}

/* Counting-sort build: bucket every half-edge (a, b) by b, then scatter into
 * rows by a. Scanning the b-buckets in increasing b fills each row in
 * increasing neighbour order, so a duplicate (from an anti-parallel pair or a
 * directed edge beside an undirected one) is adjacent to its twin and dropped
 * in the compaction pass. O(n + m), no comparison sort. */
GalgUG* galg_ug_from_edges(int n, long m, const int* eu, const int* ev) {
    size_t nn = (size_t)(n > 0 ? n : 0) + 1;
    size_t hm = (size_t)(2 * (m > 0 ? m : 0)) + 1;
    GalgUG* u = calloc(1, sizeof(GalgUG));
    int* cnt  = calloc(nn, sizeof(int));
    int* bst  = calloc(nn, sizeof(int));
    int* tmpa = malloc(hm * sizeof(int));
    int* tmpb = malloc(hm * sizeof(int));
    if (!u || !cnt || !bst || !tmpa || !tmpb) goto fail;
    u->n = n;
    u->off = calloc(nn, sizeof(int));
    u->adj = malloc(hm * sizeof(int));
    if (!u->off || !u->adj) goto fail;

    /* Stage 1: bucket half-edges by target b. */
    for (long k = 0; k < m; k++) {
        if (eu[k] == ev[k]) continue;
        cnt[ev[k]]++; cnt[eu[k]]++;
    }
    {
        int s = 0;
        for (int i = 0; i < n; i++) { bst[i] = s; s += cnt[i]; }
    }
    for (long k = 0; k < m; k++) {
        int a = eu[k], b = ev[k];
        if (a == b) continue;
        tmpa[bst[b]] = a; tmpb[bst[b]++] = b;
        tmpa[bst[a]] = b; tmpb[bst[a]++] = a;
    }
    /* Stage 2: scatter into rows by source a, in increasing b. */
    {
        int s = 0;
        for (int i = 0; i < n; i++) { u->off[i] = s; s += cnt[i]; }
        u->off[n] = s;
        for (int i = 0; i < n; i++) bst[i] = u->off[i];
        for (int j = 0; j < s; j++) u->adj[bst[tmpa[j]]++] = tmpb[j];
    }
    /* Stage 3: compact each row, dropping adjacent duplicates. */
    {
        int w = 0;
        for (int i = 0; i < n; i++) {
            int lo = u->off[i], hi = u->off[i + 1];
            u->off[i] = w;
            for (int j = lo; j < hi; j++)
                if (j == lo || u->adj[j] != u->adj[j - 1]) u->adj[w++] = u->adj[j];
        }
        u->off[n] = w;
        u->m = w / 2;
    }
    free(cnt); free(bst); free(tmpa); free(tmpb);
    return u;
fail:
    free(cnt); free(bst); free(tmpa); free(tmpb);
    galg_ug_free(u);
    return NULL;
}

GalgUG* galg_ug_from_graph(const Expr* g) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;
    return galg_ug_from_edges(galg_nv(g), galg_ne(g), eu, ev);
}

Expr* galg_vertex_list(const Expr* g, const int* idx, int k, int sort) {
    const Expr* verts = g->data.function.args[0];
    int* tmp = NULL;
    if (sort && k > 1) {
        /* Counting-free insertion sort is quadratic; use a mark array instead:
         * O(n) to emit in vertex order. */
        int n = galg_nv(g);
        char* mark = calloc((size_t)n + 1, 1);
        tmp = malloc((size_t)k * sizeof(int));
        if (!mark || !tmp) { free(mark); free(tmp); return NULL; }
        for (int i = 0; i < k; i++) mark[idx[i]] = 1;
        int w = 0;
        for (int v = 0; v < n && w < k; v++) if (mark[v]) tmp[w++] = v;
        free(mark);
        k = w;
        idx = tmp;
    }
    Expr** items = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr*));
    if (!items) { free(tmp); return NULL; }
    for (int i = 0; i < k; i++) items[i] = expr_copy(verts->data.function.args[idx[i]]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
    free(items); free(tmp);
    return out;
}

Expr* galg_edge_list(const Expr* g, const long* eidx, long k) {
    const Expr* edges = g->data.function.args[1];
    Expr** items = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr*));
    if (!items) return NULL;
    for (long i = 0; i < k; i++) items[i] = expr_copy(edges->data.function.args[eidx[i]]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
    free(items);
    return out;
}

Expr* galg_make_edge(const Expr* g, const Expr* head, int u, int v) {
    const Expr* verts = g->data.function.args[0];
    Expr* a[2];
    a[0] = expr_copy(verts->data.function.args[u]);
    a[1] = expr_copy(verts->data.function.args[v]);
    return expr_new_function(expr_copy((Expr*)head), a, 2);
}

int galg_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

int galg_is_string(const Expr* e, const char* s) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, s) == 0;
}

const Expr* galg_plain_list(const Expr* e, Expr** owned) {
    *owned = NULL;
    if (!e) return NULL;
    if (e->type == EXPR_NDARRAY) {
        if (!is_packed_list(e)) return NULL;
        *owned = pack_unpack(e);
        e = *owned;
    }
    return graph_is_list(e) ? e : NULL;
}

int galg_is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
            || e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

int galg_split_options(const Expr* res, size_t from, const char* const* names,
                       const Expr** values) {
    size_t argc = res->data.function.arg_count;
    for (int i = 0; names[i]; i++) values[i] = NULL;
    for (size_t k = from; k < argc; k++) {
        const Expr* r = res->data.function.args[k];
        if (!galg_is_rule(r)) return 0;
        const Expr* key = r->data.function.args[0];
        int hit = 0;
        for (int i = 0; names[i]; i++) {
            if (galg_is_symbol(key, names[i]) || galg_is_string(key, names[i])) {
                values[i] = r->data.function.args[1];
                hit = 1;
                break;
            }
        }
        if (!hit) return 0;
    }
    return 1;
}

size_t galg_first_option(const Expr* res, size_t from) {
    size_t argc = res->data.function.arg_count;
    size_t k = from;
    while (k < argc && !galg_is_rule(res->data.function.args[k])) k++;
    return k;
}

int galg_vertex_arg(const Expr* g, const Expr* v) {
    int p = graph_vertex_position(g, v);
    return p >= 0 ? p : -1;
}

Expr* galg_truth(int b) { return expr_new_symbol(b ? SYM_True : SYM_False); }
