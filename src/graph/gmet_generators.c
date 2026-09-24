/* gmet_generators.c - named graph families.
 *
 *   WheelGraph[n]              hub 1 joined to the cycle 2..n (n = 1 or n >= 4)
 *   HypercubeGraph[n]          2^n vertices; i ~ j when i-1, j-1 differ in one bit
 *   GridGraph[{n1, ..., nk}]   k-dimensional grid; the FIRST coordinate varies
 *                              fastest (vertex 1 + x1 + n1 x2 + n1 n2 x3 ...)
 *   KaryTree[n] / [n, k]       n vertices in heap order, k children each (k = 2)
 *   CompleteKaryTree[n] / [n, k]  the complete k-ary tree with n levels
 *   CirculantGraph[n, j] / [n, {j1, ...}]   i ~ i +- j (mod n)
 *   PetersenGraph[] / [n, k]   generalized Petersen graph (default 5, 2): inner
 *                              star 1..n (i ~ i + k mod n), outer cycle
 *                              n+1..2n, spokes i ~ n + i
 *   TuranGraph[n, k]           complete k-partite graph with parts as equal as
 *                              possible, larger parts first
 *   HararyGraph[k, n]          the minimal k-connected graph on n vertices
 *   CompleteGraph[{n1, ...}]   complete multipartite graph (CompleteGraph[n] is
 *                              delegated to the existing builtin_complete_graph)
 *
 * Every family is undirected on the vertices 1..N, and -- matching
 * Mathematica 15's EdgeList for each of them, checked with wolframscript --
 * the edge list is the set of pairs {i, j}, i < j, in lexicographic order.
 * Arguments Wolfram rejects (non-positive sizes, a Petersen step k = 0 mod n,
 * the multigraphs WheelGraph[2|3] would be) leave the call unevaluated.
 * Options (DirectedEdges, VertexCoordinates, ...) are not supported.
 *
 * Large graphs are cheap: edges are emitted directly in sorted order where the
 * family allows it (grids, hypercubes, trees, wheels) and otherwise sorted as
 * packed 64-bit keys; every edge shares one UndirectedEdge head node and the
 * vertex Integer nodes by reference.
 */

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#define GEN_MAX_VERTICES 100000000LL
#define GEN_MAX_EDGES    200000000LL

/* ---- pair accumulator ------------------------------------------------------ */

typedef struct { uint64_t* k; int64_t len, cap; int sorted; } Pairs;

static int pairs_add(Pairs* p, int64_t a, int64_t b) {     /* 0-based, a != b */
    if (a > b) { int64_t t = a; a = b; b = t; }
    if (p->len == p->cap) {
        int64_t nc = p->cap ? p->cap * 2 : 64;
        uint64_t* nk = realloc(p->k, (size_t)nc * sizeof(uint64_t));
        if (!nk) return 0;
        p->k = nk; p->cap = nc;
    }
    uint64_t key = ((uint64_t)a << 32) | (uint64_t)b;
    if (p->len > 0 && p->k[p->len - 1] >= key) p->sorted = 0;
    p->k[p->len++] = key;
    return 1;
}

static int u64_cmp(const void* x, const void* y) {
    uint64_t a = *(const uint64_t*)x, b = *(const uint64_t*)y;
    return (a > b) - (a < b);
}

/* Builds Graph[Range[n], {sorted, deduplicated UndirectedEdges}] and frees p. */
static Expr* pairs_graph(int64_t n, Pairs* p) {
    if (!p->sorted && p->len > 1) qsort(p->k, (size_t)p->len, sizeof(uint64_t), u64_cmp);
    int64_t m = 0;
    for (int64_t i = 0; i < p->len; i++)
        if (i == 0 || p->k[i] != p->k[i - 1]) p->k[m++] = p->k[i];
    Expr** vs = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    Expr** es = malloc((size_t)(m > 0 ? m : 1) * sizeof(Expr*));
    if (!vs || !es) { free(vs); free(es); free(p->k); return NULL; }
    for (int64_t i = 0; i < n; i++) vs[i] = expr_new_integer(i + 1);
    Expr* head = expr_new_symbol(SYM_UndirectedEdge);
    for (int64_t i = 0; i < m; i++) {
        Expr* ab[2] = { expr_copy(vs[p->k[i] >> 32]), expr_copy(vs[p->k[i] & 0xffffffffu]) };
        es[i] = expr_new_function(expr_copy(head), ab, 2);
    }
    expr_free(head);
    free(p->k);
    Expr* ga[2] = { expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)n),
                    expr_new_function(expr_new_symbol(SYM_List), es, (size_t)m) };
    free(vs); free(es);
    return expr_new_function(expr_new_symbol(SYM_Graph), ga, 2);
}

/* Positive machine integer argument (>= lo), or -1. */
static int64_t int_arg(const Expr* e, int64_t lo) {
    if (!e || e->type != EXPR_INTEGER || e->data.integer < lo) return -1;
    return e->data.integer;
}

/* ---- families -------------------------------------------------------------- */

Expr* builtin_wheel_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int64_t n = int_arg(res->data.function.args[0], 1);
    if (n < 0 || n > GEN_MAX_VERTICES || n == 2 || n == 3) return NULL;
    Pairs p = { NULL, 0, 0, 1 };
    for (int64_t j = 1; j < n; j++) if (!pairs_add(&p, 0, j)) goto fail;
    if (n >= 4) {
        if (!pairs_add(&p, 1, 2) || !pairs_add(&p, 1, n - 1)) goto fail;
        for (int64_t i = 2; i < n - 1; i++) if (!pairs_add(&p, i, i + 1)) goto fail;
    }
    return pairs_graph(n, &p);
fail:
    free(p.k);
    return NULL;
}

Expr* builtin_hypercube_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int64_t d = int_arg(res->data.function.args[0], 0);
    if (d < 0 || d > 24) return NULL;
    int64_t n = (int64_t)1 << d;
    Pairs p = { NULL, 0, 0, 1 };
    for (int64_t v = 0; v < n; v++)
        for (int b = 0; b < d; b++)
            if (!(v & ((int64_t)1 << b)) && !pairs_add(&p, v, v | ((int64_t)1 << b))) {
                free(p.k); return NULL;
            }
    return pairs_graph(n, &p);
}

Expr* builtin_grid_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* spec = res->data.function.args[0];
    int64_t dims[32];
    int k;
    if (spec->type == EXPR_INTEGER) { dims[0] = int_arg(spec, 1); k = 1; if (dims[0] < 0) return NULL; }
    else {
        if (!graph_is_list(spec)) return NULL;
        k = (int)spec->data.function.arg_count;
        if (k < 1 || k > 32) return NULL;
        for (int i = 0; i < k; i++) {
            dims[i] = int_arg(spec->data.function.args[i], 1);
            if (dims[i] < 0) return NULL;
        }
    }
    int64_t n = 1, stride[32];
    for (int i = 0; i < k; i++) {
        stride[i] = n;
        if (dims[i] > GEN_MAX_VERTICES / n) return NULL;
        n *= dims[i];
    }
    Pairs p = { NULL, 0, 0, 1 };
    for (int64_t v = 0; v < n; v++) {
        int64_t r = v;
        for (int i = 0; i < k; i++) {
            int64_t x = r % dims[i];
            r /= dims[i];
            if (x + 1 < dims[i] && !pairs_add(&p, v, v + stride[i])) { free(p.k); return NULL; }
        }
    }
    return pairs_graph(n, &p);
}

static Expr* kary_tree(int64_t n, int64_t k) {
    Pairs p = { NULL, 0, 0, 1 };
    for (int64_t i = 0; i < n; i++)
        for (int64_t c = 1; c <= k; c++) {
            int64_t ch = k * i + c;
            if (ch >= n) break;
            if (!pairs_add(&p, i, ch)) { free(p.k); return NULL; }
        }
    return pairs_graph(n, &p);
}

Expr* builtin_kary_tree(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    int64_t n = int_arg(res->data.function.args[0], 1);
    int64_t k = argc == 2 ? int_arg(res->data.function.args[1], 1) : 2;
    if (n < 0 || k < 0 || n > GEN_MAX_VERTICES) return NULL;
    return kary_tree(n, k);
}

Expr* builtin_complete_kary_tree(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    int64_t lv = int_arg(res->data.function.args[0], 1);
    int64_t k = argc == 2 ? int_arg(res->data.function.args[1], 1) : 2;
    if (lv < 0 || k < 0) return NULL;
    int64_t n = 0, level = 1;
    for (int64_t i = 0; i < lv; i++) {
        n += level;
        if (n > GEN_MAX_VERTICES) return NULL;
        if (level > GEN_MAX_VERTICES / k) { if (i + 1 < lv) return NULL; }
        else level *= k;
    }
    return kary_tree(n, k);
}

/* Circulant edges i ~ i + s (mod n) for each offset in offs (already reduced). */
static int circulant_add(Pairs* p, int64_t n, const int64_t* offs, int noffs) {
    for (int64_t i = 0; i < n; i++)
        for (int j = 0; j < noffs; j++) {
            int64_t s = offs[j];
            if (s == 0) continue;
            if (!pairs_add(p, i, (i + s) % n)) return 0;
        }
    return 1;
}

Expr* builtin_circulant_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    int64_t n = int_arg(res->data.function.args[0], 1);
    if (n < 0 || n > GEN_MAX_VERTICES) return NULL;
    const Expr* js = res->data.function.args[1];
    int nj;
    const Expr* const* items;
    const Expr* single[1];
    if (js->type == EXPR_INTEGER) { single[0] = js; items = single; nj = 1; }
    else if (graph_is_list(js)) {
        nj = (int)js->data.function.arg_count;
        items = (const Expr* const*)js->data.function.args;
    } else return NULL;
    int64_t* offs = malloc((size_t)(nj > 0 ? nj : 1) * sizeof(int64_t));
    if (!offs) return NULL;
    for (int j = 0; j < nj; j++) {
        if (items[j]->type != EXPR_INTEGER) { free(offs); return NULL; }
        int64_t s = items[j]->data.integer % n;
        if (s < 0) s += n;
        offs[j] = s;
    }
    Pairs p = { NULL, 0, 0, 1 };
    int ok = circulant_add(&p, n, offs, nj);
    free(offs);
    if (!ok) { free(p.k); return NULL; }
    return pairs_graph(n, &p);
}

Expr* builtin_petersen_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    int64_t n = 5, k = 2;
    if (argc == 2) {
        n = int_arg(res->data.function.args[0], 1);
        k = int_arg(res->data.function.args[1], 1);
    } else if (argc != 0) return NULL;
    if (n < 0 || k < 0 || n > GEN_MAX_VERTICES / 2 || k % n == 0) return NULL;
    Pairs p = { NULL, 0, 0, 1 };
    int ok = 1;
    for (int64_t i = 0; i < n && ok; i++) {
        ok = pairs_add(&p, i, (i + k) % n)                    /* inner star  */
          && pairs_add(&p, n + i, n + (i + 1) % n)            /* outer cycle */
          && pairs_add(&p, i, n + i);                         /* spoke       */
    }
    if (!ok) { free(p.k); return NULL; }
    /* n = 2 doubles the inner and outer edges; pairs_graph deduplicates.
     * (n = 1 cannot reach here: k % 1 == 0.) */
    return pairs_graph(2 * n, &p);
}

/* Complete multipartite graph on consecutive parts of the given sizes. */
static Expr* multipartite(const int64_t* sz, int np) {
    int64_t n = 0;
    for (int i = 0; i < np; i++) { n += sz[i]; if (n > GEN_MAX_VERTICES) return NULL; }
    int64_t* part = malloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));
    int64_t* end = malloc((size_t)(np > 0 ? np : 1) * sizeof(int64_t));
    if (!part || !end) { free(part); free(end); return NULL; }
    int64_t v = 0;
    for (int i = 0; i < np; i++) { for (int64_t j = 0; j < sz[i]; j++) part[v++] = i; end[i] = v; }
    Pairs p = { NULL, 0, 0, 1 };
    int ok = 1;
    for (int64_t a = 0; a < n && ok; a++)
        for (int64_t b = end[part[a]]; b < n && ok; b++) ok = pairs_add(&p, a, b);
    free(part); free(end);
    if (!ok) { free(p.k); return NULL; }
    return pairs_graph(n, &p);
}

Expr* builtin_turan_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    int64_t n = int_arg(res->data.function.args[0], 1);
    int64_t k = int_arg(res->data.function.args[1], 1);
    if (n < 0 || k < 0 || n > GEN_MAX_VERTICES) return NULL;
    if (k > n) k = n;
    int64_t* sz = malloc((size_t)k * sizeof(int64_t));
    if (!sz) return NULL;
    for (int64_t i = 0; i < k; i++) sz[i] = n / k + (i < n % k ? 1 : 0);
    Expr* out = multipartite(sz, (int)k);
    free(sz);
    return out;
}

Expr* builtin_gmet_complete_graph(Expr* res) {
    if (res->data.function.arg_count != 1 || !graph_is_list(res->data.function.args[0]))
        return builtin_complete_graph(res);
    const Expr* spec = res->data.function.args[0];
    int np = (int)spec->data.function.arg_count;
    if (np < 1) return NULL;
    int64_t* sz = malloc((size_t)np * sizeof(int64_t));
    if (!sz) return NULL;
    for (int i = 0; i < np; i++) {
        sz[i] = int_arg(spec->data.function.args[i], 1);
        if (sz[i] < 0) { free(sz); return NULL; }
    }
    if (np == 1) {
        /* CompleteGraph[{n}] is K_n (Wolfram), not n isolated vertices. */
        int64_t nv = sz[0];
        free(sz);
        if (nv > GEN_MAX_VERTICES / 64) return NULL;
        sz = malloc((size_t)nv * sizeof(int64_t));
        if (!sz) return NULL;
        for (int64_t i = 0; i < nv; i++) sz[i] = 1;
        np = (int)nv;
    }
    Expr* out = multipartite(sz, np);
    free(sz);
    return out;
}

Expr* builtin_harary_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    int64_t k = int_arg(res->data.function.args[0], 2);
    int64_t n = int_arg(res->data.function.args[1], 1);
    if (k < 0 || n < 0 || n <= k || n > GEN_MAX_VERTICES) return NULL;
    int64_t r = k / 2;
    int64_t* offs = malloc((size_t)(r + 1) * sizeof(int64_t));
    if (!offs) return NULL;
    for (int64_t i = 0; i < r; i++) offs[i] = i + 1;
    int no = (int)r;
    if (k % 2 == 1 && n % 2 == 0) offs[no++] = n / 2;       /* diameters */
    Pairs p = { NULL, 0, 0, 1 };
    int ok = circulant_add(&p, n, offs, no);
    free(offs);
    if (ok && k % 2 == 1 && n % 2 == 1) {
        int64_t h = (n - 1) / 2;
        ok = pairs_add(&p, 0, h) && pairs_add(&p, 0, h + 1);
        for (int64_t i = 1; i < h && ok; i++) ok = pairs_add(&p, i, i + h + 1);
    }
    if (!ok) { free(p.k); return NULL; }
    return pairs_graph(n, &p);
}
