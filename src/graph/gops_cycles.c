/* gops_cycles.c - Eulerian cycles, cycles and paths.
 *
 *   FindEulerianCycle[g]       {cycle} -- a closed walk using every edge exactly
 *   FindEulerianCycle[g, 1]    once, as a list of edges -- or {} if g is not
 *                              Eulerian. An edgeless graph with >= 1 vertex gives
 *                              {{}}; the null graph {}. Hierholzer's algorithm,
 *                              iterative, O(V + E): the walk starts at the first
 *                              vertex (VertexList order) with an edge and takes
 *                              edges in EdgeList order. Undirected cycles are
 *                              reported in Hierholzer's pop order, directed ones
 *                              forwards; an undirected edge is written in the
 *                              direction it is walked (Mathematica does the
 *                              same). Mixed graphs, and n > 1 / All, are left
 *                              unevaluated.
 *   FindCycle[g]               {cycle} for some cycle of g, or {}. Linear-time
 *                              DFS in Mathematica's own search order, so the
 *                              cycle reported is Mathematica's.
 *   FindCycle[g, k]            a cycle of length <= k (k may be Infinity)
 *   FindCycle[g, {k}]          a cycle of length exactly k
 *   FindCycle[g, {kmin, kmax}] a cycle with kmin <= length <= kmax
 *   FindCycle[g, kspec, n]     at most n cycles (n a positive integer or All)
 *   FindCycle[{g, v}, ...]     the same, restricted to cycles through vertex v
 *                              (the plain form by BFS from v, linear time)
 *                              Each cycle is found once (not once per rotation
 *                              or, undirected, per direction). The k-forms and
 *                              n > 1 enumerate simple cycles by backtracking
 *                              from each vertex as the cycle's lowest-positioned
 *                              vertex; exponential in the worst case (finding a
 *                              cycle of exact length k is NP-hard), so the
 *                              search polls TimeConstrained and gives up --
 *                              returning the expression unevaluated -- after
 *                              GOPS_SEARCH_MAX_STEPS extension steps. Cycle
 *                              length counts edges; undirected cycles have
 *                              length >= 3, directed ones >= 2. Mixed graphs and
 *                              weighted graphs with a length spec are left
 *                              unevaluated.
 *   FindPath[g, s, t]          {path} for a path from s to t (as a vertex list)
 *                              or {}: the first path depth-first search meets,
 *                              neighbours in EdgeList order -- linear time, and
 *                              the same path Mathematica returns. s == t gives
 *                              {} (as Mathematica).
 *   FindPath[g, s, t, kspec]   a path of length (edge count) within kspec
 *   FindPath[g, s, t, kspec, n]  at most n paths (n or All): the first n simple
 *                              paths found depth-first, reported shortest
 *                              first, as Mathematica. Backtracking, bounded as
 *                              for FindCycle. Weighted graphs with a kspec are
 *                              left unevaluated (Mathematica measures kspec in
 *                              total weight there).
 *
 * Deviation: for the length-bounded / enumerating FindCycle forms, which
 * cycles are reported and in what order is Mathilda's own deterministic choice;
 * Mathematica's differs in general (every reported cycle is a genuine cycle of
 * the requested length).
 *
 * Memory (SPEC section 4): results are fresh; res is borrowed. A
 * TimeConstrained abort (tc_check_deadline) longjmps out of a search, leaking
 * its scratch arrays, as elsewhere in src/graph/.
 */

#include "graph_ops.h"
#include "core.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#define GOPS_SEARCH_MAX_STEPS 50000000L

static Expr* empty_list(void) { return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }

static Expr* wrap1(Expr* x) {
    return expr_new_function(expr_new_symbol(SYM_List), &x, 1);
}

/* Edge node for walking edge e from vertex x to y: the original node when it
 * already reads x -> y (always, for directed edges), else a reoriented copy. */
static Expr* walked_edge(const GopsView* v, GopsHeads* hs, int e, int x, int y) {
    if (v->edir[e] || v->eu[e] == x) return expr_copy(v->edges[e]);
    return gops_edge(hs, 0, expr_copy(v->verts[x]), expr_copy(v->verts[y]));
}

/* ---- FindEulerianCycle ---------------------------------------------------- */
Expr* builtin_find_eulerian_cycle(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    if (argc == 2) {
        const Expr* n = res->data.function.args[1];
        if (n->type != EXPR_INTEGER || n->data.integer != 1) return NULL;
    }
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    int r = gops_eulerian(&v);
    if (r == -1 || r == -2) return NULL;
    if (r == 0) return empty_list();
    if (v.ne == 0) return wrap1(empty_list());

    int directed = (v.ndir == v.ne);
    const GopsInc* ip = gops_inc_cached(&v, directed ? GOPS_INC_OUT : GOPS_INC_ALL);
    if (!ip) return NULL;
    const GopsInc inc = *ip;             /* borrowed arrays, not freed here */
    size_t n = v.nv, ne = v.ne;
    int* ptr = gops_malloc(n, sizeof(int));
    unsigned char* used = gops_calloc(ne, 1);
    int* sv = gops_malloc(ne + 1, sizeof(int));     /* stack: vertex        */
    int* se = gops_malloc(ne + 1, sizeof(int));     /*        entering edge */
    int* cv = gops_malloc(ne + 1, sizeof(int));     /* pop order: vertex    */
    int* ce = gops_malloc(ne + 1, sizeof(int));     /*            edge      */
    Expr** out = gops_malloc(ne, sizeof(Expr*));
    Expr* result = NULL;
    if (!ptr || !used || !sv || !se || !cv || !ce || !out) goto done;
    for (size_t i = 0; i < n; i++) ptr[i] = inc.start[i];
    int s = 0;
    while (inc.start[s] == inc.start[s + 1]) s++;   /* first vertex with an edge */
    size_t top = 0, m = 0;
    sv[top] = s; se[top] = -1; top++;
    while (top > 0) {
        int x = sv[top - 1];
        int* p = &ptr[x];
        while (*p < inc.start[x + 1] && used[inc.eid[*p]]) (*p)++;
        if (*p < inc.start[x + 1]) {
            int e = inc.eid[*p], y = inc.nbr[*p];
            (*p)++;
            used[e] = 1;
            sv[top] = y; se[top] = e; top++;
        } else {
            top--;
            cv[m] = x; ce[m] = se[top]; m++;
        }
    }
    /* m == ne + 1 pops: cv[0..ne] with ce[i] joining cv[i] and cv[i+1]. */
    GopsHeads hs = { { NULL, NULL } };
    if (directed) {
        for (size_t i = 0; i < ne; i++) out[i] = expr_copy(v.edges[ce[ne - 1 - i]]);
    } else {
        for (size_t i = 0; i < ne; i++) out[i] = walked_edge(&v, &hs, ce[i], cv[i], cv[i + 1]);
    }
    gops_heads_free(&hs);
    result = wrap1(gops_list_take(out, ne));
    out = NULL;
done:
    free(ptr); free(used); free(sv); free(se); free(cv); free(ce); free(out);
    return result;
}

/* ---- Length / count specifications ---------------------------------------- */

/* kspec -> [*lo, *hi] (hi = LONG_MAX for Infinity). 0 if malformed. */
static int parse_kspec(const Expr* k, long* lo, long* hi) {
    if (gops_is_infinity(k)) { *lo = 0; *hi = -1; return 1; }
    if (k->type == EXPR_INTEGER) {
        if (k->data.integer < 0) return 0;
        *lo = 0; *hi = (long)k->data.integer; return 1;
    }
    if (!graph_is_list(k)) return 0;
    size_t c = k->data.function.arg_count;
    const Expr* a = c >= 1 ? k->data.function.args[0] : NULL;
    if (c == 1 && a->type == EXPR_INTEGER && a->data.integer >= 0) {
        *lo = *hi = (long)a->data.integer; return 1;
    }
    if (c == 2 && a->type == EXPR_INTEGER && a->data.integer >= 0) {
        const Expr* b = k->data.function.args[1];
        *lo = (long)a->data.integer;
        if (gops_is_infinity(b)) { *hi = -1; return 1; }
        if (b->type == EXPR_INTEGER && b->data.integer >= a->data.integer) {
            *hi = (long)b->data.integer; return 1;
        }
    }
    return 0;
}

/* n -> *count (-1 = All). 0 if malformed. */
static int parse_count(const Expr* n, long* count) {
    if (n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_All) { *count = -1; return 1; }
    if (n->type == EXPR_INTEGER && n->data.integer >= 1) { *count = (long)n->data.integer; return 1; }
    return 0;
}

/* Growable list of result Exprs. */
typedef struct { Expr** a; size_t n, cap; } XList;
static int xl_push(XList* l, Expr* x) {
    if (l->n == l->cap) {
        size_t cap = l->cap ? 2 * l->cap : 8;
        Expr** a = realloc(l->a, cap * sizeof(Expr*));
        if (!a) { expr_free(x); return 0; }
        l->a = a; l->cap = cap;
    }
    l->a[l->n++] = x;
    return 1;
}
static void xl_free(XList* l) {
    for (size_t i = 0; i < l->n; i++) expr_free(l->a[i]);
    free(l->a);
    memset(l, 0, sizeof(*l));
}
static Expr* xl_list(XList* l) {
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), l->a, l->n);
    free(l->a);
    memset(l, 0, sizeof(*l));
    return r;
}

/* The cycle path[0] -> ... -> path[len-1] -> path[0], through edges via[1..len]
 * (via[i] enters path[i]; via[len] closes back to path[0]). */
static Expr* cycle_expr(const GopsView* v, GopsHeads* hs, const int* path,
                        const int* via, long len) {
    Expr** es = gops_malloc((size_t)len, sizeof(Expr*));
    if (!es) return NULL;
    for (long i = 0; i < len; i++) {
        int x = path[i], y = path[(i + 1) % len];
        es[i] = walked_edge(v, hs, via[i + 1], x, y);
    }
    return gops_list_take(es, (size_t)len);
}

/* ---- FindCycle: linear DFS ------------------------------------------------ *
 * Mathematica's search order, reproduced: a stack DFS that scans ALL of a
 * vertex's edges when it is visited -- reporting the first edge back to a
 * vertex on the current root path, pushing the unvisited endpoints -- and then
 * visits the most recently pushed vertex next. It is a genuine DFS (a popped
 * vertex's parent is its latest pusher), so every back edge to an ancestor is
 * seen when its tail is visited, and a graph has a cycle iff one is reported.
 * O(V + E): at most one stack entry per edge. The cycle is written from the
 * ancestor w: directed, forwards along the tree path and back; undirected,
 * across the back edge first and then up the tree path (Mathematica's form). */
static int dfs_any_cycle(const GopsView* v, int directed, XList* out) {
    const GopsInc* ip = gops_inc_cached(v, directed ? GOPS_INC_OUT : GOPS_INC_ALL);
    if (!ip) return 0;
    const GopsInc inc = *ip;             /* borrowed arrays, not freed here */
    size_t n = v->nv, cap = (size_t)inc.start[n] + n + 1;
    int* depth = gops_malloc(n, sizeof(int));       /* -1 = unvisited          */
    int* chain = gops_malloc(n + 1, sizeof(int));   /* chain[d]: path vertex   */
    int* inedge = gops_malloc(n, sizeof(int));
    int* sv = gops_malloc(cap, sizeof(int));        /* pending: vertex         */
    int* se = gops_malloc(cap, sizeof(int));        /*          entering edge  */
    int* sd = gops_malloc(cap, sizeof(int));        /*          depth          */
    int ok = 1;
    if (!depth || !chain || !inedge || !sv || !se || !sd) { ok = 0; goto done; }
    for (size_t i = 0; i < n; i++) depth[i] = -1;
    for (size_t r = 0; r < n && out->n == 0; r++) {
        if (depth[r] >= 0) continue;
        size_t top = 0;
        sv[top] = (int)r; se[top] = -1; sd[top] = 0; top++;
        while (top > 0 && out->n == 0) {
            top--;
            int x = sv[top], d = sd[top];
            if (depth[x] >= 0) continue;               /* stale entry */
            depth[x] = d; chain[d] = x; inedge[x] = se[top];
            for (int j = inc.start[x]; j < inc.start[x + 1]; j++) {
                int y = inc.nbr[j], e = inc.eid[j];
                if (!directed && e == inedge[x]) continue;
                if (depth[y] >= 0 && depth[y] <= d && chain[depth[y]] == y) {
                    /* back edge x -> y, y an ancestor: cycle chain[depth y .. d] */
                    long len = (long)d - depth[y] + 1;
                    int* path = gops_malloc((size_t)len, sizeof(int));
                    int* via = gops_malloc((size_t)len + 1, sizeof(int));
                    if (!path || !via) { free(path); free(via); ok = 0; break; }
                    if (directed) {
                        for (long i = 0; i < len; i++) {
                            path[i] = chain[depth[y] + i];
                            via[i] = i > 0 ? inedge[path[i]] : -1;
                        }
                        via[len] = e;
                    } else {                           /* y, x, parent(x), ... */
                        path[0] = y; via[1] = e;
                        for (long i = 1; i < len; i++) {
                            path[i] = chain[d - (i - 1)];
                            via[i + 1] = inedge[path[i]];
                        }
                    }
                    GopsHeads hs = { { NULL, NULL } };
                    Expr* c = cycle_expr(v, &hs, path, via, len);
                    gops_heads_free(&hs);
                    free(path); free(via);
                    if (!c || !xl_push(out, c)) ok = 0;
                    break;
                }
                if (depth[y] < 0) {
                    sv[top] = y; se[top] = e; sd[top] = d + 1; top++;
                }
            }
            if (!ok) break;
        }
        if (!ok) break;
    }
done:
    free(depth); free(chain); free(inedge); free(sv); free(se); free(sd);
    return ok;
}

/* A shortest-found cycle through vertex r, by BFS from r: directed, the first
 * arc back into r; undirected, the first non-tree edge joining two different
 * branches of the BFS tree at r (or a second edge from a branch to r). Written
 * from r. O(V + E). Returns 0 on allocation failure; out stays empty if r lies
 * on no cycle. */
static int bfs_cycle_through(const GopsView* v, int directed, int r, XList* out) {
    const GopsInc* ip = gops_inc_cached(v, directed ? GOPS_INC_OUT : GOPS_INC_ALL);
    if (!ip) return 0;
    const GopsInc inc = *ip;             /* borrowed arrays, not freed here */
    size_t n = v->nv;
    int* par = gops_malloc(n, sizeof(int));        /* tree parent (-1 root)   */
    int* pe = gops_malloc(n, sizeof(int));         /* tree edge into vertex   */
    int* br = gops_malloc(n, sizeof(int));         /* branch (child of r)     */
    int* q = gops_malloc(n, sizeof(int));
    int* path = gops_malloc(n + 1, sizeof(int));
    int* via = gops_malloc(n + 2, sizeof(int));
    int ok = 1;
    if (!par || !pe || !br || !q || !path || !via) { ok = 0; goto done; }
    for (size_t i = 0; i < n; i++) par[i] = -2;
    size_t h = 0, t = 0;
    q[t++] = r; par[r] = -1; pe[r] = -1; br[r] = -1;
    int fx = -1, fy = -1, fe = -1;                 /* closing edge x - y      */
    while (h < t && fx < 0) {
        int x = q[h++];
        for (int j = inc.start[x]; j < inc.start[x + 1]; j++) {
            int y = inc.nbr[j], e = inc.eid[j];
            if (e == pe[x]) continue;
            if (y == r && x != r) { fx = x; fy = r; fe = e; break; }
            if (par[y] == -2) {
                par[y] = x; pe[y] = e; br[y] = (x == r) ? y : br[x];
                q[t++] = y;
            } else if (!directed && y != r && br[y] != br[x] && x != r) {
                fx = x; fy = y; fe = e; break;
            }
        }
    }
    if (fx >= 0) {
        /* r ... x along the tree, then x - fy, then fy ... r back up the tree */
        long len = 0;
        int* up = q;                               /* reuse: x's ancestors */
        long nu = 0;
        for (int z = fx; z != -1; z = par[z]) up[nu++] = z;
        for (long i = nu - 1; i >= 0; i--) {
            path[len] = up[i];
            via[len] = (len > 0) ? pe[up[i]] : -1;
            len++;
        }
        via[len] = fe;
        if (fy != r) {
            for (int z = fy; z != r; z = par[z]) {
                path[len] = z;
                via[len + 1] = pe[z];
                len++;
            }
        }
        GopsHeads hs = { { NULL, NULL } };
        Expr* c = cycle_expr(v, &hs, path, via, len);
        gops_heads_free(&hs);
        if (!c || !xl_push(out, c)) ok = 0;
    }
done:
    free(par); free(pe); free(br); free(q); free(path); free(via);
    return ok;
}

/* ---- Backtracking enumeration --------------------------------------------- */
typedef struct {
    const GopsView* v;
    GopsInc inc;
    int directed;
    long lo, hi;          /* length bounds (hi < 0: unbounded)               */
    long want;            /* results wanted (-1: all)                        */
    long steps;
    int* path; int* via; int* it; unsigned char* onpath;
    XList out;
    int oom, gave_up;
} Search;

static int search_init(Search* s, const GopsView* v, int directed, long lo, long hi, long want) {
    memset(s, 0, sizeof(*s));
    s->v = v; s->directed = directed; s->lo = lo; s->hi = hi; s->want = want;
    if (!gops_inc_build(v, directed ? GOPS_INC_OUT : GOPS_INC_ALL, &s->inc)) return 0;
    size_t n = v->nv;
    s->path = gops_malloc(n + 1, sizeof(int));
    s->via = gops_malloc(n + 2, sizeof(int));
    s->it = gops_malloc(n + 1, sizeof(int));
    s->onpath = gops_calloc(n, 1);
    if (!s->path || !s->via || !s->it || !s->onpath) return 0;
    return 1;
}

static void search_free(Search* s) {
    gops_inc_free(&s->inc);
    free(s->path); free(s->via); free(s->it); free(s->onpath);
    xl_free(&s->out);
}

static int search_tick(Search* s) {
    if ((++s->steps & 0xFFF) == 0) tc_check_deadline();
    if (s->steps > GOPS_SEARCH_MAX_STEPS) { s->gave_up = 1; return 0; }
    return 1;
}

static int search_done(const Search* s) {
    return s->oom || s->gave_up || (s->want >= 0 && (long)s->out.n >= s->want);
}

/* Simple cycles through `root`: with lowest_only, just those whose
 * lowest-positioned vertex is root (so enumerating every root finds each cycle
 * once); otherwise every cycle through root. */
static void cycles_from(Search* s, int root, int lowest_only) {
    const GopsInc* inc = &s->inc;
    int minlen = s->directed ? 2 : 3;
    long d = 0;                           /* path[0..d], edges via[1..d] */
    s->path[0] = root; s->onpath[root] = 1; s->it[0] = inc->start[root]; s->via[0] = -1;
    while (d >= 0 && !search_done(s)) {
        int x = s->path[d];
        if (s->it[d] == inc->start[x + 1]) {           /* exhausted: backtrack */
            s->onpath[x] = 0; d--;
            continue;
        }
        if (!search_tick(s)) break;
        int j = s->it[d]++;
        int y = inc->nbr[j], e = inc->eid[j];
        if (!s->directed && e == s->via[d]) continue;
        if (y == root) {
            long len = d + 1;
            if (len < minlen || len < s->lo) continue;
            if (s->hi >= 0 && len > s->hi) continue;
            /* undirected: each cycle once, in the direction whose second vertex
             * precedes its last */
            if (!s->directed && s->path[1] > s->path[d]) continue;
            s->via[d + 1] = e;
            GopsHeads hs = { { NULL, NULL } };
            Expr* c = cycle_expr(s->v, &hs, s->path, s->via, len);
            gops_heads_free(&hs);
            if (!c || !xl_push(&s->out, c)) { s->oom = 1; break; }
            continue;
        }
        if ((lowest_only && y < root) || s->onpath[y]) continue;
        if (s->hi >= 0 && d + 1 >= s->hi) continue;     /* no room to close */
        d++;
        s->path[d] = y; s->via[d] = e; s->onpath[y] = 1; s->it[d] = inc->start[y];
    }
    for (long i = 0; i <= d; i++) s->onpath[s->path[i]] = 0;
}

Expr* builtin_find_cycle(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;
    long lo = 0, hi = -1, want = 1;
    if (argc >= 2 && !parse_kspec(res->data.function.args[1], &lo, &hi)) return NULL;
    if (argc == 3 && !parse_count(res->data.function.args[2], &want)) return NULL;
    /* FindCycle[{g, v}, ...]: cycles through the vertex v. */
    const Expr* g = res->data.function.args[0];
    int through = -1;
    if (graph_is_list(g) && g->data.function.arg_count == 2
        && graph_is_valid(g->data.function.args[0])) {
        through = graph_vertex_position(g->data.function.args[0], g->data.function.args[1]);
        if (through < 0) return NULL;
        g = g->data.function.args[0];
    }
    GopsView v;
    if (!gops_view(g, &v)) return NULL;
    if (v.ndir != 0 && v.ndir != v.ne) return NULL;              /* mixed */
    if (argc >= 2 && v.weights) return NULL;
    int directed = (v.ndir == v.ne && v.ne > 0);
    int minlen = directed ? 2 : 3;

    if (want == 1 && hi < 0 && lo <= minlen) {                   /* any cycle */
        XList out = { NULL, 0, 0 };
        int ok = through >= 0 ? bfs_cycle_through(&v, directed, through, &out)
                              : dfs_any_cycle(&v, directed, &out);
        if (!ok) { xl_free(&out); return NULL; }
        return xl_list(&out);
    }
    Search s;
    if (!search_init(&s, &v, directed, lo, hi, want)) { search_free(&s); return NULL; }
    if (through >= 0) cycles_from(&s, through, 0);
    else for (size_t r = 0; r < v.nv && !search_done(&s); r++) cycles_from(&s, (int)r, 1);
    Expr* out = (s.oom || s.gave_up) ? NULL : xl_list(&s.out);
    search_free(&s);
    return out;
}

/* ---- FindPath ------------------------------------------------------------- */

static Expr* path_expr(const GopsView* v, const int* path, long len) {
    Expr** a = gops_malloc((size_t)len, sizeof(Expr*));
    if (!a) return NULL;
    for (long i = 0; i < len; i++) a[i] = expr_copy(v->verts[path[i]]);
    return gops_list_take(a, (size_t)len);
}

/* First DFS path s -> t (visited marking; neighbours in EdgeList order). The
 * forward-arc CSR comes from the per-graph cache, so a repeated query on one
 * graph is a bare DFS; the stack keeps (vertex, cursor) pairs contiguous. */
static Expr* dfs_path(const GopsView* v, int s, int t) {
    const GopsInc* inc = gops_inc_cached(v, GOPS_INC_OUT_NBR);
    if (!inc) return NULL;
    size_t n = v->nv;
    unsigned char* seen = gops_calloc(n, 1);
    int* sx = gops_malloc(n, sizeof(int));
    int* sp = gops_malloc(n, sizeof(int));
    Expr* out = NULL;
    if (!seen || !sx || !sp) goto done;
    const int* start = inc->start;
    const int* nbr = inc->nbr;
    size_t top = 0;
    sx[0] = s; sp[0] = start[s]; seen[s] = 1; top = 1;
    int found = (s == t);
    while (top > 0 && !found) {
        int x = sx[top - 1], p = sp[top - 1], end = start[x + 1];
        while (p < end && seen[nbr[p]]) p++;
        if (p == end) { top--; continue; }
        int y = nbr[p];
        sp[top - 1] = p + 1;
        seen[y] = 1;
        sx[top] = y; sp[top] = start[y]; top++;
        if (y == t) found = 1;
    }
    if (found) out = wrap1(path_expr(v, sx, (long)top));
    else out = empty_list();
done:
    free(seen); free(sx); free(sp);
    return out;
}

static int cmp_len(const void* a, const void* b) {
    /* stable via the index stored alongside */
    const long* x = a; const long* y = b;
    if (x[0] != y[0]) return (x[0] > y[0]) - (x[0] < y[0]);
    return (x[1] > y[1]) - (x[1] < y[1]);
}

Expr* builtin_find_path(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3 || argc > 5) return NULL;
    const Expr* g = res->data.function.args[0];
    long lo = 0, hi = -1, want = 1;
    if (argc >= 4 && !parse_kspec(res->data.function.args[3], &lo, &hi)) return NULL;
    if (argc == 5 && !parse_count(res->data.function.args[4], &want)) return NULL;
    int s = graph_vertex_position(g, res->data.function.args[1]);
    int t = graph_vertex_position(g, res->data.function.args[2]);
    if (s < 0 || t < 0) return NULL;
    GopsView v;
    if (!gops_view(g, &v)) return NULL;
    if (argc >= 4 && v.weights) return NULL;
    if (s == t) return empty_list();
    if (want == 1 && hi < 0 && lo <= 1) return dfs_path(&v, s, t);

    Search S;
    if (!search_init(&S, &v, 1, lo, hi, want)) { search_free(&S); return NULL; }
    const GopsInc* inc = &S.inc;
    long d = 0;
    S.path[0] = s; S.onpath[s] = 1; S.it[0] = inc->start[s];
    while (d >= 0 && !search_done(&S)) {
        int x = S.path[d];
        if (S.it[d] == inc->start[x + 1]) { S.onpath[x] = 0; d--; continue; }
        if (!search_tick(&S)) break;
        int y = inc->nbr[S.it[d]++];
        if (S.onpath[y]) continue;
        if (S.hi >= 0 && d + 1 > S.hi) continue;
        if (y == t) {
            if (d + 1 >= S.lo) {
                S.path[d + 1] = t;
                Expr* p = path_expr(&v, S.path, d + 2);
                if (!p || !xl_push(&S.out, p)) { S.oom = 1; break; }
            }
            continue;
        }
        if (S.hi >= 0 && d + 2 > S.hi) continue;        /* no room to reach t */
        d++;
        S.path[d] = y; S.onpath[y] = 1; S.it[d] = inc->start[y];
    }
    Expr* out = NULL;
    if (!S.oom && !S.gave_up) {
        /* discovery order reversed, then stably shortest first */
        size_t m = S.out.n;
        long* key = gops_malloc(2 * m, sizeof(long));
        Expr** sorted = gops_malloc(m, sizeof(Expr*));
        if (key && sorted) {
            for (size_t i = 0; i < m; i++) {
                key[2 * i] = (long)S.out.a[m - 1 - i]->data.function.arg_count;
                key[2 * i + 1] = (long)i;
            }
            qsort(key, m, 2 * sizeof(long), cmp_len);
            for (size_t i = 0; i < m; i++) sorted[i] = S.out.a[m - 1 - key[2 * i + 1]];
            free(S.out.a);
            S.out.a = NULL; S.out.n = 0; S.out.cap = 0;
            out = gops_list_take(sorted, m);
            sorted = NULL;
        }
        free(key); free(sorted);
    }
    search_free(&S);
    return out;
}
