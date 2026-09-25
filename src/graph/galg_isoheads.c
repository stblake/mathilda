/* galg_isoheads.c - IsomorphicGraphQ, FindGraphIsomorphism, CanonicalGraph,
 * GraphAutomorphismGroup: the Expr-facing side of the isomorphism engine in
 * galg_iso.c.
 *
 * A Graph is reduced to a vertex-coloured simple structure with up to three
 * relations (undirected neighbours, directed out- and in-neighbours):
 *   - self-loops become part of the vertex colour (undirected and directed
 *     loop counts, paired into one int);
 *   - an edge class of multiplicity 1 is an ordinary adjacency;
 *   - an edge class (u, v) of multiplicity k > 1 is subdivided: a new vertex w
 *     joined to u and v (u -> w -> v when directed), coloured by (k, kind) in
 *     a colour range disjoint from the original vertices.
 * Isomorphisms of the reduced structures restrict bijectively to isomorphisms
 * of the multigraphs (each subdivision vertex is determined by the images of
 * its two ends), so mixed graphs, loops and multigraphs are all supported --
 * Mathematica leaves those unevaluated ("not implemented"). Edge weights and
 * other properties are ignored, as in Mathematica.
 *
 * Every search is budgeted (GI_BUDGET nodes) and polls the TimeConstrained
 * deadline inside the engine; on budget exhaustion the head stays
 * unevaluated rather than guessing.
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define GI_BUDGET 50000000L
#define GI_SUBBASE (1 << 30)
#define GI_MAX_ENUM_CELLS (1L << 27)   /* ints held by FindGraphIsomorphism[.., All] */

typedef struct {
    GalgIsoGraph G;
    int n0;           /* original vertices (0..n0-1)          */
    long mu, md;      /* undirected / directed adjacencies     */
    int *uoff, *uadj, *ooff, *oadj, *ioff, *iadj, *vcol;
} GiBuilt;

static void gi_free(GiBuilt* B) {
    free(B->uoff); free(B->uadj); free(B->ooff); free(B->oadj);
    free(B->ioff); free(B->iadj); free(B->vcol);
    memset(B, 0, sizeof(*B));
}

/* Rearrange keys x*n + y (0 <= x, y < n) so that equal keys are contiguous,
 * in O(n + cnt): counting sort by x, then within each x-row group by y with
 * stamped per-y counters (no comparison sort; rows need no internal order).
 * Returns 0 on allocation failure. */
static int gi_group(int64_t* k, long cnt, int n) {
    if (cnt < 2) return 1;
    int* off = calloc((size_t)n + 1, sizeof(int));
    int* cy = malloc((size_t)n * sizeof(int));
    int* py = malloc((size_t)n * sizeof(int));
    int* st = malloc((size_t)n * sizeof(int));
    int64_t* tmp = malloc((size_t)cnt * sizeof(int64_t));
    int ok = off && cy && py && st && tmp;
    if (ok) {
        for (long i = 0; i < cnt; i++) off[(int)(k[i] / n) + 1]++;
        for (int x = 0; x < n; x++) off[x + 1] += off[x];
        for (long i = 0; i < cnt; i++) tmp[off[(int)(k[i] / n)]++] = k[i];
        for (int x = n; x > 0; x--) off[x] = off[x - 1];
        off[0] = 0;
        for (int y = 0; y < n; y++) st[y] = -1;
        for (int x = 0; x < n; x++) {
            int a = off[x], b = off[x + 1];
            if (b - a < 2) { for (int i = a; i < b; i++) k[i] = tmp[i]; continue; }
            for (int i = a; i < b; i++) {
                int y = (int)(tmp[i] % n);
                if (st[y] != x) { st[y] = x; cy[y] = 0; }
                cy[y]++;
            }
            int w = a;
            for (int i = a; i < b; i++) {
                int y = (int)(tmp[i] % n);
                if (st[y] == x) { st[y] = n + x; py[y] = w; w += cy[y]; }
                k[py[y]++] = tmp[i];
            }
        }
    }
    free(off); free(cy); free(py); free(st); free(tmp);
    return ok;
}

/* CSR from arc list (a[i] -> b[i]); sym adds the reverse arc too. */
static int gi_csr(int n, long m, const int* a, const int* b, int sym,
                  int** offp, int** adjp) {
    long h = sym ? 2 * m : m;
    int* off = calloc((size_t)n + 1, sizeof(int));
    int* adj = malloc((size_t)(h > 0 ? h : 1) * sizeof(int));
    int* pos = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!off || !adj || !pos) { free(off); free(adj); free(pos); return 0; }
    for (long i = 0; i < m; i++) { off[a[i] + 1]++; if (sym) off[b[i] + 1]++; }
    for (int v = 0; v < n; v++) off[v + 1] += off[v];
    for (int v = 0; v < n; v++) pos[v] = off[v];
    for (long i = 0; i < m; i++) {
        adj[pos[a[i]]++] = b[i];
        if (sym) adj[pos[b[i]]++] = a[i];
    }
    free(pos);
    *offp = off; *adjp = adj;
    return 1;
}

/* Build the reduced structure of a valid graph. Returns 1, or 0 on
 * allocation failure / absurd loop counts. */
static int gi_build(const Expr* g, GiBuilt* B) {
    memset(B, 0, sizeof(*B));
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return 0;
    int n = galg_nv(g);
    long m = galg_ne(g);
    int ok = 0;
    int* lu = calloc((size_t)n + 1, sizeof(int));
    int* ld = calloc((size_t)n + 1, sizeof(int));
    int64_t* ku = malloc((size_t)(m > 0 ? m : 1) * sizeof(int64_t));
    int64_t* kd = malloc((size_t)(m > 0 ? m : 1) * sizeof(int64_t));
    int *ua = NULL, *ub = NULL, *da = NULL, *db = NULL, *col = NULL;
    if (!lu || !ld || !ku || !kd) goto done;
    long nu = 0, nd = 0;
    int anyloop = 0;
    for (long k = 0; k < m; k++) {
        int a = eu[k], b = ev[k];
        if (a == b) {
            if (edir[k]) ld[a]++; else lu[a]++;
            anyloop = 1;
        } else if (edir[k]) {
            kd[nd++] = (int64_t)a * n + b;
        } else {
            int x = a < b ? a : b, y = a < b ? b : a;
            ku[nu++] = (int64_t)x * n + y;
        }
    }
    if (!gi_group(ku, nu, n) || !gi_group(kd, nd, n)) goto done;
    /* count subdivision vertices */
    long extra = 0;
    for (long i = 0; i < nu; ) { long j = i; while (j < nu && ku[j] == ku[i]) j++; if (j - i > 1) extra++; i = j; }
    for (long i = 0; i < nd; ) { long j = i; while (j < nd && kd[j] == kd[i]) j++; if (j - i > 1) extra++; i = j; }
    if ((long)n + extra > 0x7ffffff0L) goto done;
    int nt = n + (int)extra;
    size_t cap = (size_t)(nu + 2 * nd + 2 * extra + 1);
    ua = malloc(cap * sizeof(int)); ub = malloc(cap * sizeof(int));
    da = malloc(cap * sizeof(int)); db = malloc(cap * sizeof(int));
    if (!ua || !ub || !da || !db) goto done;
    int needcol = anyloop || extra > 0;
    if (needcol) {
        col = calloc((size_t)nt + 1, sizeof(int));
        if (!col) goto done;
        for (int v = 0; v < n; v++) {
            long s = (long)lu[v] + ld[v];
            long c = s * (s + 1) / 2 + ld[v];     /* Cantor pairing */
            if (c >= GI_SUBBASE) goto done;
            col[v] = (int)c;
        }
    }
    long mu = 0, md = 0;
    int w = n;
    for (long i = 0; i < nu; ) {
        long j = i; while (j < nu && ku[j] == ku[i]) j++;
        int x = (int)(ku[i] / n), y = (int)(ku[i] % n);
        long mult = j - i;
        if (mult == 1) { ua[mu] = x; ub[mu++] = y; }
        else {
            if (mult > (INT32_MAX - GI_SUBBASE) / 2 - 1) goto done;
            col[w] = GI_SUBBASE + 2 * (int)mult;
            ua[mu] = x; ub[mu++] = w; ua[mu] = y; ub[mu++] = w;
            w++;
        }
        i = j;
    }
    for (long i = 0; i < nd; ) {
        long j = i; while (j < nd && kd[j] == kd[i]) j++;
        int x = (int)(kd[i] / n), y = (int)(kd[i] % n);
        long mult = j - i;
        if (mult == 1) { da[md] = x; db[md++] = y; }
        else {
            if (mult > (INT32_MAX - GI_SUBBASE) / 2 - 1) goto done;
            col[w] = GI_SUBBASE + 2 * (int)mult + 1;
            da[md] = x; db[md++] = w; da[md] = w; db[md++] = y;
            w++;
        }
        i = j;
    }
    B->n0 = n; B->mu = mu; B->md = md;
    B->G.n = nt;
    if (mu > 0) {
        if (!gi_csr(nt, mu, ua, ub, 1, &B->uoff, &B->uadj)) goto done;
        B->G.uoff = B->uoff; B->G.uadj = B->uadj;
    }
    if (md > 0) {
        if (!gi_csr(nt, md, da, db, 0, &B->ooff, &B->oadj)) goto done;
        if (!gi_csr(nt, md, db, da, 0, &B->ioff, &B->iadj)) goto done;
        B->G.ooff = B->ooff; B->G.oadj = B->oadj;
        B->G.ioff = B->ioff; B->G.iadj = B->iadj;
    }
    B->vcol = col; col = NULL;
    B->G.vcol = B->vcol;
    ok = 1;
done:
    free(lu); free(ld); free(ku); free(kd);
    free(ua); free(ub); free(da); free(db); free(col);
    if (!ok) gi_free(B);
    return ok;
}

static int gi_cmpint(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

/* Cheap invariants: 1 if the pair may be isomorphic, 0 if certainly not,
 * -1 on allocation failure. */
static int gi_compatible(const GiBuilt* a, const GiBuilt* b) {
    if (a->n0 != b->n0 || a->G.n != b->G.n || a->mu != b->mu || a->md != b->md) return 0;
    if ((a->vcol == NULL) != (b->vcol == NULL)) return 0;
    int n = a->G.n;
    int* da = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* dbb = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!da || !dbb) { free(da); free(dbb); return -1; }
    int r = 1;
    /* undirected degree sequences: histogram, O(n) */
    if (a->uoff) {
        memset(da, 0, (size_t)n * sizeof(int));
        for (int v = 0; v < n; v++) da[a->uoff[v + 1] - a->uoff[v]]++;
        for (int v = 0; v < n && r; v++) if (--da[b->uoff[v + 1] - b->uoff[v]] < 0) r = 0;
    }
    /* colours and (out, in) degree pairs: sorted multisets, only when present */
    for (int pass = 0; pass < 2 && r; pass++) {
        if (pass == 0 && !a->vcol) continue;
        if (pass == 1 && !a->ooff) continue;
        const GiBuilt* s[2] = { a, b };
        int* d[2] = { da, dbb };
        for (int t = 0; t < 2; t++) {
            for (int v = 0; v < n; v++) {
                if (pass == 0) d[t][v] = s[t]->vcol[v];
                else d[t][v] = (s[t]->ooff[v + 1] - s[t]->ooff[v]) * 65536
                             + (s[t]->ioff[v + 1] - s[t]->ioff[v]);
            }
            qsort(d[t], (size_t)n, sizeof(int), gi_cmpint);
        }
        if (memcmp(da, dbb, (size_t)n * sizeof(int)) != 0) r = 0;
    }
    free(da); free(dbb);
    return r;
}

/* ---- IsomorphicGraphQ ------------------------------------------------------ */

Expr* builtin_isomorphic_graph_q(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    for (size_t i = 0; i < argc; i++)
        if (!graph_is_valid(res->data.function.args[i])) return galg_truth(0);
    GiBuilt A;
    if (!gi_build(res->data.function.args[0], &A)) return NULL;
    int* map = malloc((size_t)(A.G.n > 0 ? A.G.n : 1) * sizeof(int));
    if (!map) { gi_free(&A); return NULL; }
    int answer = 1;
    for (size_t i = 1; i < argc && answer == 1; i++) {
        GiBuilt H;
        if (!gi_build(res->data.function.args[i], &H)) { answer = -1; break; }
        int c = gi_compatible(&A, &H);
        if (c == 1 && A.G.n > 0) c = galg_iso_find(&A.G, &H.G, map, GI_BUDGET);
        gi_free(&H);
        answer = c;
    }
    free(map);
    gi_free(&A);
    return answer < 0 ? NULL : galg_truth(answer);
}

/* ---- FindGraphIsomorphism --------------------------------------------------- */

typedef struct {
    int n0;
    int* maps;       /* count * n0 */
    long count, cap;
    int overflow;
} GiCollect;

static int gi_collect_cb(const int* map, void* ctx) {
    GiCollect* C = (GiCollect*)ctx;
    if (C->count == C->cap) {
        long nc = C->cap ? 2 * C->cap : 16;
        if (nc * (long)(C->n0 > 0 ? C->n0 : 1) > GI_MAX_ENUM_CELLS) { C->overflow = 1; return 0; }
        int* nm = realloc(C->maps, (size_t)nc * (size_t)(C->n0 > 0 ? C->n0 : 1) * sizeof(int));
        if (!nm) { C->overflow = 1; return 0; }
        C->maps = nm; C->cap = nc;
    }
    memcpy(C->maps + (size_t)C->count * (size_t)C->n0, map, (size_t)C->n0 * sizeof(int));
    C->count++;
    return 1;
}

static Expr* gi_assoc(const Expr* g, const Expr* h, const int* map, int n) {
    const Expr* gv = g->data.function.args[0];
    const Expr* hv = h->data.function.args[0];
    Expr** rules = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    if (!rules) return NULL;
    for (int v = 0; v < n; v++) {
        Expr* kv[2];
        kv[0] = expr_copy(gv->data.function.args[v]);
        kv[1] = expr_copy(hv->data.function.args[map[v]]);
        rules[v] = expr_new_function(expr_new_symbol(SYM_Rule), kv, 2);
    }
    Expr* a = expr_new_function(expr_new_symbol(SYM_Association), rules, (size_t)n);
    free(rules);
    return a;
}

Expr* builtin_find_graph_isomorphism(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* h = res->data.function.args[1];
    if (!graph_is_valid(g) || !graph_is_valid(h)) return NULL;
    long count = 1;
    if (argc == 3 && !galg_parse_count(res->data.function.args[2], &count)) return NULL;
    GiBuilt A, H;
    if (!gi_build(g, &A)) return NULL;
    if (!gi_build(h, &H)) { gi_free(&A); return NULL; }
    Expr* out = NULL;
    int c = gi_compatible(&A, &H);
    if (c == 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else if (c == 1) {
        GiCollect C;
        memset(&C, 0, sizeof(C));
        C.n0 = A.n0;
        long r;
        static const int empty_map[1] = { 0 };
        if (A.G.n == 0) r = gi_collect_cb(empty_map, &C);   /* the empty map */
        else if (count == 1) {
            /* one map: a plain search (enumerate would also prove there is
             * no second one to skip) */
            int* map = malloc((size_t)A.G.n * sizeof(int));
            int f = map ? galg_iso_find(&A.G, &H.G, map, GI_BUDGET) : -1;
            r = f < 0 ? -1 : (f == 1 ? gi_collect_cb(map, &C) : 0);
            free(map);
        }
        else r = galg_iso_enumerate(&A.G, &H.G, count, gi_collect_cb, &C, GI_BUDGET);
        if (r >= 0 && !C.overflow) {
            Expr** items = malloc((size_t)(C.count > 0 ? C.count : 1) * sizeof(Expr*));
            if (items) {
                long k = 0;
                for (; k < C.count; k++) {
                    items[k] = gi_assoc(g, h, C.maps + (size_t)k * (size_t)C.n0, C.n0);
                    if (!items[k]) break;
                }
                if (k == C.count)
                    out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
                else
                    for (long j = 0; j < k; j++) expr_free(items[j]);
                free(items);
            }
        }
        free(C.maps);
    }
    gi_free(&A); gi_free(&H);
    return out;
}

/* ---- CanonicalGraph --------------------------------------------------------- */

typedef struct { int a, b, d; } GiEdge;

static int gi_cmpedge(const void* x, const void* y) {
    const GiEdge* p = (const GiEdge*)x; const GiEdge* q = (const GiEdge*)y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    return (p->d > q->d) - (p->d < q->d);
}

/* rank[v] (0-based canonical position of original vertex v), or NULL. */
static int* gi_canonical_rank(const Expr* g, int* nout) {
    GiBuilt A;
    if (!gi_build(g, &A)) return NULL;
    int nt = A.G.n, n = A.n0;
    int* lab = malloc((size_t)(nt > 0 ? nt : 1) * sizeof(int));
    int* rank = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int ok = lab && rank && (nt == 0 || galg_iso_canon(&A.G, lab, GI_BUDGET) == 1);
    if (ok) {
        int r = 0;
        for (int i = 0; i < nt; i++) if (lab[i] < n) rank[lab[i]] = r++;
        ok = (r == n);
    }
    free(lab);
    gi_free(&A);
    if (!ok) { free(rank); return NULL; }
    *nout = n;
    return rank;
}

Expr* builtin_canonical_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = 0;
    int* rank = gi_canonical_rank(g, &n);
    if (!rank) return NULL;
    const int *eu, *ev;
    const unsigned char* edir;
    long m = galg_ne(g);
    GiEdge* E = malloc((size_t)(m > 0 ? m : 1) * sizeof(GiEdge));
    Expr** vitems = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    Expr** eitems = malloc((size_t)(m > 0 ? m : 1) * sizeof(Expr*));
    Expr* out = NULL;
    if (E && vitems && eitems && graph_edge_indices(g, &eu, &ev, &edir)) {
        for (long k = 0; k < m; k++) {
            int a = rank[eu[k]] + 1, b = rank[ev[k]] + 1;
            if (!edir[k] && a > b) { int t = a; a = b; b = t; }
            E[k].a = a; E[k].b = b; E[k].d = edir[k] ? 1 : 0;
        }
        qsort(E, (size_t)m, sizeof(GiEdge), gi_cmpedge);
        for (int i = 0; i < n; i++) vitems[i] = expr_new_integer(i + 1);
        Expr* hu = expr_new_symbol(SYM_UndirectedEdge);
        Expr* hd = expr_new_symbol(SYM_DirectedEdge);
        for (long k = 0; k < m; k++) {
            Expr* ab[2] = { expr_new_integer(E[k].a), expr_new_integer(E[k].b) };
            eitems[k] = expr_new_function(expr_copy(E[k].d ? hd : hu), ab, 2);
        }
        expr_free(hu); expr_free(hd);
        Expr* ge[2];
        ge[0] = expr_new_function(expr_new_symbol(SYM_List), vitems, (size_t)n);
        ge[1] = expr_new_function(expr_new_symbol(SYM_List), eitems, (size_t)m);
        out = expr_new_function(expr_new_symbol(SYM_Graph), ge, 2);
    }
    free(E); free(vitems); free(eitems); free(rank);
    return out;
}

/* ---- GraphAutomorphismGroup ------------------------------------------------- */

/* Cycles[{{..}, ..}] of a permutation of 0..n-1 (1-based points, each cycle
 * starting at its least point, cycles by least point, fixed points dropped);
 * NULL for the identity. */
static Expr* gi_cycles(const int* p, int n, char* seen) {
    memset(seen, 0, (size_t)n);
    long ncyc = 0;
    for (int v = 0; v < n; v++) if (!seen[v] && p[v] != v) {
        for (int w = v; !seen[w]; w = p[w]) seen[w] = 1;
        ncyc++;
    }
    if (ncyc == 0) return NULL;
    memset(seen, 0, (size_t)n);
    Expr** cyc = malloc((size_t)ncyc * sizeof(Expr*));
    int* buf = malloc((size_t)n * sizeof(int));
    if (!cyc || !buf) { free(cyc); free(buf); return NULL; }
    long c = 0;
    for (int v = 0; v < n; v++) if (!seen[v] && p[v] != v) {
        int len = 0;
        for (int w = v; !seen[w]; w = p[w]) { seen[w] = 1; buf[len++] = w; }
        Expr** pts = malloc((size_t)len * sizeof(Expr*));
        if (!pts) { for (long j = 0; j < c; j++) expr_free(cyc[j]); free(cyc); free(buf); return NULL; }
        for (int i = 0; i < len; i++) pts[i] = expr_new_integer(buf[i] + 1);
        cyc[c++] = expr_new_function(expr_new_symbol(SYM_List), pts, (size_t)len);
        free(pts);
    }
    free(buf);
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), cyc, (size_t)c);
    free(cyc);
    return expr_new_function(expr_new_symbol("Cycles"), &lst, 1);
}

Expr* builtin_graph_automorphism_group(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    GiBuilt A;
    if (!gi_build(g, &A)) return NULL;
    int* gens = NULL;
    long ng = A.G.n > 0 ? galg_iso_automorphisms(&A.G, &gens, GI_BUDGET) : 0;
    Expr* out = NULL;
    if (ng >= 0) {
        int n = A.n0, nt = A.G.n;
        Expr** items = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(Expr*));
        char* seen = malloc((size_t)(n > 0 ? n : 1));
        if (items && seen) {
            long k = 0;
            for (long i = 0; i < ng; i++) {
                Expr* cy = gi_cycles(gens + (size_t)i * (size_t)nt, n, seen);
                if (cy) items[k++] = cy;
            }
            Expr* lst = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
            out = expr_new_function(expr_new_symbol("PermutationGroup"), &lst, 1);
        }
        free(items); free(seen);
    }
    free(gens);
    gi_free(&A);
    return out;
}
