/* galg_planar.c - linear-time planarity test (left-right criterion).
 *
 *   int galg_planar_test(const GalgUG* u);   1 planar, 0 not, -1 alloc failure
 *
 * Algorithm
 * ---------
 * The left-right planarity test of de Fraysseix and Rosenstiehl in the
 * formulation of
 *
 *   U. Brandes, "The Left-Right Planarity Test", manuscript, 2009
 *   (algorithms 2-5; the same algorithm networkx's check_planarity uses).
 *
 * A graph is planar iff the back edges of a DFS tree can be split into a left
 * and a right class such that, at every fork of the tree, return edges that
 * would cross go to different sides.  The test runs in two DFS passes over
 * each connected component:
 *
 *  1. Orientation.  A DFS orients every edge away from the root (tree edges
 *     down, back edges up) and computes, per oriented edge e = (v,w):
 *       lowpt(e)   the lowest height reachable from w via tree edges below e
 *                  followed by one back edge (height(v) if none),
 *       lowpt2(e)  the second lowest such height (height(v) if none),
 *       nesting(e) = 2*lowpt(e) + [lowpt2(e) < height(v)]   ("chordal" +1).
 *     Ordering each vertex's outgoing edges by nesting depth puts the edge
 *     with the lowest return point first, and among equal lowpoints the
 *     non-chordal ones (a single return height) before the chordal ones,
 *     which is the order in which constraints can be merged greedily.
 *
 *  2. Testing.  A second DFS in nesting order maintains a stack S of
 *     conflict pairs (L,R).  L and R are intervals of return (back) edges,
 *     represented by their lowest and highest element; the edges inside an
 *     interval are chained from high to low through ref[].  All edges of one
 *     interval must be on the same side, and L and R on opposite sides.
 *     After returning from the subtree of an outgoing edge e_i of v that
 *     has return edges (lowpt(e_i) < height(v)):
 *       - for the first edge e_1 there is nothing to merge; its lowest
 *         return edge becomes lowpt_edge of v's parent edge;
 *       - for later e_i, add_constraints() merges all conflict pairs pushed
 *         while exploring e_i into one pair P: their R intervals are fused
 *         into P.R (a pair with both sides non-empty here is a
 *         contradiction), except intervals whose lowest return point is
 *         lowpt(e) itself: those are "aligned" with lowpt_edge(e), the
 *         lowest return edge of e_1, on its side.  Then every earlier pair
 *         (from e_1..e_{i-1}) whose interval reaches above lowpt(e_i)
 *         conflicts with e_i and is folded into P.L (if both of its sides
 *         conflict, the graph is not planar).
 *     When the DFS retreats over the tree edge e = (u,v), back edges ending
 *     at u are no longer relevant: remove_back_edges() pops the pairs whose
 *     lowest return point is u and trims the u-ending tops of the next pair
 *     by following ref[] down the interval.
 *     If no contradiction arises the graph is planar.  The embedding phase
 *     of Brandes' paper (sides, signs, final rotation system) is not needed
 *     for the test and is omitted.
 *
 * Complexity
 * ----------
 * O(n + m) time and memory.  Nesting depths lie in 0..2n+1, so the adjacency
 * lists are ordered by one global counting sort.  Pairs are pushed once per
 * back edge and at most once per add_constraints() call, so there are O(m)
 * pushes and hence O(m) pops; every trim step discards an edge from an
 * interval for good, so the testing phase is linear too.  Before any of this the
 * Euler bound rejects m > 3n - 6, so m = O(n) for the main phase.
 *
 * Iterative design
 * ----------------
 * Both DFS passes use an explicit vertex stack plus a per-vertex cursor into
 * its adjacency list, so path-like graphs with 10^6 vertices (DFS depth n)
 * need no call stack.  The work Brandes performs "after the recursive call"
 * for tree edge (v,w) is done when w is popped.  The conflict-pair stack is
 * a flat array of four edge ids per pair (at most m pairs are ever live);
 * stack_bottom[e] records the stack height when e was entered (comparing
 * heights is equivalent to Brandes' comparison of stack elements, since S
 * never drops below that height before e is finished).
 * Phase 1 numbers vertices in DFS preorder and phase 2 runs entirely on
 * those numbers, so DFS subtrees are contiguous index ranges; this roughly
 * halves the cache misses of phase 2 on graphs with scattered vertex ids.
 *
 * Memory: a single malloc'd block of ints, O(n + m), freed on every return
 * path.  A NULL graph also yields -1.
 * No Expr/CAS dependency; graph_algos.h is included for GalgUG only.
 */

#include "graph.h"
#include "graph_algos.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define NONE (-1)

/* One conflict pair: intervals L = [ll..lh], R = [rl..rh] (edge ids, NONE
 * when absent).  An interval is empty iff both ends are NONE. */
typedef struct { int ll, lh, rl, rh; } Pair;

typedef struct {
    const int* height;     /* per vertex: DFS depth                          */
    const int* lowpt;      /* per edge                                       */
    const int* esrc;       /* per edge: tail of the oriented edge            */
    const int* edst;       /* per edge: head                                 */
    int*       ref;        /* per edge: next lower edge in its interval      */
    int*       lowpt_edge; /* per edge: a return edge realizing lowpt        */
    int*       stack_bottom;
    Pair*      S;          /* conflict-pair stack                            */
    int        sp;         /* its height                                     */
} LRState;

static int iv_empty(int lo, int hi) { return lo == NONE && hi == NONE; }

/* Interval [lo..hi] conflicts with edge b if it returns above lowpt(b).
 * (A non-empty interval always has both ends; the hi check is defensive.) */
static int iv_conflicting(const LRState* st, int lo, int hi, int b) {
    return !iv_empty(lo, hi) && hi != NONE && st->lowpt[hi] > st->lowpt[b];
}

/* Lowest return point of a (non-empty) pair. */
static int pair_lowest(const LRState* st, const Pair* p) {
    if (iv_empty(p->ll, p->lh)) return st->lowpt[p->rl];
    if (iv_empty(p->rl, p->rh)) return st->lowpt[p->ll];
    return st->lowpt[p->ll] < st->lowpt[p->rl] ? st->lowpt[p->ll]
                                               : st->lowpt[p->rl];
}

static void pair_swap(Pair* p) {
    int t;
    t = p->ll; p->ll = p->rl; p->rl = t;
    t = p->lh; p->lh = p->rh; p->rh = t;
}

/* Merge the constraints of e_i (a non-first outgoing edge of the head of e)
 * with those of its elder siblings.  Returns 0 on a contradiction. */
static int add_constraints(LRState* st, int ei, int e) {
    const int* lowpt = st->lowpt;
    int* ref = st->ref;
    Pair P = { NONE, NONE, NONE, NONE };

    /* Pairs above stack_bottom[ei] stem from the subtree of e_i.  All of
     * them must be one-sided (a two-sided one certifies non-planarity, see
     * Brandes); fuse their R intervals into
     * P.R, linking the old low end to the next interval's high end through
     * ref[].  An interval whose lowest return is lowpt(e) cannot conflict
     * with e_1 and is aligned with lowpt_edge(e) instead of joining P. */
    do {
        Pair Q = st->S[--st->sp];
        if (!iv_empty(Q.ll, Q.lh)) pair_swap(&Q);
        if (!iv_empty(Q.ll, Q.lh)) return 0;           /* two-sided  */
        if (lowpt[Q.rl] > lowpt[e]) {
            if (iv_empty(P.rl, P.rh)) {                 /* topmost interval */
                P.rh = Q.rh;
            } else {
                ref[P.rl] = Q.rh;                       /* chain below P.R */
            }
            P.rl = Q.rl;
        } else {                                        /* align */
            ref[Q.rl] = st->lowpt_edge[e];
        }
    } while (st->sp != st->stack_bottom[ei]);

    /* Pairs of e_1..e_{i-1} that return above lowpt(e_i) conflict with it
     * (their return edges would cross e_i's lowest one): the conflicting
     * side is chained under P.L, the other side (which must not conflict)
     * under P.R.  As in Brandes, the scan stops at the first
     * non-conflicting pair. */
    while (st->sp > 0) {
        Pair* T = &st->S[st->sp - 1];
        Pair Q;
        if (!iv_conflicting(st, T->ll, T->lh, ei) &&
            !iv_conflicting(st, T->rl, T->rh, ei))
            break;
        Q = *T;
        st->sp--;
        if (iv_conflicting(st, Q.rl, Q.rh, ei)) pair_swap(&Q);
        if (iv_conflicting(st, Q.rl, Q.rh, ei)) return 0; /* both sides */
        if (P.rl != NONE) ref[P.rl] = Q.rh;
        if (Q.rl != NONE) P.rl = Q.rl;
        if (iv_empty(P.ll, P.lh)) {
            P.lh = Q.lh;
        } else {
            ref[P.ll] = Q.lh;
        }
        P.ll = Q.ll;
    }

    if (!iv_empty(P.ll, P.lh) || !iv_empty(P.rl, P.rh)) st->S[st->sp++] = P;
    return 1;
}

/* Retreating over tree edge e = (u,v): discard back edges that end at u. */
static void remove_back_edges(LRState* st, int e) {
    int* ref = st->ref;
    int u = st->esrc[e];
    int hu = st->height[u];

    /* Whole pairs whose lowest return point is u. */
    while (st->sp > 0 && pair_lowest(st, &st->S[st->sp - 1]) == hu)
        st->sp--;

    /* The next pair may still have u-ending edges at the top of its
     * intervals (the highest return points); walk them off through ref[].
     * An interval trimmed to nothing hands its ref to the opposite side's
     * low end, as in Brandes (only relevant for the embedding sides). */
    if (st->sp > 0) {
        Pair* P = &st->S[st->sp - 1];
        while (P->lh != NONE && st->edst[P->lh] == u) P->lh = ref[P->lh];
        if (P->lh == NONE && P->ll != NONE) {           /* just emptied */
            ref[P->ll] = P->rl;
            P->ll = NONE;
        }
        while (P->rh != NONE && st->edst[P->rh] == u) P->rh = ref[P->rh];
        if (P->rh == NONE && P->rl != NONE) {           /* just emptied */
            ref[P->rl] = P->ll;
            P->rl = NONE;
        }
    }
    /* (Brandes also sets ref[e] to the highest return edge here; that only
     * matters for the embedding phase.) */
}

/* The lowpoint bookkeeping once edge k = (v, .) is finished (hv = height of
 * v): fold its lowpoints into the parent edge e of v, then replace lowpt2[k],
 * which is dead from here on, by the nesting depth of k.  The lowpt2 array
 * thus holds nesting depths once phase 1 is over. */
static void finish_edge(int k, int e, int hv, int* lowpt, int* lowpt2) {
    int nest = 2 * lowpt[k] + (lowpt2[k] < hv);
    if (e != NONE) {
        if (lowpt[k] < lowpt[e]) {
            lowpt2[e] = lowpt[e] < lowpt2[k] ? lowpt[e] : lowpt2[k];
            lowpt[e] = lowpt[k];
        } else if (lowpt[k] > lowpt[e]) {
            if (lowpt[k] < lowpt2[e]) lowpt2[e] = lowpt[k];
        } else {
            if (lowpt2[k] < lowpt2[e]) lowpt2[e] = lowpt2[k];
        }
    }
    lowpt2[k] = nest;
}

int galg_planar_test(const GalgUG* u) {
    int n, m, i, v, w, k, sp;
    const int* off = u ? u->off : NULL;
    const int* adj = u ? u->adj : NULL;
    size_t nn, mm, nb, total;
    int *block, *p;
    int *hp, *height, *parent_edge, *vstack, *cur, *ooff;
    int *esrc, *edst, *lowpt, *lowpt2, *nesting, *oadj, *ref, *lowpt_edge,
        *stack_bottom, *bucket;
    Pair* S;
    LRState st;
    int result = 1;

    if (!u) return -1;
    n = u->n;
    if (n < 5) return 1;                       /* every graph on <= 4 vertices */
    if (u->m < 9) return 1;                    /* K5 and K3,3 need >= 9 edges  */
    if ((long long)u->m > 3LL * n - 6) return 0;  /* Euler bound            */
    if (u->m > INT_MAX / 2) return -1;         /* cannot happen with int CSR   */
    m = (int)u->m;

    /* One block of ints:
     *   by input vertex id:  hp[2n] (height, preorder number)
     *   by preorder number:  height[n], parent_edge[n], ooff[n+1]
     *   cur[n] (adjacency cursors: by input id in phase 1, by preorder
     *   number in phase 2), vstack[n] (DFS stack), bucket[2n+2]
     *   by edge id:          esrc, edst, lowpt, lowpt2/nesting, oadj, ref,
     *                        lowpt_edge, stack_bottom
     *   m conflict pairs (4 ints each). */
    nn = (size_t)n;
    mm = (size_t)m;
    nb = 2 * nn + 2;
    total = 7 * nn + 1 + nb + 8 * mm + 4 * mm;
    if (total > (size_t)-1 / sizeof(int)) return -1;
    block = (int*)malloc(total * sizeof(int));
    if (!block) return -1;
    p = block;
    hp = p;           p += 2 * nn;
    cur = p;          p += nn;
    height = p;       p += nn;
    parent_edge = p;  p += nn;
    vstack = p;       p += nn;
    ooff = p;         p += nn + 1;
    bucket = p;       p += nb;
    esrc = p;         p += mm;
    edst = p;         p += mm;
    lowpt = p;        p += mm;
    lowpt2 = p;       p += mm;
    oadj = p;         p += mm;
    ref = p;          p += mm;
    lowpt_edge = p;   p += mm;
    stack_bottom = p; p += mm;
    S = (Pair*)(void*)p;   /* Pair is four ints: same alignment as int */

    /* ---- Phase 1: orientation, lowpoints, nesting depths ---------------
     * The DFS runs over input ids, but numbers vertices in preorder and
     * records every oriented edge by preorder numbers, so that phase 2 works
     * on a relabelled graph whose DFS subtrees are contiguous ranges (much
     * better locality than arbitrary input ids).  hp[2v] is the height of v
     * (NONE = unvisited), hp[2v+1] its preorder number; interleaved so the
     * neighbour scan touches one cache line per neighbour. */
    for (i = 0; i < n; i++) hp[2 * i] = NONE;
    k = 0;
    {
        int npre = 0;
        for (i = 0; i < n; i++) {
            if (hp[2 * i] != NONE) continue;
            hp[2 * i] = 0;
            hp[2 * i + 1] = npre;
            height[npre] = 0;
            parent_edge[npre] = NONE;
            npre++;
            cur[i] = off[i];
            sp = 0;
            vstack[sp++] = i;
            while (sp > 0) {
                int e, hv, pv, qv, j, end;
                v = vstack[sp - 1];
                hv = hp[2 * v];
                qv = hp[2 * v + 1];
                e = parent_edge[qv];
                pv = e == NONE ? NONE : esrc[e];       /* parent, preorder */
                end = off[v + 1];
                /* Scan on until the next unvisited neighbour.  A visited w
                 * above v (other than the parent) closes a new back edge
                 * v->w; w below v is either a child or a descendant whose
                 * back edge to v was oriented when it scanned v. */
                for (j = cur[v]; j < end; j++) {
                    int hw, qw;
                    w = adj[j];
                    hw = hp[2 * w];
                    if (hw == NONE) break;
                    qw = hp[2 * w + 1];
                    if (hw < hv && qw != pv) {
                        esrc[k] = qv; edst[k] = qw;
                        lowpt[k] = hw;
                        lowpt2[k] = hv;
                        finish_edge(k, e, hv, lowpt, lowpt2);
                        k++;
                    }
                }
                if (j < end) {                         /* tree edge v->w */
                    w = adj[j];
                    cur[v] = j + 1;
                    hp[2 * w] = hv + 1;
                    hp[2 * w + 1] = npre;
                    height[npre] = hv + 1;
                    parent_edge[npre] = k;
                    npre++;
                    esrc[k] = qv; edst[k] = hp[2 * w + 1];
                    lowpt[k] = lowpt2[k] = hv;
                    k++;
                    cur[w] = off[w];
                    vstack[sp++] = w;                  /* finished on pop */
                } else {                               /* retreat over e */
                    sp--;
                    if (e != NONE)
                        finish_edge(e, parent_edge[pv], height[pv], lowpt,
                                    lowpt2);
                }
            }
        }
    }
    /* k == m: in a simple graph every edge is oriented exactly once.  From
     * here on vertices are preorder numbers; cur is reused per preorder. */

    /* ---- Order outgoing edges by nesting depth (counting sort) --------- */
    nesting = lowpt2;
    memset(bucket, 0, nb * sizeof(int));
    for (k = 0; k < m; k++) bucket[nesting[k]]++;
    {
        int s = 0, t;
        for (i = 0; i < (int)nb; i++) { t = bucket[i]; bucket[i] = s; s += t; }
    }
    for (k = 0; k < m; k++) ref[bucket[nesting[k]]++] = k;  /* ref as temp */
    memset(ooff, 0, (nn + 1) * sizeof(int));
    for (k = 0; k < m; k++) ooff[esrc[k] + 1]++;
    for (i = 0; i < n; i++) ooff[i + 1] += ooff[i];
    for (i = 0; i < n; i++) cur[i] = ooff[i];
    for (i = 0; i < m; i++) { k = ref[i]; oadj[cur[esrc[k]]++] = k; }
    for (k = 0; k < m; k++) { ref[k] = NONE; lowpt_edge[k] = NONE; }

    /* ---- Phase 2: testing ---------------------------------------------- */
    st.height = height;
    st.lowpt = lowpt;
    st.esrc = esrc;
    st.edst = edst;
    st.ref = ref;
    st.lowpt_edge = lowpt_edge;
    st.stack_bottom = stack_bottom;
    st.S = S;
    st.sp = 0;

    for (i = 0; i < n && result; i++) {
        if (parent_edge[i] != NONE) continue;          /* not a DFS root */
        cur[i] = ooff[i];
        sp = 0;
        vstack[sp++] = i;
        while (sp > 0) {
            int ei, x;
            v = vstack[sp - 1];
            if (cur[v] < ooff[v + 1]) {
                ei = oadj[cur[v]];
                w = edst[ei];
                st.stack_bottom[ei] = st.sp;
                if (parent_edge[w] == ei) {            /* tree edge: descend */
                    cur[w] = ooff[w];
                    vstack[sp++] = w;
                    continue;
                }
                /* back edge: a new one-sided pair R = [ei..ei] */
                lowpt_edge[ei] = ei;
                S[st.sp].ll = NONE; S[st.sp].lh = NONE;
                S[st.sp].rl = ei;   S[st.sp].rh = ei;
                st.sp++;
                x = v;
            } else {                                   /* retreat over e */
                int e = parent_edge[v];
                sp--;
                if (e == NONE) continue;
                remove_back_edges(&st, e);
                ei = e;
                x = esrc[e];
            }
            /* ei = (x, .) is finished: integrate its return edges. */
            if (lowpt[ei] < height[x]) {
                int e = parent_edge[x];
                if (ei == oadj[ooff[x]]) {
                    lowpt_edge[e] = lowpt_edge[ei];
                } else if (!add_constraints(&st, ei, e)) {
                    result = 0;
                    break;
                }
            }
            cur[x]++;
        }
    }

    free(block);
    return result;
}

/* ---- PlanarGraphQ ----------------------------------------------------------
 * PlanarGraphQ[g]: the test on the underlying simple undirected graph (edge
 * direction and anti-parallel pairs do not affect planarity). Like every *Q
 * predicate it gives False for a non-graph. Memory: returns a fresh symbol. */
Expr* builtin_planar_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return galg_truth(0);
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int r = galg_planar_test(u);
    galg_ug_free(u);
    return r < 0 ? NULL : galg_truth(r);
}
