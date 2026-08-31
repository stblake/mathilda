/* vertexcoloring.c - FindVertexColoring[g]: a MINIMAL vertex colouring.
 *
 * Wolfram's FindVertexColoring returns a colouring whose number of distinct
 * colours equals the chromatic number. That word "minimal" is the whole
 * difficulty: computing it is NP-hard. A greedy or DSATUR-only implementation
 * returns a valid-but-frequently-larger colouring, and it fails SILENTLY -- a
 * plausible list of integers that quietly contradicts the documented meaning.
 * So the search here is exact:
 *
 *   ub = fvc_dsatur_bound()  -- a good upper bound, cheap, and a real colouring
 *   lb = fvc_clique_bound()  -- a greedy-clique lower bound
 *   if lb == ub              -- ub is proven optimal; answer with no search
 *   else                     -- DSATUR branch-and-bound (fvc_bb), improving the
 *                               incumbent until it meets lb or the tree is
 *                               exhausted
 *
 * The lower bound is not an optimisation. Without it CompleteGraph[128] -- under
 * the vertex cap, so accepted -- would search instead of answering from the
 * bounds, which for a complete graph is hopeless. With it lb == ub and the
 * answer is immediate, at zero search nodes.
 *
 * Neither bound makes the search cheap in general, though: see FVC_MAX_STEPS.
 * Exactness is guaranteed by REFUSING (unevaluated) whenever it cannot be
 * proven, never by returning a merely-valid colouring.
 *
 * This is Wolfram's own "BacktrackingDS" method, so shipping only it is a
 * documented subset rather than a divergence. No Method option is offered.
 *
 * Adjacency is the UNDIRECTED neighbourhood: an edge constrains its endpoints
 * whichever way it points. GraphAdj stores successors in out[] and predecessors
 * in in[], with an UndirectedEdge contributing to both, so the neighbourhood of
 * v is out[v] together with in[v] -- the same walk graph_count_components does.
 * It is walked IN PLACE; no union structure is materialised, so there is
 * nothing beyond the GraphAdj itself for a caller to free.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; the evaluator
 * frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "core.h"     /* tc_check_deadline: cooperative TimeConstrained abort */
#include <stdlib.h>

/* An exact chromatic-number search is NP-hard, so refuse outright above a fixed
 * size rather than build an unbounded amount of scratch state. 128 is chosen on
 * TYPICAL rather than worst-case cost: sparse graphs at this order solve in well
 * under a second. Same shape as FM_MAX_CON (src/solve/reduce_fm.c:18).
 *
 * This is a guard on SIZE only. It is not what keeps the head responsive --
 * TimeConstrained is (see FVC_MAX_STEPS below), and FVC_MAX_STEPS is only the
 * backstop for when nobody is there to interrupt. */
#define FVC_MAX_VERTICES 128

/* The vertex cap bounds SIZE, but size does not bound COST: exact colouring is
 * exponential, and a dense random graph well under the cap runs unboundedly.
 * Measured on this implementation at density ~0.24: n=80 solves in 0.33s,
 * n=100 in 30s, n=128 not within a minute. So the cap alone would still hang.
 *
 * Hence a second, independent guard: a budget on branch-and-bound NODES, after
 * which the search gives up and the head returns unevaluated exactly as it does
 * above the vertex cap. A node COUNT rather than a wall clock deliberately --
 * a time-based cutoff would make the ANSWER machine-dependent (a fast host
 * proves minimality, a slow one refuses the same graph), whereas a node count
 * is deterministic: the same graph gets the same answer on every machine.
 *
 * THE BUDGET IS A BACKSTOP, NOT THE RESPONSIVENESS MECHANISM. TimeConstrained
 * is: it already interrupts long-running pure-C builtins (SIGPROF/siglongjmp,
 * src/core.c:4038), fvc_bb polls tc_check_deadline() so it is honoured even
 * where the signal is unreliable, and it is the idiomatic, user-facing way to
 * bound a call. So this budget's only job is the truly unbounded case -- an
 * unattended script with nobody to interrupt it.
 *
 * That distinction sets the number. Sized to bound the WORST case rather than
 * to keep the typical one snappy, because a budget tight enough to feel fast
 * converts correct answers into refusals: at 2M nodes a dense n=100 graph
 * refused after 14s where the unbudgeted search had ANSWERED it in 30s. Trading
 * a correct result for a guard TimeConstrained already provides better is the
 * wrong trade.
 *
 * Calibration: the search sustains roughly 90k-195k nodes/second (n=80 dense,
 * 64k nodes, 0.33s; n=100 dense, 3.9M nodes, 30s; n=128 dense, 90k nodes/s).
 * So 8M nodes is ~1.5-2 minutes at the ceiling. Dense n=100 (3.9M nodes) is
 * comfortably inside it and answers again; every case in the default test suite
 * finishes under 100k nodes, so the budget is never why a documented case
 * fails. Bound an interactive call with TimeConstrained, not by lowering this. */
#define FVC_MAX_STEPS 8000000L

/* ---- Undirected adjacency, walked in place -------------------------------- */

/* True iff u and v are joined by an edge in either direction. O(deg(u)); at
 * n <= FVC_MAX_VERTICES this is cheaper than materialising a union structure
 * that graph_adj_free would not know to release. */
static int fvc_adjacent(const GraphAdj* a, int u, int v) {
    for (int j = 0; j < a->outdeg[u]; j++) if (a->out[u][j] == v) return 1;
    for (int j = 0; j < a->indeg[u];  j++) if (a->in[u][j]  == v) return 1;
    return 0;
}

/* Undirected degree of v. Counts a mutual pair twice (an UndirectedEdge appears
 * in both out[] and in[]), which is fine: it is used only as an ordering
 * heuristic and a tie-break, never as a graph-theoretic degree. */
static int fvc_degree(const GraphAdj* a, int v) {
    return a->outdeg[v] + a->indeg[v];
}

/* Vertex indices ordered by descending fvc_degree, ties by ascending index --
 * deterministic, so repeated calls give identical colourings. Insertion sort:
 * n <= FVC_MAX_VERTICES. */
static void fvc_order_by_degree(const GraphAdj* a, int* order) {
    int n = a->n;
    for (int i = 0; i < n; i++) order[i] = i;
    for (int i = 1; i < n; i++) {
        int v = order[i], dv = fvc_degree(a, v), j = i - 1;
        while (j >= 0 && fvc_degree(a, order[j]) < dv) { order[j + 1] = order[j]; j--; }
        order[j + 1] = v;
    }
}

/* ---- Upper bound: DSATUR -------------------------------------------------- */

int fvc_dsatur_bound(const GraphAdj* a, int* colour) {
    int n = a->n;
    if (n <= 0) return 0;

    for (int i = 0; i < n; i++) colour[i] = 0;   /* 0 = uncoloured */

    /* Scratch reused across iterations: which colours a vertex's neighbours
     * already carry. Sized n+2 because no colouring of n vertices uses more
     * than n colours. */
    char* nbcol = calloc((size_t)n + 2, sizeof(char));
    if (!nbcol) return 0;

    int used = 0;
    for (int step = 0; step < n; step++) {
        /* Pick the uncoloured vertex of maximum saturation (distinct colours
         * among its neighbours), breaking ties on degree then on index. */
        int best = -1, best_sat = -1, best_deg = -1;
        for (int v = 0; v < n; v++) {
            if (colour[v]) continue;
            for (int c = 0; c <= n + 1; c++) nbcol[c] = 0;
            int sat = 0;
            for (int j = 0; j < a->outdeg[v]; j++) {
                int c = colour[a->out[v][j]];
                if (c && !nbcol[c]) { nbcol[c] = 1; sat++; }
            }
            for (int j = 0; j < a->indeg[v]; j++) {
                int c = colour[a->in[v][j]];
                if (c && !nbcol[c]) { nbcol[c] = 1; sat++; }
            }
            int deg = fvc_degree(a, v);
            if (sat > best_sat || (sat == best_sat && deg > best_deg)) {
                best = v; best_sat = sat; best_deg = deg;
            }
        }

        /* Smallest colour no neighbour of `best` carries. */
        for (int c = 0; c <= n + 1; c++) nbcol[c] = 0;
        for (int j = 0; j < a->outdeg[best]; j++) {
            int c = colour[a->out[best][j]];
            if (c) nbcol[c] = 1;
        }
        for (int j = 0; j < a->indeg[best]; j++) {
            int c = colour[a->in[best][j]];
            if (c) nbcol[c] = 1;
        }
        int c = 1;
        while (c <= n + 1 && nbcol[c]) c++;
        colour[best] = c;
        if (c > used) used = c;
    }

    free(nbcol);
    return used;
}

/* ---- Lower bound: greedy clique ------------------------------------------- */

int fvc_clique_bound(const GraphAdj* a) {
    int n = a->n;
    if (n <= 0) return 0;

    int* order  = calloc((size_t)n, sizeof(int));
    int* clique = calloc((size_t)n, sizeof(int));
    if (!order || !clique) { free(order); free(clique); return n > 0 ? 1 : 0; }
    fvc_order_by_degree(a, order);

    /* Multi-start: grow a greedy clique from each vertex in turn and keep the
     * largest. O(n^3) worst case, which at n <= 128 is a few million integer
     * comparisons -- negligible next to the search it prunes, and much tighter
     * than a single greedy pass. */
    int best = 1;
    for (int s = 0; s < n; s++) {
        int size = 0;
        clique[size++] = order[s];
        for (int i = 0; i < n; i++) {
            int v = order[i];
            if (v == order[s]) continue;
            int ok = 1;
            for (int j = 0; j < size && ok; j++)
                if (!fvc_adjacent(a, v, clique[j])) ok = 0;
            if (ok) clique[size++] = v;
        }
        if (size > best) best = size;
    }

    free(order); free(clique);
    return best;
}

/* ---- Exact search -------------------------------------------------------- */

/* Branch-and-bound state, threaded through the recursion by pointer so the
 * recursive frame stays small (n can be FVC_MAX_VERTICES deep). */
typedef struct {
    const GraphAdj* a;
    int* col;        /* working colouring, 1-based; 0 = uncoloured           */
    int* best;       /* best complete colouring found so far                 */
    int  best_k;     /* colours used by *best -- the incumbent upper bound    */
    int  lb;         /* proven lower bound; search may stop once best_k == lb */
    long steps;
    int  aborted;    /* set when steps exceeded FVC_MAX_STEPS: best_k is then
                      * only an upper bound, NOT proven minimal              */
} FvcBB;

/* DSATUR branch-and-bound. Selects the next vertex DYNAMICALLY as the
 * uncoloured vertex of maximum saturation, which is far stronger than any
 * static order: it drives the search into the most constrained region first,
 * so a dead end is reached in a few levels rather than near the leaves.
 *
 * Three prunings. (1) BOUND: a partial colouring already using best_k colours
 * cannot beat the incumbent, so cut. (2) SYMMETRY BREAKING: colours are
 * interchangeable labels, so a vertex need only try colours 1..used+1 --
 * opening colour used+2 first merely renames a branch already explored.
 * (3) OPTIMALITY: once best_k == lb the incumbent is provably optimal.
 *
 * `used` is the number of distinct colours in the current partial colouring.
 * Unlike an iterative deepening over k, this improves the incumbent as it goes,
 * so it never has to refute each smaller k from scratch. */
static void fvc_bb(FvcBB* s, int ncoloured, int used) {
    const GraphAdj* a = s->a;
    int n = a->n;

    if (s->aborted) return;                      /* budget spent; unwind     */
    if (s->best_k <= s->lb) return;              /* already provably optimal */
    if (used >= s->best_k) return;               /* bound: cannot improve    */

    if (ncoloured == n) {
        for (int i = 0; i < n; i++) s->best[i] = s->col[i];
        s->best_k = used;
        return;
    }

    if (++s->steps > FVC_MAX_STEPS) { s->aborted = 1; return; }

    /* Cooperate with TimeConstrained. The evaluator normally polls this once
     * per rewrite step, which never happens inside a single long-running
     * builtin -- so a search like this one must poll for itself, exactly as
     * the facility's contract intends (declared for this purpose in core.h).
     *
     * On hosts with working SIGPROF/ITIMER_PROF the signal interrupts us
     * anyway; this is the backstop for hosts where it does not (WSL 1). Note
     * it may siglongjmp out of this recursion and therefore leak this search's
     * scratch buffers -- the same tradeoff every abortable builtin in the tree
     * makes (there is no cleanup registry), not something specific here.
     *
     * Polled every 4096 nodes: often enough that a deadline is honoured
     * promptly, rare enough to stay off the hot path. */
    if ((s->steps & 0xFFF) == 0) tc_check_deadline();

    /* Select the uncoloured vertex of maximum saturation, breaking ties on
     * degree then index, so the choice is deterministic. */
    int v = -1, best_sat = -1, best_deg = -1;
    char seen[FVC_MAX_VERTICES + 2];
    for (int u = 0; u < n; u++) {
        if (s->col[u]) continue;
        for (int c = 0; c <= used + 1; c++) seen[c] = 0;
        int sat = 0;
        for (int j = 0; j < a->outdeg[u]; j++) {
            int c = s->col[a->out[u][j]];
            if (c && !seen[c]) { seen[c] = 1; sat++; }
        }
        for (int j = 0; j < a->indeg[u]; j++) {
            int c = s->col[a->in[u][j]];
            if (c && !seen[c]) { seen[c] = 1; sat++; }
        }
        int deg = fvc_degree(a, u);
        if (sat > best_sat || (sat == best_sat && deg > best_deg)) {
            v = u; best_sat = sat; best_deg = deg;
        }
    }

    /* One walk of v's undirected neighbourhood rules out every colour it
     * conflicts with. Stack-sized by the cap: no allocation in the hot path. */
    char forbid[FVC_MAX_VERTICES + 2];
    int limit = (used + 1 < s->best_k - 1) ? used + 1 : s->best_k - 1;
    for (int c = 0; c <= limit + 1; c++) forbid[c] = 0;
    for (int j = 0; j < a->outdeg[v]; j++) {
        int c = s->col[a->out[v][j]];
        if (c > 0 && c <= limit) forbid[c] = 1;
    }
    for (int j = 0; j < a->indeg[v]; j++) {
        int c = s->col[a->in[v][j]];
        if (c > 0 && c <= limit) forbid[c] = 1;
    }

    for (int c = 1; c <= limit; c++) {
        if (forbid[c]) continue;
        s->col[v] = c;
        fvc_bb(s, ncoloured + 1, (c > used) ? c : used);
        if (s->aborted || s->best_k <= s->lb) { s->col[v] = 0; return; }
    }
    s->col[v] = 0;
}

int fvc_search(const GraphAdj* a, int* colour, long* steps_out) {
    int n = a->n;
    if (steps_out) *steps_out = 0;
    if (n <= 0) return 0;

    int ub = fvc_dsatur_bound(a, colour);
    if (ub <= 0) return 0;                  /* allocation failure inside DSATUR */
    int lb = fvc_clique_bound(a);

    /* The clique is a subgraph needing lb distinct colours and DSATUR exhibits
     * a colouring using ub, so lb <= chi <= ub. When they meet, DSATUR's
     * colouring is already optimal -- answer with zero search steps. This is
     * what makes CompleteGraph[128] immediate rather than a hang. */
    if (lb >= ub) return ub;

    int* work = calloc((size_t)n, sizeof(int));
    int* best = calloc((size_t)n, sizeof(int));
    /* Returning `ub` here would hand back DSATUR's colouring as if it were
     * proven minimal -- precisely the silent wrong answer this head exists to
     * avoid. Give up instead; the caller leaves the expression unevaluated. */
    if (!work || !best) { free(work); free(best); return 0; }
    for (int i = 0; i < n; i++) best[i] = colour[i];   /* DSATUR is the incumbent */

    FvcBB s;
    s.a = a; s.col = work; s.best = best; s.best_k = ub; s.lb = lb;
    s.steps = 0; s.aborted = 0;
    fvc_bb(&s, 0, 0);

    /* An aborted search has an incumbent but no proof it is minimal, and a
     * valid-but-larger colouring is exactly the silent contradiction of the
     * documented semantics that motivated the exact search. So report failure
     * rather than the incumbent, and let the head return unevaluated. */
    int result = s.aborted ? 0 : s.best_k;
    if (result > 0) for (int i = 0; i < n; i++) colour[i] = best[i];

    free(work); free(best);
    if (steps_out) *steps_out = s.steps;
    return result;
}

/* ---- The head ------------------------------------------------------------- */

Expr* builtin_find_vertex_coloring(Expr* res) {
    /* Form 1 only. Forms 2 and 3 (FindVertexColoring[g, {c1, ...}] and
     * [g, l]) are a later mapping layer over this one, so a second argument
     * must leave the expression unevaluated rather than be silently ignored. */
    if (res->data.function.arg_count != 1) return NULL;

    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;                       /* not a valid graph            */

    size_t n = (size_t)a->n;
    if (n == 0) {                              /* the empty graph colours to {} */
        graph_adj_free(a);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    if (a->n > FVC_MAX_VERTICES) { graph_adj_free(a); return NULL; }

    int* colour = calloc(n, sizeof(int));
    if (!colour) { graph_adj_free(a); return NULL; }

    /* Zero means the search could not PROVE minimality (budget spent, or an
     * allocation failure inside it). Refusing is the whole point: returning the
     * incumbent would be a valid colouring that silently contradicts the
     * documented semantics. */
    if (fvc_search(a, colour, NULL) <= 0) {
        free(colour); graph_adj_free(a); return NULL;
    }

    Expr** elems = calloc(n, sizeof(Expr*));
    if (!elems) { free(colour); graph_adj_free(a); return NULL; }
    for (size_t i = 0; i < n; i++) elems[i] = expr_new_integer((int64_t)colour[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);   /* expr_new_function memcpys -- src/expr.c:257 */

    free(colour);
    graph_adj_free(a);
    return out;
}
