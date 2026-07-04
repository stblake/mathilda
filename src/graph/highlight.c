/* highlight.c - HighlightGraph[g, parts]: draw g with selected vertices/edges
 * emphasized.
 *
 * Returns a styled Graphics[...] (like GraphPlot) in which highlighted elements
 * use the accent color and everything else is dimmed, so it auto-displays as a
 * diagram. `parts` is a list whose elements may be:
 *   - a vertex               -> that vertex is highlighted;
 *   - an edge (u<->v, u->v, Directed/UndirectedEdge[u,v]) -> that edge;
 *   - a list of vertices     -> a path: its vertices AND the edges joining
 *                               consecutive vertices are highlighted.
 * Edge endpoints are matched unordered, so direction need not be specified.
 *
 * We deliberately return a Graphics rather than a Graph: the canonical
 * Graph[List,List] form is locked to simple graphs (no annotations), so the
 * highlight lives only in the rendered picture. Layout/styling options accepted
 * by GraphPlot (GraphLayout, VertexStyle, EdgeStyle, VertexLabels, VertexSize)
 * pass through unchanged.
 *
 * Memory (SPEC section 4): returns a fresh Graphics tree; the evaluator frees
 * res. Read-only over the input graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Index of the graph edge (in canonical order) whose endpoints equal {a,b}
 * as an unordered pair, or -1. */
static int find_edge(const Expr* edges, const Expr* verts, const Expr* a,
                     const Expr* b) {
    (void)verts;
    size_t ne = edges->data.function.arg_count;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) continue;
        const Expr* u = e->data.function.args[0];
        const Expr* v = e->data.function.args[1];
        if ((expr_eq((Expr*)u, (Expr*)a) && expr_eq((Expr*)v, (Expr*)b))
            || (expr_eq((Expr*)u, (Expr*)b) && expr_eq((Expr*)v, (Expr*)a)))
            return (int)k;
    }
    return -1;
}

/* True if `e` looks like a 2-endpoint edge spec (any of the four heads). */
static int is_edge_spec(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2
        || !e->data.function.head || e->data.function.head->type != EXPR_SYMBOL)
        return 0;
    const char* h = e->data.function.head->data.symbol;
    return h == SYM_DirectedEdge || h == SYM_UndirectedEdge
        || h == SYM_Rule || h == SYM_TwoWayRule;
}

static void mark_vertex(const Expr* verts, char* hv, const Expr* v) {
    int idx = graph_vertex_index(verts, v);
    if (idx >= 0) hv[idx] = 1;
}

static void mark_edge(const Expr* edges, const Expr* verts, char* he,
                      const Expr* a, const Expr* b) {
    int idx = find_edge(edges, verts, a, b);
    if (idx >= 0) he[idx] = 1;
}

Expr* builtin_highlight_graph(Expr* res) {
    if (res->data.function.arg_count < 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* parts = res->data.function.args[1];
    if (parts->type != EXPR_FUNCTION || !parts->data.function.head
        || parts->data.function.head->type != EXPR_SYMBOL
        || parts->data.function.head->data.symbol != SYM_List)
        return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* hv = calloc((size_t)(n > 0 ? n : 1), 1);
    char* he = calloc((ne > 0 ? ne : 1), 1);
    if (!hv || !he) { free(hv); free(he); return NULL; }

    for (size_t i = 0; i < parts->data.function.arg_count; i++) {
        const Expr* part = parts->data.function.args[i];
        if (is_edge_spec(part)) {
            mark_edge(edges, verts, he, part->data.function.args[0],
                      part->data.function.args[1]);
        } else if (part->type == EXPR_FUNCTION && part->data.function.head
                   && part->data.function.head->type == EXPR_SYMBOL
                   && part->data.function.head->data.symbol == SYM_List) {
            /* A path: highlight its vertices and the joining edges. */
            size_t m = part->data.function.arg_count;
            for (size_t j = 0; j < m; j++)
                mark_vertex(verts, hv, part->data.function.args[j]);
            for (size_t j = 1; j < m; j++)
                mark_edge(edges, verts, he, part->data.function.args[j - 1],
                          part->data.function.args[j]);
        } else {
            mark_vertex(verts, hv, part);
        }
    }

    /* Reuse GraphPlot's option parsing for any trailing GraphLayout/style opts
     * by carrying them into a GraphStyle here. */
    GraphStyle st = {0};
    st.show_labels = 1;
    st.hi_vert = hv;
    st.hi_edge = he;
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        const Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2
            && opt->data.function.head
            && opt->data.function.head->type == EXPR_SYMBOL
            && opt->data.function.head->data.symbol == SYM_Rule
            && opt->data.function.args[0]->type == EXPR_SYMBOL
            && opt->data.function.args[0]->data.symbol == SYM_GraphLayout
            && opt->data.function.args[1]->type == EXPR_STRING) {
            st.layout = opt->data.function.args[1]->data.string;
        }
    }

    Expr* out = graph_render(g, &st);
    free(hv); free(he);
    return out;
}
