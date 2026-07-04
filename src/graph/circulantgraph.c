/* circulantgraph.c - CirculantGraph[n, {j1, j2, ...}] (or CirculantGraph[n, j]):
 * the circulant graph on n vertices 1..n in which vertex i is joined to
 * i +/- j (mod n) for each jump j in the list.
 *
 * Every vertex has the same neighbourhood pattern, so the graph is vertex-
 * transitive and regular. Special cases: CirculantGraph[n, {1}] is the cycle
 * C_n; CirculantGraph[n, {1,...,floor(n/2)}] is the complete graph K_n; a jump of
 * exactly n/2 contributes a single (matching) edge per vertex rather than two.
 * An n x n boolean adjacency dedups the +j / -j and cross-vertex coincidences.
 * O(n * (#jumps) + n^2).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * n >= 1 is an integer and the jumps are an integer (list).
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_circulant_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* ne = res->data.function.args[0];
    const Expr* js = res->data.function.args[1];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 1) return NULL;

    /* Collect jumps: a single integer or a List of integers. */
    const Expr* const* jump_src;
    size_t njumps;
    const Expr* single[1];
    if (js->type == EXPR_INTEGER) {
        single[0] = js; jump_src = single; njumps = 1;
    } else if (graph_is_list(js)) {
        for (size_t i = 0; i < js->data.function.arg_count; i++)
            if (js->data.function.args[i]->type != EXPR_INTEGER) return NULL;
        jump_src = (const Expr* const*)js->data.function.args;
        njumps = js->data.function.arg_count;
    } else {
        return NULL;
    }

    char* adj = calloc((size_t)n * (size_t)n, 1);
    for (long i = 0; i < n; i++)
        for (size_t t = 0; t < njumps; t++) {
            long j = (long)jump_src[t]->data.integer % n;
            if (j < 0) j += n;
            long a = (i + j) % n, b = (i - j % n + n) % n;
            if (a != i) adj[(size_t)i * n + a] = adj[(size_t)a * n + i] = 1;
            if (b != i) adj[(size_t)i * n + b] = adj[(size_t)b * n + i] = 1;
        }

    Expr** vc = malloc((size_t)n * sizeof(Expr*));
    for (long i = 0; i < n; i++) vc[i] = expr_new_integer(i + 1);

    size_t cap = (size_t)n * (n > 0 ? (size_t)(n - 1) : 0) / 2;
    Expr** ec = (cap > 0) ? malloc(cap * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (long i = 0; i < n; i++)
        for (long k = i + 1; k < n; k++)
            if (adj[(size_t)i * n + k]) {
                Expr* args[2] = { expr_new_integer(i + 1), expr_new_integer(k + 1) };
                ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), args, 2);
            }
    free(adj);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
