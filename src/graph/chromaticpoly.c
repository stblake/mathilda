/* chromaticpoly.c - ChromaticPolynomial[g, k]: the chromatic polynomial of g,
 * the number of proper k-colourings as a polynomial in k. With a symbolic k it
 * returns the polynomial (e.g. k(k-1)(k-2) for a triangle); with a numeric k it
 * returns the colouring count.
 *
 * Computed from the subgraph-expansion (Whitney) form
 *
 *   P(g, k) = sum_{S subset of E} (-1)^|S| * k^{c(S)}
 *
 * where c(S) is the number of connected components of the spanning subgraph
 * (V, S). Collecting subsets by component count gives integer coefficients
 * a_j = sum_{S : c(S)=j} (-1)^|S|, and the result is sum_j a_j k^j, assembled as
 * an Expr and reduced by the evaluator -- so a symbolic k yields the polynomial
 * and a numeric k yields the integer through one code path. Edge direction is
 * ignored.
 *
 * The sum is over 2^|E| subsets, so it is exponential in the edge count; the
 * call is left unevaluated when |E| exceeds a modest bound to avoid runaway
 * work (small graphs -- the interactive/gallery case -- are unaffected).
 *
 * Memory (SPEC section 4): returns a fresh Expr; frees res. NULL on a non-graph
 * or an oversized graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

#define CHROM_MAX_EDGES 24   /* 2^24 subsets ~ 16.7M: the practical ceiling */

static int uf_find(int* p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

Expr* builtin_chromatic_polynomial(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* var = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n  = (int)verts->data.function.arg_count;
    int m  = (int)edges->data.function.arg_count;
    if (m > CHROM_MAX_EDGES) return NULL;    /* too large; stay unevaluated */

    /* Edge endpoints as vertex indices. */
    int* ea = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    int* eb = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    for (int k = 0; k < m; k++) {
        const Expr* ed = edges->data.function.args[k];
        ea[k] = graph_vertex_index(verts, ed->data.function.args[0]);
        eb[k] = graph_vertex_index(verts, ed->data.function.args[1]);
    }

    /* coeff[j] = sum over subsets S with c(S)=j of (-1)^|S|, j = 0..n. */
    long* coeff = calloc((size_t)(n + 1), sizeof(long));
    int* parent = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    unsigned long subsets = 1UL << m;
    for (unsigned long s = 0; s < subsets; s++) {
        for (int i = 0; i < n; i++) parent[i] = i;
        int comps = n, bits = 0;
        for (int k = 0; k < m; k++) {
            if (!(s & (1UL << k))) continue;
            bits++;
            int ra = uf_find(parent, ea[k]), rb = uf_find(parent, eb[k]);
            if (ra != rb) { parent[ra] = rb; comps--; }
        }
        coeff[comps] += (bits & 1) ? -1 : 1;
    }
    free(ea); free(eb); free(parent);

    /* Build sum_j coeff[j] * var^j and let the evaluator reduce it. */
    Expr** terms = calloc((size_t)(n + 1), sizeof(Expr*));
    int nt = 0;
    for (int j = 0; j <= n; j++) {
        if (coeff[j] == 0) continue;
        Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy((Expr*)var), expr_new_integer(j) }, 2);
        terms[nt++] = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(coeff[j]), pw }, 2);
    }
    free(coeff);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt);
    free(terms);
    return evaluate(sum);
}
