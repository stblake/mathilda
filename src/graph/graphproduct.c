/* graphproduct.c - GraphProduct[g1, g2, type]: a product graph on the vertex set
 * V1 x V2, with the four standard adjacency rules. Product vertices are the
 * pairs {a, b} (a in g1, b in g2); the result is undirected (products are
 * defined on the underlying undirected graphs).
 *
 * For product vertices (a1,b1) != (a2,b2), with A1/A2 the adjacencies of g1/g2:
 *   "Cartesian"     joined iff (a1=a2 and b1~b2) or (b1=b2 and a1~a2)
 *   "Tensor"        joined iff a1~a2 and b1~b2                (categorical/direct)
 *   "Strong"        Cartesian OR Tensor
 *   "Lexicographic" joined iff a1~a2, or (a1=a2 and b1~b2)    (composition)
 *
 * O((n1 n2)^2) over vertex pairs -- small-graph scale. Classic identities:
 * P2 [] P2 = C4, K2 (Tensor) K2 = 2K2, K2 (Strong) K2 = K4, C4 [] K2 = the cube.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * both arguments are valid graphs and the third is a known product-type string.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

enum { PROD_CART, PROD_TENSOR, PROD_STRONG, PROD_LEX };

/* Undirected boolean adjacency for a validated graph (caller frees). */
static char* build_adj(const Expr* g, int* np) {
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    *np = n;
    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }
    return adj;
}

Expr* builtin_graph_product(Expr* res) {
    if (res->data.function.arg_count != 3) return NULL;
    const Expr* g1 = res->data.function.args[0];
    const Expr* g2 = res->data.function.args[1];
    const Expr* ty = res->data.function.args[2];
    if (!graph_is_valid(g1) || !graph_is_valid(g2)) return NULL;
    if (ty->type != EXPR_STRING) return NULL;
    int kind;
    if      (strcmp(ty->data.string, "Cartesian") == 0)     kind = PROD_CART;
    else if (strcmp(ty->data.string, "Tensor") == 0)        kind = PROD_TENSOR;
    else if (strcmp(ty->data.string, "Strong") == 0)        kind = PROD_STRONG;
    else if (strcmp(ty->data.string, "Lexicographic") == 0) kind = PROD_LEX;
    else return NULL;

    const Expr* v1 = g1->data.function.args[0];
    const Expr* v2 = g2->data.function.args[0];
    int n1, n2;
    char* a1 = build_adj(g1, &n1);
    char* a2 = build_adj(g2, &n2);
    int N = n1 * n2;

    /* Product vertices {a, b} at index i*n2 + j. */
    Expr** vc = malloc((N > 0 ? N : 1) * sizeof(Expr*));
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            Expr* pair[2] = { expr_copy(v1->data.function.args[i]),
                              expr_copy(v2->data.function.args[j]) };
            vc[i * n2 + j] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        }

    /* Edges over unordered product-vertex pairs. */
    size_t cap = (size_t)N * (N > 0 ? (size_t)(N - 1) : 0) / 2;
    Expr** ec = (cap > 0) ? malloc(cap * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (int p = 0; p < N; p++) {
        int i1 = p / n2, j1 = p % n2;
        for (int q = p + 1; q < N; q++) {
            int i2 = q / n2, j2 = q % n2;
            int ea = (i1 != i2) && a1[(size_t)i1 * n1 + i2];   /* a1 ~ a2 */
            int eb = (j1 != j2) && a2[(size_t)j1 * n2 + j2];   /* b1 ~ b2 */
            int adj;
            switch (kind) {
                case PROD_CART:   adj = (i1 == i2 && eb) || (j1 == j2 && ea); break;
                case PROD_TENSOR: adj = ea && eb; break;
                case PROD_STRONG: adj = ((i1 == i2 && eb) || (j1 == j2 && ea)) || (ea && eb); break;
                default:          adj = ea || (i1 == i2 && eb); break;   /* Lexicographic */
            }
            if (!adj) continue;
            Expr* args[2] = { expr_copy(vc[p]), expr_copy(vc[q]) };
            ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), args, 2);
        }
    }
    free(a1); free(a2);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)N);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
