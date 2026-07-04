/* kneser.c - KneserGraph[n, k]: the Kneser graph K(n, k). Its vertices are the
 * k-element subsets of {1, ..., n}, and two are adjacent iff the subsets are
 * disjoint.
 *
 * Vertices are C(n, k) subsets, each labelled by the sorted List of its elements;
 * disjointness is tested in O(1) with element bitmasks (so n is capped at 62).
 * O(C(n,k)^2). Special cases: K(n, 1) is the complete graph K_n, K(5, 2) is the
 * Petersen graph (10 vertices, 15 edges, 3-regular), K(2k, k) is a perfect
 * matching.
 *
 * n >= 0, 0 <= k <= n integers; left unevaluated if the vertex count would exceed
 * a safety cap.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on bad
 * arguments or an oversized graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

#define KNESER_MAX_VERTICES 3000

Expr* builtin_kneser_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* ne = res->data.function.args[0];
    const Expr* ke = res->data.function.args[1];
    if (ne->type != EXPR_INTEGER || ke->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer, k = (long)ke->data.integer;
    if (n < 0 || k < 0 || k > n || n > 62) return NULL;

    /* Vertex count C(n, k), with overflow/size guard. */
    long C = 1;
    for (long i = 0; i < k; i++) {
        C = C * (n - i) / (i + 1);
        if (C > KNESER_MAX_VERTICES) return NULL;
    }

    /* Enumerate k-subsets in lexicographic order: labels + element bitmasks. */
    Expr** vc = malloc((size_t)(C > 0 ? C : 1) * sizeof(Expr*));
    unsigned long long* mask = malloc((size_t)(C > 0 ? C : 1) * sizeof(unsigned long long));
    long* idx = malloc((size_t)(k > 0 ? k : 1) * sizeof(long));
    for (long i = 0; i < k; i++) idx[i] = i;              /* 0-based elements */
    for (long v = 0; v < C; v++) {
        Expr** els = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr*));
        unsigned long long m = 0;
        for (long i = 0; i < k; i++) { els[i] = expr_new_integer(idx[i] + 1); m |= 1ULL << idx[i]; }
        vc[v] = expr_new_function(expr_new_symbol(SYM_List), els, (size_t)k);
        mask[v] = m;
        free(els);
        /* Advance to the next k-combination. */
        long p = k - 1;
        while (p >= 0 && idx[p] == n - k + p) p--;
        if (p < 0) break;
        idx[p]++;
        for (long i = p + 1; i < k; i++) idx[i] = idx[i - 1] + 1;
    }
    free(idx);

    /* Edges: disjoint subsets (mask & mask == 0), distinct vertices. */
    Expr** ec = NULL; size_t ecap = 0, me = 0;
    for (long i = 0; i < C; i++)
        for (long j = i + 1; j < C; j++)
            if ((mask[i] & mask[j]) == 0) {
                if (me == ecap) { ecap = ecap ? ecap * 2 : 16; ec = realloc(ec, ecap * sizeof(Expr*)); }
                Expr* a[2] = { expr_copy(vc[i]), expr_copy(vc[j]) };
                ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
            }
    free(mask);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)C);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
