/* Diagonal[m] / Diagonal[m, k] -- extract the k-th diagonal of a matrix.
 *
 * Diagonal[m] gives the leading diagonal {m[[1,1]], m[[2,2]], ...} (length
 * Min[rows, cols], so it works for a non-square m). Diagonal[m, k] gives the
 * k-th diagonal: k > 0 above the leading diagonal, k < 0 below. An in-range but
 * empty diagonal (|k| beyond the matrix) gives {}.
 *
 * A machine-precision matrix (a packed List or a visible NDArray) takes the
 * buffer fast path in ndstruct_diagonal (rank-2 buffer straight to rank-1). The
 * generic path below walks the nested List, so it also handles any rank >= 2:
 * the diagonal of a rank-n tensor is rank-(n-1), since each m[[i, i+k]] is itself
 * an (n-1)-tensor and is copied through verbatim.
 */

#include "linalg.h"
#include "linkedlist.h"
#include "ndstruct.h"
#include "ndarray.h"
#include "sym_names.h"
#include <stdlib.h>
#include <stdbool.h>

/* Is `e` a List[...]? */
static bool diag_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Diagonal of a matrix with a LinkedList level -- a LinkedList of rows, or a
 * List with LinkedList rows: the List walk below over element views (a
 * LinkedList is the List of its elements). A LinkedList of rows gives a
 * LinkedList, a List of rows a List. */
static Expr* diagonal_ll(const Expr* m, int64_t k) {
    int64_t nrows = m->type == EXPR_LINKEDLIST ? ll_length(m->data.llist)
                                               : (int64_t)m->data.function.arg_count;
    if (nrows == 0) return NULL;
    {
        bool own;
        Expr* r0 = ll_seq_view(m, 0, &own);
        bool ok = r0 && ll_is_seq(r0);
        if (r0) ll_elem_done(r0, own);
        if (!ok) return NULL;
    }
    int64_t start_row, start_col;
    if (k >= 0)          { start_row = 0;     start_col = k; }
    else if (k > -nrows) { start_row = -k;    start_col = 0; }
    else                 { start_row = nrows; start_col = 0; }
    LinkedList* out = ll_new();
    for (int64_t t = 0; ; t++) {
        int64_t i = start_row + t, j = start_col + t;
        if (i >= nrows) break;
        bool own;
        Expr* row = ll_seq_view(m, i, &own);
        int64_t len = !ll_is_seq(row) ? -1
                    : row->type == EXPR_LINKEDLIST ? ll_length(row->data.llist)
                    : (int64_t)row->data.function.arg_count;
        if (len < 0 || j >= len) { ll_elem_done(row, own); break; }
        ll_push_back(out, row->type == EXPR_LINKEDLIST ? ll_get(row->data.llist, j)
                                                       : expr_copy(row->data.function.args[j]));
        ll_elem_done(row, own);
    }
    Expr* r = expr_new_linkedlist(out);
    if (m->type != EXPR_LINKEDLIST) {
        /* the List of rows answers with a List */
        Expr* l = ll_call_with_elements(expr_new_symbol(SYM_List), out);
        expr_free(r);
        return l;
    }
    return r;
}

Expr* builtin_diagonal(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* m = res->data.function.args[0];

    /* Machine-array fast path (packed List or visible NDArray): rank-2 buffer
     * straight to a rank-1 buffer; see src/ndstruct.c. */
    if (is_ndarray(m)) return ndstruct_diagonal(res);

    int64_t k = 0;
    if (argc == 2) {
        Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) return NULL;
        k = ke->data.integer;
    }

    if (ll_has_ll_level1(m)) return diagonal_ll(m, k);

    /* Require a rank >= 2 structure (a List whose rows are Lists); a vector or a
     * non-list is not a matrix, so leave the call unevaluated as Mathematica does. */
    if (!diag_is_list(m)) return NULL;
    int64_t nrows = (int64_t)m->data.function.arg_count;
    if (nrows == 0 || !diag_is_list(m->data.function.args[0])) return NULL;

    int64_t start_row, start_col;
    if (k >= 0)          { start_row = 0;     start_col = k; }
    else if (k > -nrows) { start_row = -k;    start_col = 0; }
    else                 { start_row = nrows; start_col = 0; }  /* below matrix: loop is empty */

    size_t cap = 16, n = 0;
    Expr** elems = malloc(sizeof(Expr*) * cap);
    for (int64_t t = 0; ; t++) {
        int64_t i = start_row + t, j = start_col + t;
        if (i >= nrows) break;
        Expr* row = m->data.function.args[i];
        if (!diag_is_list(row) || j >= (int64_t)row->data.function.arg_count) break;
        if (n == cap) { cap *= 2; elems = realloc(elems, sizeof(Expr*) * cap); }
        elems[n++] = expr_copy(row->data.function.args[j]);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return result;
}
