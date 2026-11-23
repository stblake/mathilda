/* Identity and diagonal matrix constructors.
 *
 * Both produce rank-2 nested lists of Integers (zeros and ones for
 * IdentityMatrix; user-supplied diagonal entries for DiagonalMatrix).
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include "common.h"
#include "pack.h"
#include "sparse/sparse.h"
#include "message.h"
#include "checked_alloc.h"
#include "print.h"
#include <stdlib.h>

Expr* builtin_identitymatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* IdentityMatrix[n | {m, n}, SparseArray]: built sparse directly, one
     * stored entry per row (sparse/sparse_struct.c). */
    if (res->data.function.arg_count == 2 &&
        res->data.function.args[1]->type == EXPR_SYMBOL &&
        res->data.function.args[1]->data.symbol.name == SYM_SparseArray) {
        Expr* d = res->data.function.args[0];
        if (d->type == EXPR_INTEGER) return sp_identity_matrix(d->data.integer, d->data.integer);
        if (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL &&
            d->data.function.head->data.symbol.name == SYM_List && d->data.function.arg_count == 2 &&
            d->data.function.args[0]->type == EXPR_INTEGER && d->data.function.args[1]->type == EXPR_INTEGER)
            return sp_identity_matrix(d->data.function.args[0]->data.integer,
                                      d->data.function.args[1]->data.integer);
        return NULL;
    }
    if (res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t m = -1, n = -1;

    /* A number that is not a positive machine integer (0, -1, 2.5, 5/2, a
     * BigInt, 1 + I) is IdentityMatrix::dims; a symbol or a malformed pair
     * stays unevaluated silently. */
    bool number = arg->type == EXPR_INTEGER || arg->type == EXPR_BIGINT || arg->type == EXPR_REAL ||
                  head_is(arg, SYM_Rational) || head_is(arg, SYM_Complex);
#ifdef USE_MPFR
    number = number || arg->type == EXPR_MPFR;
#endif
    if (number && !(arg->type == EXPR_INTEGER && arg->data.integer > 0)) {
        char* s = expr_to_string_outputform(arg);
        mth_message("IdentityMatrix", "dims", "Dimension specification %s should be a positive "
                    "machine integer or a pair of positive machine integers.", s ? s : "");
        free(s);
        return NULL;
    }

    if (arg->type == EXPR_INTEGER) {
        m = arg->data.integer;
        n = m;
    } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol.name == SYM_List) {
        if (arg->data.function.arg_count == 2) {
            Expr* arg_m = arg->data.function.args[0];
            Expr* arg_n = arg->data.function.args[1];
            if (arg_m->type == EXPR_INTEGER && arg_n->type == EXPR_INTEGER) {
                m = arg_m->data.integer;
                n = arg_n->data.integer;
            }
        }
    }

    if (m <= 0 || n <= 0) return NULL;
    /* m, n are the user's: a matrix past memory, or one whose m*n overflows,
     * aborts with General::nomem instead of wrapping the buffer size
     * (audit IdentityMatrix#1).  Each entry costs at least one machine word. */
    if (!mth_alloc2_ok(m, n, sizeof(int64_t))) return mth_nomem_abort();

    /* THE BUFFER PATH. Every entry is the exact Integer 0 or 1, so the matrix
     * is uniformly int64 and no element's head can change by writing it into a
     * buffer -- an integer buffer materialises back to Integers. That is what
     * makes this legal rather than merely fast (N1/N2).
     *
     * IdentityMatrix[1000] built 10^6 Expr nodes and 1000 List wrappers for a
     * matrix that is one bit per entry: 77.3 ms against NumPy's 231 us for
     * np.eye (HPC plan item C.5). */
    {
        int64_t dims[2] = { m, n };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, NDT_INT64, &raw);
        if (packed) {
            int64_t* b = (int64_t*)raw;
            for (int64_t z = 0; z < m * n; z++) b[z] = 0;
            int64_t d = (m < n) ? m : n;
            for (int64_t i = 0; i < d; i++) b[i * n + i] = 1;
            return packed;
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int64_t i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t j = 0; j < n; j++) {
            if (i == j) {
                row_elems[j] = expr_new_integer(1);
            } else {
                row_elems[j] = expr_new_integer(0);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row_elems, n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, m);
    free(rows);
    return result;
}

Expr* builtin_diagonalmatrix(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || list->data.function.head->data.symbol.name != SYM_List) {
        return expr_copy(res);
    }

    int64_t s = list->data.function.arg_count;
    int64_t k = 0;

    if (res->data.function.arg_count >= 2) {
        Expr* k_expr = res->data.function.args[1];
        if (k_expr->type == EXPR_INTEGER) {
            k = k_expr->data.integer;
        } else {
            return expr_copy(res);
        }
    }

    int64_t m = -1, n = -1;

    if (res->data.function.arg_count == 3) {
        Expr* dim_expr = res->data.function.args[2];
        if (dim_expr->type == EXPR_INTEGER) {
            m = dim_expr->data.integer;
            n = m;
        } else if (dim_expr->type == EXPR_FUNCTION && dim_expr->data.function.head->type == EXPR_SYMBOL && dim_expr->data.function.head->data.symbol.name == SYM_List) {
            if (dim_expr->data.function.arg_count == 2) {
                Expr* arg_m = dim_expr->data.function.args[0];
                Expr* arg_n = dim_expr->data.function.args[1];
                if (arg_m->type == EXPR_INTEGER && arg_n->type == EXPR_INTEGER) {
                    m = arg_m->data.integer;
                    n = arg_n->data.integer;
                } else return expr_copy(res);
            } else return expr_copy(res);
        } else return expr_copy(res);
    } else {
        m = s + (k > 0 ? k : -k);
        n = m;
    }

    if (m < 0 || n < 0) return expr_copy(res);
    /* As IdentityMatrix: an impossible m x n aborts (audit DiagonalMatrix#1). */
    if (!mth_alloc2_ok(m, n, sizeof(int64_t))) return mth_nomem_abort();

    /* THE ZEROS TAKE THE DIAGONAL'S EXACTNESS -- see common.h's note on
     * machine-real contagion. One machine `Real` on the diagonal makes the
     * whole matrix machine-real, the invented zeros and the exact entries
     * alike, which is Mathematica's answer:
     *
     *     DiagonalMatrix[{1., 2., 3.}]  {{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}
     *     DiagonalMatrix[{1, 2, 3.}]    the same -- the Integers widen too
     *     DiagonalMatrix[{1/2, 1.}]     {{0.5, 0.}, {0., 1.}}
     *     DiagonalMatrix[{a, 1.}]       {{a, 0.}, {0., 1.}}  (symbol stays)
     *
     * This used to write the exact Integer 0 in every case, which was both
     * wrong against Mathematica and the reason a Real diagonal could not be
     * packed: a matrix of two heads fits no uniform buffer, and it cost 320x
     * against NumPy. Getting the exactness right is what makes it one dtype.
     *
     * Only MACHINE Real is contagious. An MPFR diagonal keeps exact zeros, and
     * so does an exact one -- so there are three cases, and the flags below
     * name them: all-Integer (int64 buffer), all-machine-real (float64 buffer),
     * and everything else (the boxed matrix, with `zero` set correctly). */
    bool diag_has_real = false;
    for (int64_t i = 0; i < s; i++)
        if (list->data.function.args[i]->type == EXPR_REAL) { diag_has_real = true; break; }

    bool all_int = !diag_has_real;
    for (int64_t i = 0; i < s && all_int; i++)
        if (list->data.function.args[i]->type != EXPR_INTEGER) all_int = false;

    bool all_real = diag_has_real;
    for (int64_t i = 0; i < s && all_real; i++) {
        double d;
        if (!common_machine_real_value(list->data.function.args[i], &d)) all_real = false;
    }

    if (all_int || all_real) {
        int64_t dims[2] = { m, n };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, all_int ? NDT_INT64 : NDT_FLOAT64, &raw);
        if (packed) {
            int64_t* bi = (int64_t*)raw;
            double*  bd = (double*)raw;
            for (int64_t z = 0; z < m * n; z++) { if (all_int) bi[z] = 0; else bd[z] = 0.0; }
            for (int64_t i = 0; i < m; i++) {
                int64_t j = i + k;
                if (j < 0 || j >= n) continue;
                int64_t idx = (i < j) ? i : j;
                if (idx < 0 || idx >= s) continue;
                const Expr* el = list->data.function.args[idx];
                if (all_int) bi[i * n + j] = el->data.integer;
                else         common_machine_real_value(el, &bd[i * n + j]);
            }
            return packed;
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * m);
    for (int64_t i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int64_t j = 0; j < n; j++) {
            row_elems[j] = NULL;
            if (j - i == k) {
                int64_t list_idx = (i < j) ? i : j;
                if (list_idx >= 0 && list_idx < s) {
                    Expr* el = list->data.function.args[list_idx];
                    double d;
                    /* An exact numeric entry beside a machine Real becomes one:
                     * DiagonalMatrix[{1/2, 1.}] is {{0.5, 0.}, {0., 1.}}. A
                     * symbolic (or MPFR, or Complex) entry is copied as it is,
                     * and only the zeros around it turn Real. */
                    row_elems[j] = (diag_has_real && common_machine_real_value(el, &d))
                                       ? expr_new_real(d) : expr_copy(el);
                }
            }
            if (!row_elems[j])
                row_elems[j] = diag_has_real ? expr_new_real(0.0) : expr_new_integer(0);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row_elems, n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, m);
    free(rows);
    return result;
}
