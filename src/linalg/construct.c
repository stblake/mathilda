/* Identity and diagonal matrix constructors.
 *
 * Both produce rank-2 nested lists of Integers (zeros and ones for
 * IdentityMatrix; user-supplied diagonal entries for DiagonalMatrix).
 */

#include "linalg.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_identitymatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t m = -1, n = -1;

    if (arg->type == EXPR_INTEGER) {
        m = arg->data.integer;
        n = m;
    } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol == SYM_List) {
        if (arg->data.function.arg_count == 2) {
            Expr* arg_m = arg->data.function.args[0];
            Expr* arg_n = arg->data.function.args[1];
            if (arg_m->type == EXPR_INTEGER && arg_n->type == EXPR_INTEGER) {
                m = arg_m->data.integer;
                n = arg_n->data.integer;
            }
        }
    }

    if (m < 0 || n < 0) return expr_copy(res);

    Expr** rows = malloc(sizeof(Expr*) * m);
    for (int i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
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
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || list->data.function.head->data.symbol != SYM_List) {
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
        } else if (dim_expr->type == EXPR_FUNCTION && dim_expr->data.function.head->type == EXPR_SYMBOL && dim_expr->data.function.head->data.symbol == SYM_List) {
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

    Expr** rows = malloc(sizeof(Expr*) * m);
    for (int64_t i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int64_t j = 0; j < n; j++) {
            if (j - i == k) {
                int64_t list_idx = (i < j) ? i : j;
                if (list_idx >= 0 && list_idx < s) {
                    row_elems[j] = expr_copy(list->data.function.args[list_idx]);
                } else {
                    row_elems[j] = expr_new_integer(0);
                }
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
