/* Generalised cross product.
 *
 * Cross[v1, ..., v(n-1)] for n-dimensional vectors v_i returns the vector
 * orthogonal to all inputs.  Each component is the signed determinant of
 * the minor obtained by deleting the i-th column from the matrix of
 * stacked rows.
 */

#include "linalg.h"
#include "eval.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>

Expr* builtin_cross(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t m = res->data.function.arg_count;
    if (m == 0) return NULL;

    size_t n = m + 1;
    bool valid = true;
    for (size_t i = 0; i < m; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || arg->data.function.head->type != EXPR_SYMBOL || arg->data.function.head->data.symbol != SYM_List) {
            valid = false;
            break;
        }
        if (arg->data.function.arg_count != n) {
            valid = false;
            break;
        }
    }

    if (!valid) {
        fprintf(stderr, "Cross::nonn1: The arguments are expected to be vectors of equal length, and the number of arguments is expected to be 1 less than their length.\n");
        return NULL;
    }

    Expr** result_args = malloc(sizeof(Expr*) * n);

    for (size_t i = 0; i < n; i++) {
        Expr** minor_flat = malloc(sizeof(Expr*) * m * m);
        for (size_t r = 0; r < m; r++) {
            Expr* vec = res->data.function.args[r];
            size_t c_idx = 0;
            for (size_t c = 0; c < n; c++) {
                if (c == i) continue;
                minor_flat[r * m + c_idx] = vec->data.function.args[c];
                c_idx++;
            }
        }

        int* cols = malloc(sizeof(int) * m);
        for (size_t c = 0; c < m; c++) cols[c] = (int)c;

        Expr* det_val = laplace_det(minor_flat, (int)m, (int)m, 0, cols);
        free(cols);
        free(minor_flat);

        int sign = ((m + i) % 2 == 0) ? 1 : -1;
        if (sign == -1) {
            Expr* t_args[2] = { expr_new_integer(-1), det_val };
            result_args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
        } else {
            result_args[i] = det_val;
        }
    }

    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), result_args, n);
    free(result_args);
    return final_res;
}
