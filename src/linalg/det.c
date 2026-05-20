/* Determinant by Laplace expansion.
 *
 * laplace_det is also used by builtin_cross (in cross.c) to evaluate the
 * minors of the generalised cross product, so it is exported via linalg.h.
 */

#include "linalg.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>

Expr* laplace_det(Expr** flat, int original_n, int n, int row, int* cols) {
    if (n == 1) {
        return expr_copy(flat[row * original_n + cols[0]]);
    }
    Expr** sum_args = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        int* new_cols = malloc(sizeof(int) * (n - 1));
        for (int j = 0, k = 0; j < n; j++) {
            if (j != i) new_cols[k++] = cols[j];
        }
        Expr* minor_det = laplace_det(flat, original_n, n - 1, row + 1, new_cols);
        free(new_cols);
        Expr* elem = flat[row * original_n + cols[i]];
        if (i % 2 != 0) {
            Expr* t_args[3] = { expr_new_integer(-1), expr_copy(elem), minor_det };
            sum_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 3));
        } else {
            Expr* t_args[2] = { expr_copy(elem), minor_det };
            sum_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 2));
        }
    }
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), sum_args, n));
    free(sum_args);
    return res;
}

Expr* builtin_det(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);

    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        char* arg_str = expr_to_string_fullform(arg);
        fprintf(stderr, "Det::matsq: Argument %s at position 1 is not a non-empty square matrix.\n", arg_str);
        free(arg_str);
        return NULL;
    }

    int n = (int)dims[0];
    Expr** flat = malloc(sizeof(Expr*) * n * n);
    size_t idx = 0;
    flatten_tensor(arg, flat, &idx);

    int* cols = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) cols[i] = i;

    Expr* det_val = laplace_det(flat, n, n, 0, cols);

    free(cols);
    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);

    return det_val;
}
