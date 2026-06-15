/* Vector and matrix norms.
 *
 *   Norm[x]       scalar  -> Abs[x]
 *   Norm[v]       vector  -> 2-norm
 *   Norm[v, p]    vector  -> p-norm, also Infinity
 *   Norm[m, "Frobenius"]  -> Frobenius norm (works on any tensor rank)
 *
 * Other matrix norms (SVD-based 2-norm, etc.) are not yet implemented and
 * fall through to NULL so the call stays symbolic.
 */

#include "linalg.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

Expr* builtin_norm(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* p = NULL;
    if (res->data.function.arg_count == 2) {
        p = res->data.function.args[1];
    }

    int64_t dims[64];
    int rank = get_tensor_dims(expr, dims);
    if (rank < 0) return NULL; // jagged array

    if (rank == 0) {
        // Scalar: Norm[x] -> Abs[x]
        if (!p) {
            Expr* args[1] = { expr_copy(expr) };
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), args, 1));
        }
        return NULL;
    }

    if (rank == 1 || (rank >= 2 && p && p->type == EXPR_STRING && strcmp(p->data.string, "Frobenius") == 0)) {
        int64_t N = 1;
        for (int i = 0; i < rank; i++) N *= dims[i];

        Expr** flat = NULL;
        if (N > 0) {
            flat = malloc(sizeof(Expr*) * N);
            size_t idx = 0;
            flatten_tensor(expr, flat, &idx);
        }

        Expr* result = NULL;

        if (p && p->type == EXPR_SYMBOL && p->data.symbol == SYM_Infinity) {
            if (N == 0) {
                result = expr_new_integer(0);
            } else {
                Expr** max_args = malloc(sizeof(Expr*) * N);
                for (int64_t i = 0; i < N; i++) {
                    max_args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(flat[i])}, 1));
                }
                result = eval_and_free(expr_new_function(expr_new_symbol(SYM_Max), max_args, N));
                free(max_args);
            }
        } else {
            Expr* norm_p = NULL;
            if (!p || (p->type == EXPR_STRING && strcmp(p->data.string, "Frobenius") == 0)) {
                norm_p = expr_new_integer(2);
            } else {
                norm_p = expr_copy(p);
            }

            if (N == 0) {
                result = expr_new_integer(0);
            } else {
                Expr** plus_args = malloc(sizeof(Expr*) * N);
                for (int64_t i = 0; i < N; i++) {
                    Expr* abs_val = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(flat[i])}, 1));
                    plus_args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){abs_val, expr_copy(norm_p)}, 2));
                }
                Expr* sum = NULL;
                if (N == 1) {
                    sum = plus_args[0];
                } else {
                    sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), plus_args, N));
                }
                free(plus_args);

                Expr* inv_p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(norm_p), expr_new_integer(-1)}, 2));
                result = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){sum, inv_p}, 2));
            }
            expr_free(norm_p);
        }

        if (flat) {
            for (int64_t i = 0; i < N; i++) expr_free(flat[i]);
            free(flat);
        }
        return result;
    }

    // Fallback for unhandled matrix norm (e.g. SVD max singular value)
    return NULL;
}
