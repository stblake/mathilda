/* Tensor contraction (Dot) and the shared dot2 helper.
 *
 * dot2 contracts the trailing axis of a with the leading axis of b and is
 * called from MatrixPower and the eigenvalue solver as well as from
 * builtin_dot itself.
 */

#include "linalg.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

/* Build a nested-List tensor of shape dims[0..rank-1] by consuming leaves
 * from `flat` in row-major order.  Each leaf is deep-copied. */
static Expr* build_tensor(int64_t* dims, int rank, Expr** flat, size_t* idx) {
    if (rank == 0) {
        return expr_copy(flat[(*idx)++]);
    }
    int64_t len = dims[0];
    Expr** args = NULL;
    if (len > 0) args = malloc(sizeof(Expr*) * len);
    for (int64_t i = 0; i < len; i++) {
        args[i] = build_tensor(dims + 1, rank - 1, flat, idx);
    }
    Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, len);
    if (args) free(args);
    return res;
}

Expr* dot2(Expr* a, Expr* b, bool* error_printed) {
    /* Fast path: both operands are dense Matrix ndarrays of rank <= 2 —
     * contract directly over the flat double buffers, no symbolic
     * Times/Plus per element. Higher-rank Matrix operands fall through to
     * the generic tensor path below via a nested-List conversion. */
    if (a->type == EXPR_MATRIX && b->type == EXPR_MATRIX) {
        bool shape_error = false;
        Expr* fast = matrix_dot2(a, b, &shape_error);
        if (fast) return fast;
        if (shape_error) {
            if (!*error_printed) {
                char* a_str = expr_to_string_fullform(a);
                char* b_str = expr_to_string_fullform(b);
                fprintf(stderr, "Dot::dotsh: Tensors %s and %s have incompatible shapes.\n", a_str, b_str);
                free(a_str);
                free(b_str);
                *error_printed = true;
            }
            return NULL;
        }
    }

    Expr* conv_a = NULL;
    Expr* conv_b = NULL;
    if (a->type == EXPR_MATRIX) { conv_a = matrix_to_nested_list(a); a = conv_a; }
    if (b->type == EXPR_MATRIX) { conv_b = matrix_to_nested_list(b); b = conv_b; }

    int64_t dimsA[64];
    int rankA = get_tensor_dims(a, dimsA);
    if (rankA <= 0) { // Not a tensor, or jagged
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t dimsB[64];
    int rankB = get_tensor_dims(b, dimsB);
    if (rankB <= 0) {
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t K = dimsA[rankA - 1];
    if (K != dimsB[0]) {
        if (!*error_printed) {
            char* a_str = expr_to_string_fullform(a);
            char* b_str = expr_to_string_fullform(b);
            fprintf(stderr, "Dot::dotsh: Tensors %s and %s have incompatible shapes.\n", a_str, b_str);
            free(a_str);
            free(b_str);
            *error_printed = true;
        }
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t N_A = 1; for(int i=0; i<rankA; i++) N_A *= dimsA[i];
    int64_t N_B = 1; for(int i=0; i<rankB; i++) N_B *= dimsB[i];

    Expr** flatA = NULL; if (N_A > 0) flatA = malloc(sizeof(Expr*) * N_A);
    Expr** flatB = NULL; if (N_B > 0) flatB = malloc(sizeof(Expr*) * N_B);
    size_t idxA = 0; if (N_A > 0) flatten_tensor(a, flatA, &idxA);
    size_t idxB = 0; if (N_B > 0) flatten_tensor(b, flatB, &idxB);

    int64_t R = K == 0 ? N_A : N_A / K;
    int64_t S = K == 0 ? N_B : N_B / K;

    Expr** flatC = NULL;
    if (R * S > 0) flatC = malloc(sizeof(Expr*) * (R * S));

    for (int64_t r = 0; r < R; r++) {
        for (int64_t s = 0; s < S; s++) {
            Expr** sum_args = NULL;
            if (K > 0) sum_args = malloc(sizeof(Expr*) * K);
            for (int64_t k = 0; k < K; k++) {
                Expr* a_elem = flatA[r * K + k];
                Expr* b_elem = flatB[k * S + s];
                Expr* t_args[2] = { expr_copy(a_elem), expr_copy(b_elem) };
                sum_args[k] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
            }
            if (K == 0) {
                flatC[r * S + s] = expr_new_integer(0);
            } else if (K == 1) {
                flatC[r * S + s] = sum_args[0];
            } else {
                flatC[r * S + s] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), sum_args, K));
            }
            if (sum_args) free(sum_args);
        }
    }

    if (flatA) { for(size_t i=0; i<idxA; i++) expr_free(flatA[i]); free(flatA); }
    if (flatB) { for(size_t i=0; i<idxB; i++) expr_free(flatB[i]); free(flatB); }

    int64_t dimsC[64];
    int rankC = rankA + rankB - 2;
    for (int i = 0; i < rankA - 1; i++) dimsC[i] = dimsA[i];
    for (int i = 0; i < rankB - 1; i++) dimsC[rankA - 1 + i] = dimsB[i + 1];

    size_t c_idx = 0;
    Expr* result = build_tensor(dimsC, rankC, flatC, &c_idx);
    if (flatC) {
        for (int64_t i = 0; i < R * S; i++) expr_free(flatC[i]);
        free(flatC);
    }

    if (conv_a) expr_free(conv_a);
    if (conv_b) expr_free(conv_b);
    return result;
}

Expr* builtin_dot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return NULL;
    if (count == 1) return expr_copy(res->data.function.args[0]);

    Expr** new_args = malloc(sizeof(Expr*) * count);
    for (size_t i=0; i<count; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    size_t new_count = count;

    bool changed = false;
    bool error_printed = false;
    for (size_t i = 0; i < new_count - 1; i++) {
        Expr* a = new_args[i];
        Expr* b = new_args[i+1];

        int64_t dA[64], dB[64];
        bool a_is_tensor = (a->type == EXPR_MATRIX) || get_tensor_dims(a, dA) > 0;
        bool b_is_tensor = (b->type == EXPR_MATRIX) || get_tensor_dims(b, dB) > 0;
        if (a_is_tensor && b_is_tensor) {
            Expr* d = dot2(a, b, &error_printed);
            if (d) {
                expr_free(a);
                expr_free(b);
                new_args[i] = d;
                for (size_t j = i + 2; j < new_count; j++) {
                    new_args[j - 1] = new_args[j];
                }
                new_count--;
                changed = true;
                i--; // re-check the new element at index i with i+1 (which shifted)
                if (i == (size_t)-1) i = 0; // In case we backed up past 0 (wait, i-- makes i=2^64-1, then i++ in loop makes it 0, correct)
            } else if (error_printed) {
                for (size_t j=0; j<new_count; j++) expr_free(new_args[j]);
                free(new_args);
                return NULL;
            }
        }
    }

    if (!changed) {
        for (size_t j=0; j<new_count; j++) expr_free(new_args[j]);
        free(new_args);
        return NULL;
    }

    Expr* final_res;
    if (new_count == 1) {
        final_res = new_args[0];
    } else {
        final_res = expr_new_function(expr_new_symbol(SYM_Dot), new_args, new_count);
    }

    free(new_args);
    return final_res;
}
