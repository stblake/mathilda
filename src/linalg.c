#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "expand.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static int get_tensor_dims(Expr* e, int64_t* dims) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL || strcmp(e->data.function.head->data.symbol, "List") != 0) {
        return 0; // rank 0
    }
    int64_t len = e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return 1;

    int sub_rank = get_tensor_dims(e->data.function.args[0], dims + 1);
    for (int64_t i = 1; i < len; i++) {
        int64_t cur_dims[64];
        int cur_rank = get_tensor_dims(e->data.function.args[i], cur_dims);
        if (cur_rank != sub_rank) return -1; // jagged
        for (int j = 0; j < sub_rank; j++) {
            if (cur_dims[j] != dims[j + 1]) return -1; // jagged
        }
    }
    return sub_rank + 1;
}

static void flatten_tensor(Expr* e, Expr** flat, size_t* idx) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL && strcmp(e->data.function.head->data.symbol, "List") == 0) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_tensor(e->data.function.args[i], flat, idx);
        }
    } else {
        flat[(*idx)++] = expr_copy(e);
    }
}

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
    Expr* res = expr_new_function(expr_new_symbol("List"), args, len);
    if (args) free(args);
    return res;
}

static Expr* dot2(Expr* a, Expr* b, bool* error_printed) {
    int64_t dimsA[64];
    int rankA = get_tensor_dims(a, dimsA);
    if (rankA <= 0) return NULL; // Not a tensor, or jagged

    int64_t dimsB[64];
    int rankB = get_tensor_dims(b, dimsB);
    if (rankB <= 0) return NULL;

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
                sum_args[k] = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 2));
            }
            if (K == 0) {
                flatC[r * S + s] = expr_new_integer(0);
            } else if (K == 1) {
                flatC[r * S + s] = sum_args[0];
            } else {
                flatC[r * S + s] = eval_and_free(expr_new_function(expr_new_symbol("Plus"), sum_args, K));
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
        if (get_tensor_dims(a, dA) > 0 && get_tensor_dims(b, dB) > 0) {
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
        final_res = expr_new_function(expr_new_symbol("Dot"), new_args, new_count);
    }
    
    free(new_args);
    return final_res;
}

static Expr* laplace_det(Expr** flat, int original_n, int n, int row, int* cols) {
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

Expr* builtin_cross(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t m = res->data.function.arg_count;
    if (m == 0) return NULL;

    size_t n = m + 1;
    bool valid = true;
    for (size_t i = 0; i < m; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || arg->data.function.head->type != EXPR_SYMBOL || strcmp(arg->data.function.head->data.symbol, "List") != 0) {
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
            result_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 2));
        } else {
            result_args[i] = det_val;
        }
    }

    Expr* final_res = expr_new_function(expr_new_symbol("List"), result_args, n);
    free(result_args);
    return final_res;
}

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
            return eval_and_free(expr_new_function(expr_new_symbol("Abs"), args, 1));
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
        
        if (p && p->type == EXPR_SYMBOL && strcmp(p->data.symbol, "Infinity") == 0) {
            if (N == 0) {
                result = expr_new_integer(0);
            } else {
                Expr** max_args = malloc(sizeof(Expr*) * N);
                for (int64_t i = 0; i < N; i++) {
                    max_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Abs"), (Expr*[]){expr_copy(flat[i])}, 1));
                }
                result = eval_and_free(expr_new_function(expr_new_symbol("Max"), max_args, N));
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
                    Expr* abs_val = eval_and_free(expr_new_function(expr_new_symbol("Abs"), (Expr*[]){expr_copy(flat[i])}, 1));
                    plus_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){abs_val, expr_copy(norm_p)}, 2));
                }
                Expr* sum = NULL;
                if (N == 1) {
                    sum = plus_args[0];
                } else {
                    sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"), plus_args, N));
                }
                free(plus_args);
                
                Expr* inv_p = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(norm_p), expr_new_integer(-1)}, 2));
                result = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){sum, inv_p}, 2));
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

static int64_t get_default_trace_depth(Expr* list) {
    int64_t depth = 0;
    Expr* curr = list;
    while (curr->type == EXPR_FUNCTION && curr->data.function.head->type == EXPR_SYMBOL && strcmp(curr->data.function.head->data.symbol, "List") == 0) {
        depth++;
        if (curr->data.function.arg_count == 0) break;
        curr = curr->data.function.args[0];
    }
    return depth;
}

static Expr* extract_diagonal_element(Expr* list, int64_t n, size_t index) {
    Expr* curr = list;
    for (int64_t level = 0; level < n; level++) {
        if (curr->type != EXPR_FUNCTION || curr->data.function.head->type != EXPR_SYMBOL || strcmp(curr->data.function.head->data.symbol, "List") != 0) {
            return NULL; // Not a list at this level
        }
        if (index >= curr->data.function.arg_count) {
            return NULL; // Index out of bounds
        }
        curr = curr->data.function.args[index];
    }
    return expr_copy(curr);
}

Expr* builtin_tr(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count < 1 || count > 3) return NULL;

    Expr* list = res->data.function.args[0];

    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || strcmp(list->data.function.head->data.symbol, "List") != 0) {
        return expr_copy(list);
    }

    Expr* f = NULL;
    bool free_f = false;
    if (count >= 2) {
        f = res->data.function.args[1];
    } else {
        f = expr_new_symbol("Plus");
        free_f = true;
    }

    int64_t n = 0;
    if (count == 3) {
        Expr* n_expr = res->data.function.args[2];
        if (n_expr->type == EXPR_INTEGER) {
            n = n_expr->data.integer;
            if (n < 0) n = 0; 
        } else {
            if (free_f) expr_free(f);
            return NULL;
        }
    } else {
        n = get_default_trace_depth(list);
        if (n == 0) n = 1;
    }

    if (n == 0) {
        if (free_f) expr_free(f);
        return expr_copy(list);
    }

    size_t cap = 16;
    size_t num_elems = 0;
    Expr** elements = malloc(sizeof(Expr*) * cap);

    size_t i = 0;
    while (true) {
        Expr* elem = extract_diagonal_element(list, n, i);
        if (!elem) {
            break;
        }
        if (num_elems == cap) {
            cap *= 2;
            elements = realloc(elements, sizeof(Expr*) * cap);
        }
        elements[num_elems++] = elem;
        i++;
    }

    Expr* result;
    if (num_elems == 0) {
        result = eval_and_free(expr_new_function(expr_copy(f), NULL, 0));
    } else {
        result = eval_and_free(expr_new_function(expr_copy(f), elements, num_elems));
    }
    
    free(elements);
    if (free_f) expr_free(f);

    return result;
}

static Expr* exact_div_wrapper(Expr* num, Expr* den) {
    if (is_zero_poly(num)) return expr_new_integer(0);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return expr_expand(num);

    Expr* exp_num = expr_expand(num);
    Expr* exp_den = expr_expand(den);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(exp_num, &vars, &v_count, &v_cap);
    collect_variables(exp_den, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    Expr* res = exact_poly_div(exp_num, exp_den, vars, v_count);

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);

    if (res) {
        expr_free(exp_num);
        expr_free(exp_den);
        return res;
    } else {
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){exp_num, t}, 2));
        return r;
    }
}

Expr* builtin_rowreduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        return expr_copy(arg);
    }
    
    int m = (int)dims[0];
    int n = (int)dims[1];
    
    Expr** matrix = malloc(sizeof(Expr*) * m * n);
    size_t idx = 0;
    flatten_tensor(arg, matrix, &idx);
    
    Expr* P = expr_new_integer(1);
    int r = 0;
    
    for (int c = 0; c < n && r < m; c++) {
        int pivot_row = -1;
        for (int i = r; i < m; i++) {
            if (!is_zero_poly(matrix[i * n + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;
        
        if (pivot_row != r) {
            for (int j = 0; j < n; j++) {
                Expr* tmp = matrix[r * n + j];
                matrix[r * n + j] = matrix[pivot_row * n + j];
                matrix[pivot_row * n + j] = tmp;
            }
        }
        
        Expr* pivot = matrix[r * n + c];
        
        for (int i = 0; i < m; i++) {
            if (i == r) continue;
            Expr* M_ic = matrix[i * n + c];
            if (is_zero_poly(M_ic)) {
                for (int j = 0; j < n; j++) {
                    if (j == c) continue;
                    Expr* num_eval = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * n + j])}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = new_val;
                }
            } else {
                for (int j = 0; j < n; j++) {
                    if (j == c) continue;
                    Expr* t1 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * n + j])}, 2));
                    Expr* t2 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(M_ic), expr_copy(matrix[r * n + j])}, 2));
                    Expr* t2_neg = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), t2}, 2));
                    Expr* num_eval = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){t1, t2_neg}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = new_val;
                }
            }
            expr_free(matrix[i * n + c]);
            matrix[i * n + c] = expr_new_integer(0);
        }
        
        expr_free(P);
        P = expr_copy(pivot);
        r++;
    }
    expr_free(P);
    
    for (int i = 0; i < m; i++) {
        Expr* leading = NULL;
        int lead_j = -1;
        for (int j = 0; j < n; j++) {
            if (!is_zero_poly(matrix[i * n + j])) {
                leading = expr_copy(matrix[i * n + j]);
                lead_j = j;
                break;
            }
        }
        if (leading) {
            for (int j = lead_j; j < n; j++) {
                if (j == lead_j) {
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = expr_new_integer(1);
                } else if (!is_zero_poly(matrix[i * n + j])) {
                    Expr* num_val = matrix[i * n + j];
                    Expr* den_val = expr_copy(leading);
                    
                    Expr* g = eval_and_free(expr_new_function(expr_new_symbol("PolynomialGCD"), (Expr*[]){expr_copy(num_val), expr_copy(den_val)}, 2));
                    Expr* new_num = exact_div_wrapper(num_val, g);
                    Expr* new_den = exact_div_wrapper(den_val, g);
                    expr_free(g);
                    
                    if (new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
                        Expr* t = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, expr_new_integer(-1)}, 2));
                        new_num = expr_expand(t);
                        expr_free(t);
                        /* Replace, don't mutate: the integer atom may be
                         * shared (M3 atom-sharing). */
                        int64_t v = -new_den->data.integer;
                        expr_free(new_den);
                        new_den = expr_new_integer(v);
                    }
                    
                    if (new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
                        expr_free(new_den);
                        matrix[i * n + j] = new_num;
                    } else {
                        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){new_den, expr_new_integer(-1)}, 2));
                        Expr* final_val = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, inv_den}, 2));
                        matrix[i * n + j] = expr_expand(final_val);
                        expr_free(final_val);
                    }
                    expr_free(num_val);
                    expr_free(den_val);
                }
            }
            expr_free(leading);
        }
    }
    
    Expr** rows = malloc(sizeof(Expr*) * m);
    for (int i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * n + j];
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, n);
        free(row_elems);
    }
    
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, m);
    free(rows);
    free(matrix);
    
    return result;
}

Expr* builtin_identitymatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t m = -1, n = -1;

    if (arg->type == EXPR_INTEGER) {
        m = arg->data.integer;
        n = m;
    } else if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && strcmp(arg->data.function.head->data.symbol, "List") == 0) {
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
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol("List"), rows, m);
    free(rows);
    return result;
}

Expr* builtin_diagonalmatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || strcmp(list->data.function.head->data.symbol, "List") != 0) {
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
        } else if (dim_expr->type == EXPR_FUNCTION && dim_expr->data.function.head->type == EXPR_SYMBOL && strcmp(dim_expr->data.function.head->data.symbol, "List") == 0) {
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
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol("List"), rows, m);
    free(rows);
    return result;
}

Expr* builtin_inverse(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);

    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr, "Inverse::matsq: Argument %s at position 1 is not a non-empty square matrix.\n", arg_str);
        free(arg_str);
        return NULL;
    }

    int n = (int)dims[0];
    int cols = 2 * n;

    /* Build augmented matrix [A | I] */
    Expr** matrix = malloc(sizeof(Expr*) * n * cols);
    size_t idx = 0;
    Expr** flat_a = malloc(sizeof(Expr*) * n * n);
    flatten_tensor(arg, flat_a, &idx);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i * cols + j] = flat_a[i * n + j]; /* take ownership */
        }
        for (int j = 0; j < n; j++) {
            matrix[i * cols + n + j] = expr_new_integer(i == j ? 1 : 0);
        }
    }
    free(flat_a); /* elements transferred to matrix */

    /* Fraction-free Gauss-Jordan elimination (Bareiss-like, same as RowReduce) */
    Expr* P = expr_new_integer(1);
    int r = 0;

    for (int c = 0; c < n && r < n; c++) {
        /* Find pivot */
        int pivot_row = -1;
        for (int i = r; i < n; i++) {
            if (!is_zero_poly(matrix[i * cols + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) {
            /* Singular matrix */
            char* arg_str = expr_to_string(arg);
            fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
            free(arg_str);
            expr_free(P);
            for (int i = 0; i < n * cols; i++) expr_free(matrix[i]);
            free(matrix);
            return NULL;
        }

        /* Swap rows if needed */
        if (pivot_row != r) {
            for (int j = 0; j < cols; j++) {
                Expr* tmp = matrix[r * cols + j];
                matrix[r * cols + j] = matrix[pivot_row * cols + j];
                matrix[pivot_row * cols + j] = tmp;
            }
        }

        Expr* pivot = matrix[r * cols + c];

        /* Eliminate column c in all other rows */
        for (int i = 0; i < n; i++) {
            if (i == r) continue;
            Expr* M_ic = matrix[i * cols + c];
            if (is_zero_poly(M_ic)) {
                for (int j = 0; j < cols; j++) {
                    if (j == c) continue;
                    Expr* num_eval = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * cols + j])}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * cols + j]);
                    matrix[i * cols + j] = new_val;
                }
            } else {
                for (int j = 0; j < cols; j++) {
                    if (j == c) continue;
                    Expr* t1 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * cols + j])}, 2));
                    Expr* t2 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(M_ic), expr_copy(matrix[r * cols + j])}, 2));
                    Expr* t2_neg = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_new_integer(-1), t2}, 2));
                    Expr* num_eval = eval_and_free(expr_new_function(
                        expr_new_symbol("Plus"),
                        (Expr*[]){t1, t2_neg}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * cols + j]);
                    matrix[i * cols + j] = new_val;
                }
            }
            expr_free(matrix[i * cols + c]);
            matrix[i * cols + c] = expr_new_integer(0);
        }

        expr_free(P);
        P = expr_copy(pivot);
        r++;
    }
    expr_free(P);

    /* If we didn't get n pivots, matrix is singular */
    if (r < n) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
        free(arg_str);
        for (int i = 0; i < n * cols; i++) expr_free(matrix[i]);
        free(matrix);
        return NULL;
    }

    /* Normalize: divide each row of the right half by its diagonal element */
    for (int i = 0; i < n; i++) {
        Expr* diag = matrix[i * cols + i];
        for (int j = n; j < cols; j++) {
            if (is_zero_poly(matrix[i * cols + j])) continue;

            Expr* num_val = matrix[i * cols + j];
            Expr* den_val = expr_copy(diag);

            /* Cancel common factors via PolynomialGCD */
            Expr* g = eval_and_free(expr_new_function(
                expr_new_symbol("PolynomialGCD"),
                (Expr*[]){expr_copy(num_val), expr_copy(den_val)}, 2));
            Expr* new_num = exact_div_wrapper(num_val, g);
            Expr* new_den = exact_div_wrapper(den_val, g);
            expr_free(g);

            /* Normalize sign so denominator is positive */
            if (new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
                Expr* t = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){new_num, expr_new_integer(-1)}, 2));
                new_num = expr_expand(t);
                expr_free(t);
                /* Replace, don't mutate: the integer atom may be shared. */
                int64_t v = -new_den->data.integer;
                expr_free(new_den);
                new_den = expr_new_integer(v);
            }

            if (new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
                expr_free(new_den);
                matrix[i * cols + j] = new_num;
            } else {
                Expr* inv_den = eval_and_free(expr_new_function(
                    expr_new_symbol("Power"),
                    (Expr*[]){new_den, expr_new_integer(-1)}, 2));
                Expr* final_val = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){new_num, inv_den}, 2));
                matrix[i * cols + j] = expr_expand(final_val);
                expr_free(final_val);
            }
            expr_free(num_val);
            expr_free(den_val);
        }
    }

    /* Extract the right half [I | A^{-1}] -> A^{-1} */
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * cols + n + j]; /* take ownership */
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, n);
        free(row_elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);

    /* Free the left half and the flat matrix array */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            expr_free(matrix[i * cols + j]);
        }
        /* right half elements were transferred to result */
    }
    free(matrix);

    return result;
}

/* Matrix multiplication helper: computes A.B for two square n x n matrices
 * represented as List[List[...], ...] expressions. Returns a new Expr*. */
static Expr* mat_dot(Expr* a, Expr* b) {
    bool err = false;
    return dot2(a, b, &err);
}

Expr* builtin_matrixpower(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* m = res->data.function.args[0];
    Expr* exp_arg = res->data.function.args[1];

    /* Validate matrix: must be rank 2, square, non-empty */
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        if (rank >= 0) {
            char* m_str = expr_to_string(m);
            fprintf(stderr, "MatrixPower::matsq: Argument %s at position 1 is not a non-empty square matrix.\n", m_str);
            free(m_str);
        }
        return NULL;
    }
    int n = (int)dims[0];

    /* Validate exponent: must be an integer */
    bool is_int = (exp_arg->type == EXPR_INTEGER);
    bool is_bigint = (exp_arg->type == EXPR_BIGINT);
    bool is_rational = false;
    bool is_real = (exp_arg->type == EXPR_REAL);

    if (exp_arg->type == EXPR_FUNCTION && exp_arg->data.function.head->type == EXPR_SYMBOL
        && strcmp(exp_arg->data.function.head->data.symbol, "Rational") == 0) {
        is_rational = true;
    }

    /* Fractional powers: warn and return unevaluated */
    if (is_rational || is_real) {
        char* exp_str = expr_to_string(exp_arg);
        fprintf(stderr, "MatrixPower::fract: Fractional matrix powers are not currently supported. Exponent %s is not an integer.\n", exp_str);
        free(exp_str);
        return NULL;
    }

    /* Symbolic or non-numeric exponent: return unevaluated */
    if (!is_int && !is_bigint) return NULL;

    /* Get the exponent value; for bigint, check it fits in int64_t */
    int64_t exp_val = 0;
    if (is_int) {
        exp_val = exp_arg->data.integer;
    } else {
        /* EXPR_BIGINT: check if it fits in int64_t range */
        if (mpz_fits_slong_p(exp_arg->data.bigint)) {
            exp_val = mpz_get_si(exp_arg->data.bigint);
        } else {
            /* Exponent too large to compute */
            return NULL;
        }
    }

    /* If argc == 3, validate vector argument */
    Expr* vec = NULL;
    if (argc == 3) {
        vec = res->data.function.args[2];
        int64_t vdims[64];
        int vrank = get_tensor_dims(vec, vdims);
        if (vrank != 1 || vdims[0] != n) {
            char* v_str = expr_to_string(vec);
            fprintf(stderr, "MatrixPower::vecsh: Vector %s has incompatible length for matrix of size %d.\n", v_str, n);
            free(v_str);
            return NULL;
        }
    }

    /* For negative exponents, compute inverse first */
    Expr* base = NULL;
    int64_t abs_exp = exp_val;
    if (exp_val < 0) {
        abs_exp = -exp_val;
        /* Compute Inverse[m] */
        Expr* inv_call = expr_new_function(expr_new_symbol("Inverse"),
            (Expr*[]){expr_copy(m)}, 1);
        Expr* inv_result = evaluate(inv_call);
        expr_free(inv_call);

        /* Check if Inverse returned unevaluated (singular matrix) */
        if (inv_result->type == EXPR_FUNCTION && inv_result->data.function.head->type == EXPR_SYMBOL
            && strcmp(inv_result->data.function.head->data.symbol, "Inverse") == 0) {
            expr_free(inv_result);
            return NULL; /* Singular: Inverse already printed warning */
        }
        base = inv_result;
    } else {
        base = expr_copy(m);
    }

    /* Compute matrix power via binary exponentiation */
    Expr* result = NULL;

    if (abs_exp == 0) {
        /* M^0 = IdentityMatrix[n] */
        expr_free(base);
        Expr* id_args[1];
        id_args[0] = expr_new_integer(n);
        Expr* id_call = expr_new_function(expr_new_symbol("IdentityMatrix"), id_args, 1);
        result = evaluate(id_call);
        expr_free(id_call);
    } else {
        /* Binary exponentiation: square-and-multiply */
        result = NULL;
        Expr* sq = base; /* current square, takes ownership of base */

        int64_t e = abs_exp;
        while (e > 0) {
            if (e & 1) {
                if (result == NULL) {
                    result = expr_copy(sq);
                } else {
                    Expr* new_result = mat_dot(result, sq);
                    if (!new_result) {
                        /* Should not happen for valid square matrices */
                        expr_free(result);
                        expr_free(sq);
                        return NULL;
                    }
                    Expr* evaluated = evaluate(new_result);
                    expr_free(new_result);
                    expr_free(result);
                    result = evaluated;
                }
            }
            e >>= 1;
            if (e > 0) {
                Expr* new_sq = mat_dot(sq, sq);
                if (!new_sq) {
                    if (result) expr_free(result);
                    expr_free(sq);
                    return NULL;
                }
                Expr* evaluated = evaluate(new_sq);
                expr_free(new_sq);
                expr_free(sq);
                sq = evaluated;
            }
        }
        expr_free(sq);
    }

    /* If vector argument provided, compute result . v */
    if (vec && result) {
        Expr* dotted = mat_dot(result, vec);
        if (dotted) {
            Expr* evaluated = evaluate(dotted);
            expr_free(dotted);
            expr_free(result);
            result = evaluated;
        }
        /* If dot fails, just return matrix power without applying to vector */
    }

    /* Note: evaluator frees res after we return non-NULL */
    return result;
}

void linalg_init(void) {
    symtab_add_builtin("Dot", builtin_dot);
    symtab_get_def("Dot")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_add_builtin("Det", builtin_det);
    symtab_get_def("Det")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Cross", builtin_cross);
    symtab_get_def("Cross")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Norm", builtin_norm);
    symtab_get_def("Norm")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Tr", builtin_tr);
    symtab_get_def("Tr")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("RowReduce", builtin_rowreduce);
    symtab_get_def("RowReduce")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("IdentityMatrix", builtin_identitymatrix);
    symtab_get_def("IdentityMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DiagonalMatrix", builtin_diagonalmatrix);
    symtab_get_def("DiagonalMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Inverse", builtin_inverse);
    symtab_get_def("Inverse")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatrixPower", builtin_matrixpower);
    symtab_get_def("MatrixPower")->attributes |= ATTR_PROTECTED;
}
