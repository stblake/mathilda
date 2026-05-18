#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

int get_tensor_dims(Expr* e, int64_t* dims) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL || e->data.function.head->data.symbol != SYM_List) {
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

void flatten_tensor(Expr* e, Expr** flat, size_t* idx) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol == SYM_List) {
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
        
        if (p && p->type == EXPR_SYMBOL && p->data.symbol == SYM_Infinity) {
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
    while (curr->type == EXPR_FUNCTION && curr->data.function.head->type == EXPR_SYMBOL && curr->data.function.head->data.symbol == SYM_List) {
        depth++;
        if (curr->data.function.arg_count == 0) break;
        curr = curr->data.function.args[0];
    }
    return depth;
}

static Expr* extract_diagonal_element(Expr* list, int64_t n, size_t index) {
    Expr* curr = list;
    for (int64_t level = 0; level < n; level++) {
        if (curr->type != EXPR_FUNCTION || curr->data.function.head->type != EXPR_SYMBOL || curr->data.function.head->data.symbol != SYM_List) {
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

    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || list->data.function.head->data.symbol != SYM_List) {
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

Expr* exact_div_wrapper(Expr* num, Expr* den) {
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
        && exp_arg->data.function.head->data.symbol == SYM_Rational) {
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
            && inv_result->data.function.head->data.symbol == SYM_Inverse) {
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

/* ============================================================ *
 *  Eigenvalues and Eigenvectors                                  *
 *                                                                 *
 *  Strategy: compute the characteristic polynomial p(λ) =         *
 *  Det[m - λ*I] (or Det[m - λ*a] for generalised eigenvalues),    *
 *  then route to the public `Solve` builtin so its rationalise -> *
 *  solve -> numericalize pipeline handles approximate matrices    *
 *  automatically.                                                 *
 *                                                                 *
 *  The internal-only lambda symbol is context-qualified           *
 *  ("Eigenvalues`Lambda") so it does not collide with user names. *
 * ============================================================ */

/* Canonical interned name of the lambda symbol used as the char-poly
 * variable.  Pointer is stable across all calls (interner is idempotent). */
static const char* eigen_lambda_name(void) {
    static const char* s = NULL;
    if (!s) s = intern_symbol("Eigenvalues`Lambda");
    return s;
}

/* True iff `m` is a non-empty n×n list-of-lists matrix.  Writes n to
 * *n_out on success. */
static bool eigen_is_square_matrix(Expr* m, int64_t* n_out) {
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) return false;
    *n_out = dims[0];
    return true;
}

/* Classify `arg`:
 *   {m, a} with both m and a square n×n  -> generalised eigenvalue case
 *   single square n×n                    -> ordinary eigenvalue case
 * On success, *m_out, *a_out (NULL for ordinary) and *n_out are set.
 * Returns false when `arg` is not a recognisable matrix shape. */
static bool eigen_extract_matrix_pair(Expr* arg, Expr** m_out, Expr** a_out,
                                       int64_t* n_out) {
    if (arg->type == EXPR_FUNCTION
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol == SYM_List
        && arg->data.function.arg_count == 2) {
        Expr* m = arg->data.function.args[0];
        Expr* a = arg->data.function.args[1];
        int64_t nm, na;
        if (eigen_is_square_matrix(m, &nm)
            && eigen_is_square_matrix(a, &na) && nm == na) {
            *m_out = m; *a_out = a; *n_out = nm;
            return true;
        }
    }
    int64_t n;
    if (eigen_is_square_matrix(arg, &n)) {
        *m_out = arg; *a_out = NULL; *n_out = n;
        return true;
    }
    return false;
}

/* Build the n×n matrix m - λ*a (or m - λ*I when a_or_null == NULL) where
 * each entry is an evaluated polynomial-in-lambda expression.  Caller
 * owns the returned List-of-Lists. */
static Expr* eigen_build_lambda_matrix(Expr* m, Expr* a_or_null,
                                        const char* lambda_name, int64_t n) {
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int64_t i = 0; i < n; i++) {
        Expr* m_row = m->data.function.args[i];
        Expr* a_row = a_or_null ? a_or_null->data.function.args[i] : NULL;
        Expr** elems = malloc(sizeof(Expr*) * n);
        for (int64_t j = 0; j < n; j++) {
            Expr* mij = m_row->data.function.args[j];
            Expr* sub;
            if (a_or_null) {
                Expr* aij = a_row->data.function.args[j];
                Expr* neg_lam_a = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1),
                               expr_new_symbol(lambda_name),
                               expr_copy(aij) }, 3));
                sub = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij), neg_lam_a }, 2));
            } else if (i == j) {
                Expr* neg_lam = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1),
                               expr_new_symbol(lambda_name) }, 2));
                sub = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij), neg_lam }, 2));
            } else {
                sub = expr_copy(mij);
            }
            elems[j] = sub;
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), elems, n);
        free(elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return result;
}

/* Compute Det of the matrix `matrix` (n×n) using existing Laplace expansion.
 * Returns a freshly allocated Expr*. */
static Expr* eigen_compute_det(Expr* matrix, int n) {
    Expr** flat = malloc(sizeof(Expr*) * n * n);
    size_t idx = 0;
    flatten_tensor(matrix, flat, &idx);
    int* cols = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) cols[i] = i;
    Expr* det_val = laplace_det(flat, n, n, 0, cols);
    free(cols);
    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);
    return det_val;
}

/* ---- Faddeev-Leverrier char-poly fast path for the ordinary case ----
 *
 * The Laplace-expansion det used by eigen_compute_det is O(n!) -- usable
 * up to about n = 8 but quickly intolerable beyond.  When the user just
 * wants the eigenvalues of an n×n matrix (no generalised `a`), we can run
 * Faddeev-Leverrier instead: it computes the coefficients of
 * det(λI - A) directly in O(n^4) matrix operations on the constant
 * matrix A, never touching a polynomial-in-λ entry.  The resulting
 * polynomial has the same roots as det(A - λI) (they differ only by the
 * (-1)^n sign, which Solve ignores).
 */

/* Trace of an n×n matrix expressed as List of Lists. */
static Expr* eigen_mat_trace(Expr* M, int n) {
    Expr** terms = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        terms[i] = expr_copy(
            M->data.function.args[i]->data.function.args[i]);
    }
    Expr* res = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), terms, n));
    free(terms);
    return res;
}

/* Return M - s*I (entrywise) as a freshly allocated matrix. */
static Expr* eigen_mat_minus_scalar_id(Expr* M, Expr* s, int n) {
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        Expr** row = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
            Expr* mij = M->data.function.args[i]->data.function.args[j];
            if (i == j) {
                row[j] = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij),
                               eval_and_free(expr_new_function(
                                   expr_new_symbol("Times"),
                                   (Expr*[]){ expr_new_integer(-1),
                                              expr_copy(s) }, 2)) }, 2));
            } else {
                row[j] = expr_copy(mij);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row, n);
        free(row);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return result;
}

/* Matrix multiply via mat_dot + evaluate (so the entries are normalised). */
static Expr* eigen_mat_mul(Expr* A, Expr* B) {
    bool err = false;
    Expr* prod = dot2(A, B, &err);
    if (!prod) return NULL;
    return eval_and_free(prod);
}

/* Faddeev-Leverrier-Souriau characteristic polynomial.
 *
 *     M_1 = A,    p_1 = Tr(M_1)
 *     for k = 2, ..., n:
 *         M_k = A . (M_{k-1} − p_{k-1} I)
 *         p_k = Tr(M_k) / k
 *
 * Yields p(λ) = det(λI − A) = λ^n − p_1 λ^{n-1} − p_2 λ^{n-2} − … − p_n.
 * Every coefficient other than the leading one is −p_k (no alternating sign).
 *
 * Returns the polynomial in the lambda variable.  Caller owns the result. */
static Expr* eigen_char_poly_faddeev(Expr* A, const char* lambda_name, int n) {
    Expr* M = expr_copy(A);                              /* M_1 = A */
    Expr* p_prev = eigen_mat_trace(M, n);                /* p_1 */

    Expr** coeffs = malloc(sizeof(Expr*) * (n + 1));
    coeffs[n] = expr_new_integer(1);
    coeffs[n - 1] = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), expr_copy(p_prev) }, 2));

    for (int k = 2; k <= n; k++) {
        Expr* shifted = eigen_mat_minus_scalar_id(M, p_prev, n);
        expr_free(M);
        M = eigen_mat_mul(A, shifted);
        expr_free(shifted);
        if (!M) {
            expr_free(p_prev);
            expr_free(coeffs[n]);
            expr_free(coeffs[n - 1]);
            free(coeffs);
            return NULL;
        }
        Expr* tr = eigen_mat_trace(M, n);
        Expr* p_k = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ tr,
                       eval_and_free(expr_new_function(
                           expr_new_symbol("Power"),
                           (Expr*[]){ expr_new_integer(k),
                                      expr_new_integer(-1) }, 2)) }, 2));
        expr_free(p_prev);
        p_prev = p_k;

        /* coefficient of λ^{n-k} in det(λI − A) is −p_k. */
        coeffs[n - k] = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(p_k) }, 2));
    }

    expr_free(p_prev);
    expr_free(M);

    /* Build polynomial in lambda_name. */
    Expr** terms = malloc(sizeof(Expr*) * (n + 1));
    size_t tcount = 0;
    for (int k = 0; k <= n; k++) {
        if (k == 0) {
            terms[tcount++] = coeffs[0];
            continue;
        }
        Expr* lam_pow;
        if (k == 1) {
            lam_pow = expr_new_symbol(lambda_name);
        } else {
            lam_pow = eval_and_free(expr_new_function(
                expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol(lambda_name),
                           expr_new_integer(k) }, 2));
        }
        terms[tcount++] = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ coeffs[k], lam_pow }, 2));
    }
    free(coeffs);
    Expr* poly = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), terms, tcount));
    free(terms);
    return poly;
}

/* Solve poly == 0 for the lambda variable via the public Solve builtin so
 * the inexact-preprocessing pipeline (rationalise -> solve -> numericalize)
 * runs automatically.  Cubics/Quartics options are forwarded so the user
 * can request held Root[] objects via Cubics -> False / Quartics -> False.
 * Returns the solution List, or NULL on failure. */
static Expr* eigen_solve_poly(Expr* poly, const char* lambda_name,
                              bool cubics_radical, bool quartics_radical) {
    Expr* eq = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ expr_copy(poly), expr_new_integer(0) }, 2);
    Expr* lam = expr_new_symbol(lambda_name);
    Expr* opt_cubics = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol("Cubics"),
                   expr_new_symbol(cubics_radical ? "True" : "False") }, 2);
    Expr* opt_quartics = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol("Quartics"),
                   expr_new_symbol(quartics_radical ? "True" : "False") }, 2);
    Expr* solve_call = expr_new_function(expr_new_symbol("Solve"),
        (Expr*[]){ eq, lam, opt_cubics, opt_quartics }, 4);
    return eval_and_free(solve_call);
}

/* Extract eigenvalues from Solve's output:
 *   {{λ -> v1}, {λ -> v2}, ...}  ->  freshly-owned array [v1, v2, ...]
 * Empty solution `{{}}` (tautology -- e.g. the input was the zero matrix
 * and m == a generalised case) yields a NULL value placeholder so callers
 * can pad with Indeterminate. */
static Expr** eigen_extract_values(Expr* solutions, size_t* count_out) {
    *count_out = 0;
    if (!solutions || solutions->type != EXPR_FUNCTION
        || solutions->data.function.head->type != EXPR_SYMBOL
        || solutions->data.function.head->data.symbol != SYM_List) {
        return NULL;
    }
    size_t n = solutions->data.function.arg_count;
    if (n == 0) return NULL;
    Expr** out = malloc(sizeof(Expr*) * n);
    size_t out_count = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = solutions->data.function.args[i];
        if (sol->type != EXPR_FUNCTION
            || sol->data.function.head->type != EXPR_SYMBOL
            || sol->data.function.head->data.symbol != SYM_List
            || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (rule->type != EXPR_FUNCTION
            || rule->data.function.head->type != EXPR_SYMBOL
            || rule->data.function.head->data.symbol != SYM_Rule
            || rule->data.function.arg_count != 2) continue;
        out[out_count++] = expr_copy(rule->data.function.args[1]);
    }
    *count_out = out_count;
    return out;
}

/* Chop very small inexact imaginary / real parts.  Threshold is relative
 * to the magnitude of the value: anything below 1e-10 * |val| (or below
 * 1e-12 absolute when |val| is near zero) is dropped.  Used to clean up
 * numerical noise introduced by Cardano-style closed-form root formulas
 * applied to real polynomials. */
static double eigen_part_to_double(Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_BIGINT) return mpz_get_d(e->data.bigint);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        double p = eigen_part_to_double(e->data.function.args[0]);
        double q = eigen_part_to_double(e->data.function.args[1]);
        return q == 0 ? 0 : p / q;
    }
    return NAN;
}

static Expr* eigen_chop(Expr* val) {
    if (!val) return NULL;
    /* Complex[re, im] - drop im (or re) if small. */
    if (val->type == EXPR_FUNCTION
        && val->data.function.head->type == EXPR_SYMBOL
        && val->data.function.head->data.symbol == SYM_Complex
        && val->data.function.arg_count == 2) {
        Expr* re = val->data.function.args[0];
        Expr* im = val->data.function.args[1];
        double rd = eigen_part_to_double(re);
        double id = eigen_part_to_double(im);
        if (!isnan(rd) && !isnan(id)) {
            double mag = fabs(rd) + fabs(id);
            double thresh = 1e-10 * mag + 1e-12;
            bool drop_im = fabs(id) < thresh;
            bool drop_re = fabs(rd) < thresh;
            if (drop_im && drop_re) return expr_new_real(0.0);
            if (drop_im) return expr_copy(re);
            if (drop_re) {
                /* Pure imaginary: i * im */
                return eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(im), expr_new_symbol("I") }, 2));
            }
        }
    }
    return expr_copy(val);
}

/* Reduce val to a concrete-real double via N[Abs[val]].  Returns false
 * when the reduction does not collapse to a single Real / Integer /
 * Rational / MPFR / BigInt -- i.e. the value carries symbolic terms. */
static bool eigen_abs_to_double(Expr* val, double* out) {
    Expr* abs_e = eval_and_free(expr_new_function(
        expr_new_symbol("Abs"), (Expr*[]){ expr_copy(val) }, 1));
    Expr* n_abs = eval_and_free(expr_new_function(
        expr_new_symbol("N"), (Expr*[]){ abs_e }, 1));
    bool ok = true;
    double d = 0;
    if (n_abs->type == EXPR_REAL) d = n_abs->data.real;
    else if (n_abs->type == EXPR_INTEGER) d = (double)n_abs->data.integer;
    else if (n_abs->type == EXPR_BIGINT) d = mpz_get_d(n_abs->data.bigint);
#ifdef USE_MPFR
    else if (n_abs->type == EXPR_MPFR)
        d = mpfr_get_d(n_abs->data.mpfr, MPFR_RNDN);
#endif
    else ok = false;
    expr_free(n_abs);
    *out = d;
    return ok;
}

typedef struct {
    Expr*  val;       /* borrowed during sort */
    double abs_d;
    size_t orig_idx;  /* tiebreaker to keep stable order */
} EigenSortKey;

static int eigen_sort_cmp_desc(const void* a, const void* b) {
    const EigenSortKey* ka = (const EigenSortKey*)a;
    const EigenSortKey* kb = (const EigenSortKey*)b;
    if (ka->abs_d > kb->abs_d) return -1;
    if (ka->abs_d < kb->abs_d) return 1;
    if (ka->orig_idx < kb->orig_idx) return -1;
    if (ka->orig_idx > kb->orig_idx) return 1;
    return 0;
}

/* Sort vals[] by descending |λ| if every entry reduces to a concrete real
 * via N[Abs[...]] (and is not Infinity/Indeterminate).  Otherwise leave
 * vals[] in Solve's natural order. */
static void eigen_sort_by_abs_desc(Expr** vals, size_t n) {
    if (n <= 1) return;
    EigenSortKey* keys = malloc(sizeof(EigenSortKey) * n);
    bool all_numeric = true;
    for (size_t i = 0; i < n; i++) {
        keys[i].val = vals[i];
        keys[i].orig_idx = i;
        /* Infinity should sort first regardless. */
        if (vals[i]->type == EXPR_SYMBOL
            && vals[i]->data.symbol == SYM_Infinity) {
            keys[i].abs_d = INFINITY;
            continue;
        }
        if (vals[i]->type == EXPR_SYMBOL
            && vals[i]->data.symbol == SYM_Indeterminate) {
            keys[i].abs_d = -INFINITY; /* sort last */
            continue;
        }
        if (!eigen_abs_to_double(vals[i], &keys[i].abs_d)) {
            all_numeric = false;
            break;
        }
    }
    if (all_numeric) {
        qsort(keys, n, sizeof(EigenSortKey), eigen_sort_cmp_desc);
        for (size_t i = 0; i < n; i++) vals[i] = keys[i].val;
    }
    free(keys);
}

/* Compute null-space basis of n×n matrix M.  Uses RowReduce then constructs
 * a basis vector for each free column.  Returns an array of `*count_out`
 * Expr* vectors (each a List of n elements). */
static Expr** eigen_null_space(Expr* M, int n, size_t* count_out) {
    *count_out = 0;
    Expr* rr_call = expr_new_function(expr_new_symbol("RowReduce"),
        (Expr*[]){ expr_copy(M) }, 1);
    Expr* R = eval_and_free(rr_call);
    if (!R || R->type != EXPR_FUNCTION
        || R->data.function.head->type != EXPR_SYMBOL
        || R->data.function.head->data.symbol != SYM_List
        || (int)R->data.function.arg_count != n) {
        if (R) expr_free(R);
        return NULL;
    }

    bool* pivot_col = calloc(n, sizeof(bool));
    int* row_pivot = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        row_pivot[i] = -1;
        Expr* row = R->data.function.args[i];
        for (int j = 0; j < n; j++) {
            if (!is_zero_poly(row->data.function.args[j])) {
                row_pivot[i] = j;
                pivot_col[j] = true;
                break;
            }
        }
    }

    Expr** basis = malloc(sizeof(Expr*) * n);
    size_t bc = 0;
    for (int free_col = 0; free_col < n; free_col++) {
        if (pivot_col[free_col]) continue;
        Expr** vec = malloc(sizeof(Expr*) * n);
        for (int k = 0; k < n; k++) vec[k] = expr_new_integer(0);
        expr_free(vec[free_col]);
        vec[free_col] = expr_new_integer(1);
        for (int i = 0; i < n; i++) {
            int pc = row_pivot[i];
            if (pc < 0) continue;
            Expr* r_val =
                R->data.function.args[i]->data.function.args[free_col];
            Expr* neg = eval_and_free(expr_new_function(
                expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(r_val) }, 2));
            expr_free(vec[pc]);
            vec[pc] = neg;
        }
        basis[bc++] = expr_new_function(expr_new_symbol("List"), vec, n);
        free(vec);
    }

    free(pivot_col);
    free(row_pivot);
    expr_free(R);
    *count_out = bc;
    return basis;
}

/* Detect inexact (Real / MPFR) leaves anywhere in a matrix. */
static bool eigen_matrix_is_inexact(Expr* m) {
    if (!m) return false;
    if (m->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (m->type == EXPR_MPFR) return true;
#endif
    if (m->type != EXPR_FUNCTION) return false;
    if (eigen_matrix_is_inexact(m->data.function.head)) return true;
    for (size_t i = 0; i < m->data.function.arg_count; i++) {
        if (eigen_matrix_is_inexact(m->data.function.args[i])) return true;
    }
    return false;
}

/* Common option/positional argument parsing for Eigenvalues / Eigenvectors. */
typedef struct {
    Expr* arg0;        /* m or {m, a}                          */
    Expr* k_spec;      /* Integer k, or UpTo[k], or NULL       */
    bool  cubics_radical;
    bool  quartics_radical;
} EigenOpts;

/* True iff `e` is the symbol True. */
static bool eigen_is_true(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* Parse Eigenvalues/Eigenvectors arguments.  Returns false on shape error. */
static bool eigen_parse_args(Expr* res, EigenOpts* opts) {
    if (!res || res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return false;

    opts->arg0 = res->data.function.args[0];
    opts->k_spec = NULL;
    /* Default Cubics/Quartics: True so radicals are emitted by default --
     * essential for the numerical path where Root[] objects cannot be
     * numericalised.  Spec lists False as the default Solve option;
     * Eigenvalues overrides to keep the closed-form pipeline functional. */
    opts->cubics_radical = true;
    opts->quartics_radical = true;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 1) {
        Expr* a = res->data.function.args[pos_end - 1];
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol == SYM_Rule
                || a->data.function.head->data.symbol == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol;
            Expr* rhs = a->data.function.args[1];
            if (name == SYM_Cubics) {
                opts->cubics_radical = eigen_is_true(rhs);
                pos_end--; continue;
            }
            if (name == SYM_Quartics) {
                opts->quartics_radical = eigen_is_true(rhs);
                pos_end--; continue;
            }
            if (name == SYM_Method) {
                /* Reserved.  Currently no Method dispatch -- ignore. */
                pos_end--; continue;
            }
        }
        break;
    }
    if (pos_end > 2) return false;
    if (pos_end == 2) opts->k_spec = res->data.function.args[1];
    return true;
}

/* Apply k-spec (Integer k, -k, or UpTo[k]) to vals[0..count].  Returns a
 * freshly allocated trimmed array; the caller frees the originals it owns
 * for the values not selected.  *out_count holds the new size. */
static Expr** eigen_apply_k_spec(Expr** vals, size_t count, Expr* k_spec,
                                  size_t* out_count) {
    size_t result_count = count;
    bool from_end = false;
    if (k_spec) {
        if (k_spec->type == EXPR_INTEGER) {
            int64_t k = k_spec->data.integer;
            if (k >= 0) {
                result_count = ((size_t)k < count) ? (size_t)k : count;
            } else {
                int64_t abs_k = -k;
                result_count = ((size_t)abs_k < count) ? (size_t)abs_k : count;
                from_end = true;
            }
        } else if (k_spec->type == EXPR_FUNCTION
            && k_spec->data.function.head->type == EXPR_SYMBOL
            && k_spec->data.function.head->data.symbol == SYM_UpTo
            && k_spec->data.function.arg_count == 1
            && k_spec->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t k = k_spec->data.function.args[0]->data.integer;
            result_count = ((size_t)k < count) ? (size_t)k : count;
        }
    }
    Expr** result = result_count ? malloc(sizeof(Expr*) * result_count) : NULL;
    if (from_end) {
        size_t start = count - result_count;
        for (size_t i = 0; i < start; i++) expr_free(vals[i]);
        for (size_t i = 0; i < result_count; i++) result[i] = vals[start + i];
    } else {
        for (size_t i = 0; i < result_count; i++) result[i] = vals[i];
        for (size_t i = result_count; i < count; i++) expr_free(vals[i]);
    }
    *out_count = result_count;
    return result;
}

/* Compute the eigenvalue list (padded to n entries with Infinity for the
 * generalised degree-drop case).  Returns a freshly allocated array; the
 * caller owns the entries and the array.  *out_count is set to n on success
 * or to 0 if the calculation cannot be completed (caller should leave the
 * outer Eigenvalues call unevaluated). */
static Expr** eigen_compute_eigenvalues_full(Expr* m, Expr* a,
                                              int64_t n,
                                              bool cubics_radical,
                                              bool quartics_radical,
                                              size_t* out_count) {
    *out_count = 0;
    const char* lam = eigen_lambda_name();
    Expr* poly;
    if (a == NULL) {
        /* Ordinary case: Faddeev-Leverrier in O(n^4) matrix multiplications.
         * Far cheaper than Laplace expansion of the polynomial-entry
         * matrix det(m − λI) (which is O(n!)) once n grows past ~8. */
        poly = eigen_char_poly_faddeev(m, lam, (int)n);
        if (!poly) {
            return NULL;
        }
    } else {
        /* Generalised case: still use Laplace expansion of det(m − λa).
         * Acceptable in practice: generalised eigenproblems in the test
         * corpus are small (≤ 3×3). */
        Expr* M = eigen_build_lambda_matrix(m, a, lam, n);
        poly = eigen_compute_det(M, (int)n);
        expr_free(M);
    }

    Expr* sols = eigen_solve_poly(poly, lam,
                                  cubics_radical, quartics_radical);
    expr_free(poly);
    if (!sols) return NULL;

    size_t val_count = 0;
    Expr** vals = eigen_extract_values(sols, &val_count);
    expr_free(sols);

    /* For generalised eigenvalues, pad short result with Infinity. */
    if ((size_t)n > val_count) {
        Expr** padded = realloc(vals, sizeof(Expr*) * n);
        if (padded) vals = padded;
        else if (!vals) {
            vals = malloc(sizeof(Expr*) * n);
        }
        for (size_t i = val_count; i < (size_t)n; i++) {
            vals[i] = expr_new_symbol("Infinity");
        }
        val_count = (size_t)n;
    }
    *out_count = val_count;
    return vals;
}

Expr* builtin_eigenvalues(Expr* res) {
    EigenOpts opts;
    if (!eigen_parse_args(res, &opts)) return NULL;

    Expr *m, *a; int64_t n;
    if (!eigen_extract_matrix_pair(opts.arg0, &m, &a, &n)) return NULL;

    bool inexact = eigen_matrix_is_inexact(m)
                || (a && eigen_matrix_is_inexact(a));

    size_t val_count = 0;
    Expr** vals = eigen_compute_eigenvalues_full(m, a, n,
        opts.cubics_radical, opts.quartics_radical, &val_count);
    if (!vals || val_count == 0) {
        free(vals);
        return NULL;
    }

    /* For inexact input, chop numerical-noise imaginary parts so real
     * eigenvalues of real matrices appear as real numbers. */
    if (inexact) {
        for (size_t i = 0; i < val_count; i++) {
            Expr* c = eigen_chop(vals[i]);
            expr_free(vals[i]);
            vals[i] = c;
        }
    }

    eigen_sort_by_abs_desc(vals, val_count);

    size_t result_count = val_count;
    Expr** result_vals = eigen_apply_k_spec(vals, val_count, opts.k_spec,
                                             &result_count);
    free(vals);

    Expr* out = expr_new_function(expr_new_symbol("List"),
                                  result_vals, result_count);
    free(result_vals);
    return out;
}

/* ---------------- Eigenvectors ---------------- */

/* Substitute the lambda symbol in `M` with `val`, evaluating each entry.
 * Returns a freshly allocated matrix. */
static Expr* eigen_subst_lambda(Expr* M, const char* lambda_name, Expr* val) {
    /* Use ReplaceAll: M /. lambda_name -> val */
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol(lambda_name), expr_copy(val) }, 2);
    Expr* replaced = eval_and_free(expr_new_function(
        expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(M), rule }, 2));
    return replaced;
}

/* Build the m - λ*a matrix, then substitute λ = `val`.  Returns a numerical
 * (or symbolic-residual) n×n matrix.  Used by the eigenvector routine. */
static Expr* eigen_residual_matrix(Expr* m, Expr* a, Expr* val, int64_t n) {
    /* For ordinary case: just m - val*I, computed directly.
     * For generalised case: m - val*a. */
    const char* lam = eigen_lambda_name();
    Expr* M_lambda = eigen_build_lambda_matrix(m, a, lam, n);
    Expr* M_sub = eigen_subst_lambda(M_lambda, lam, val);
    expr_free(M_lambda);
    /* Try to canonicalise: expand each entry so RowReduce sees them in a
     * normalised form.  Many tests rely on integer-like cancellation here
     * (e.g. (7/2 - 4)*v1 + (1/2)*v3 -> -1/2 v1 + 1/2 v3). */
    return M_sub;
}

/* Build a length-n vector of zeros. */
static Expr* eigen_zero_vector(int64_t n) {
    Expr** v = malloc(sizeof(Expr*) * n);
    for (int64_t i = 0; i < n; i++) v[i] = expr_new_integer(0);
    Expr* out = expr_new_function(expr_new_symbol("List"), v, (size_t)n);
    free(v);
    return out;
}

/* Normalise vector `v` (a List of n elements) by dividing by its Norm.
 * Used by the numerical eigenvector path to emit unit vectors. */
static Expr* eigen_normalize_vector(Expr* v) {
    Expr* norm = eval_and_free(expr_new_function(
        expr_new_symbol("Norm"), (Expr*[]){ expr_copy(v) }, 1));
    /* If Norm is zero or symbolic, skip normalisation. */
    if (is_zero_poly(norm)) { expr_free(norm); return expr_copy(v); }
    Expr* inv = eval_and_free(expr_new_function(
        expr_new_symbol("Power"),
        (Expr*[]){ norm, expr_new_integer(-1) }, 2));
    Expr* scaled = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){ inv, expr_copy(v) }, 2));
    return scaled;
}

Expr* builtin_eigenvectors(Expr* res) {
    EigenOpts opts;
    if (!eigen_parse_args(res, &opts)) return NULL;

    Expr *m, *a; int64_t n;
    if (!eigen_extract_matrix_pair(opts.arg0, &m, &a, &n)) return NULL;

    bool inexact = eigen_matrix_is_inexact(m)
                || (a && eigen_matrix_is_inexact(a));

    /* Inexact input: rationalise once up front so we can perform the
     * eigenvalue / null-space arithmetic in exact form -- numerical
     * RowReduce on a (m − λI) substituted with inexact λ would zero out
     * the rank-defect we need to discover the eigenvector.  Numericalize
     * + normalise the result at the end. */
    Expr* m_orig = m;
    Expr* a_orig = a;
    long prec_bits = 53;
    Expr* m_rat = NULL;
    Expr* a_rat = NULL;
    if (inexact) {
        CommonInexactInfo info = common_scan_inexact(m);
        if (a) {
            CommonInexactInfo info_a = common_scan_inexact(a);
            if (info_a.has_inexact
                && (!info.has_inexact || info_a.min_bits < info.min_bits)) {
                info = info_a;
            }
        }
        prec_bits = info.min_bits ? info.min_bits : 53;
        m_rat = common_rationalize_input(m, prec_bits);
        if (a) a_rat = common_rationalize_input(a, prec_bits);
        m = m_rat;
        if (a) a = a_rat;
    }

    /* Compute eigenvalues in the same arithmetic form as the matrix. */
    size_t val_count = 0;
    Expr** vals = eigen_compute_eigenvalues_full(m, a, n,
        opts.cubics_radical, opts.quartics_radical, &val_count);
    if (!vals || val_count == 0) {
        free(vals);
        if (m_rat) expr_free(m_rat);
        if (a_rat) expr_free(a_rat);
        return NULL;
    }
    eigen_sort_by_abs_desc(vals, val_count);

    /* Collect eigenvectors -- traverse eigenvalues in order, collapsing
     * runs of equal values into a single null-space computation that
     * yields up to `multiplicity` vectors. */
    Expr** vectors = malloc(sizeof(Expr*) * n);
    size_t vc = 0;

    size_t i = 0;
    while (i < val_count && vc < (size_t)n) {
        /* Determine run of equal eigenvalues. */
        size_t j = i + 1;
        while (j < val_count && expr_eq(vals[j], vals[i])) j++;
        int64_t mult = (int64_t)(j - i);

        Expr* val = vals[i];

        /* Special handling: Infinity eigenvalues correspond to the null
         * space of `a` (the generalised pencil's "infinite" branch). */
        bool is_inf = (val->type == EXPR_SYMBOL
                       && val->data.symbol == SYM_Infinity);

        Expr* residual = NULL;
        if (is_inf && a) {
            residual = expr_copy(a);
        } else {
            residual = eigen_residual_matrix(m, a, val, n);
        }

        size_t basis_count = 0;
        Expr** basis = eigen_null_space(residual, (int)n, &basis_count);
        expr_free(residual);

        /* Take up to `mult` vectors. */
        size_t take = (basis_count < (size_t)mult) ? basis_count : (size_t)mult;
        for (size_t k = 0; k < take && vc < (size_t)n; k++) {
            Expr* v = basis[k];
            if (inexact) v = eigen_normalize_vector(v);
            else v = expr_copy(v);
            vectors[vc++] = v;
        }
        /* If null-space gives fewer vectors than multiplicity, the matrix
         * is defective for this eigenvalue: pad the shortfall in-line with
         * zero vectors so the i-th eigenvector still lines up positionally
         * with the i-th eigenvalue. */
        for (size_t k = take; k < (size_t)mult && vc < (size_t)n; k++) {
            vectors[vc++] = eigen_zero_vector(n);
        }
        for (size_t k = 0; k < basis_count; k++) expr_free(basis[k]);
        free(basis);

        i = j;
    }

    /* Pad with zero vectors. */
    while (vc < (size_t)n) {
        vectors[vc++] = eigen_zero_vector(n);
    }

    /* Free unused eigenvalues. */
    for (size_t k = 0; k < val_count; k++) expr_free(vals[k]);
    free(vals);

    size_t result_count = vc;
    Expr** result_vecs = eigen_apply_k_spec(vectors, vc, opts.k_spec,
                                             &result_count);
    free(vectors);

    Expr* out = expr_new_function(expr_new_symbol("List"),
                                  result_vecs, result_count);
    free(result_vecs);

    /* Inexact path: numericalize the eigenvectors back to the original
     * precision.  We did all of the null-space and normalisation work in
     * exact (rational) form so the rank-deficient direction was correctly
     * captured.  Free the rationalised matrices we created. */
    if (inexact) {
        Expr* numeric = common_numericalize_result(out, prec_bits);
        expr_free(out);
        out = numeric;
        if (m_rat) expr_free(m_rat);
        if (a_rat) expr_free(a_rat);
    }
    (void)m_orig; (void)a_orig;
    return out;
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
    /* RowReduce and LinearSolve live in src/matsol.c; registered via matsol_init(). */
    symtab_add_builtin("IdentityMatrix", builtin_identitymatrix);
    symtab_get_def("IdentityMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DiagonalMatrix", builtin_diagonalmatrix);
    symtab_get_def("DiagonalMatrix")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Inverse", builtin_inverse);
    symtab_get_def("Inverse")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatrixPower", builtin_matrixpower);
    symtab_get_def("MatrixPower")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Eigenvalues", builtin_eigenvalues);
    symtab_get_def("Eigenvalues")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Eigenvectors", builtin_eigenvectors);
    symtab_get_def("Eigenvectors")->attributes |= ATTR_PROTECTED;
}
