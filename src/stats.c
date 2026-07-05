#include "print.h"
#include "stats.h"
#include "list.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "complex.h"
#include "sym_names.h"
#include "assoc.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

static bool is_numeric(Expr* e, double* val, bool* complex) {
    if (e->type == EXPR_INTEGER) {
        if (val) *val = (double)e->data.integer;
        if (complex) *complex = false;
        return true;
    }
    if (e->type == EXPR_REAL) {
        if (val) *val = e->data.real;
        if (complex) *complex = false;
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        if (val) *val = (double)n / (double)d;
        if (complex) *complex = false;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (complex) *complex = true;
        return true;
    }
    return false;
}

static Expr* apply_columnwise(const char* func_name, Expr* matrix) {
    // Result is Map[func_name, Transpose[matrix]]
    Expr* transpose_args[1] = { expr_copy(matrix) };
    Expr* transpose_call = expr_new_function(expr_new_symbol(SYM_Transpose), transpose_args, 1);
    Expr* transposed = evaluate(transpose_call);
    expr_free(transpose_call);

    if (transposed->type != EXPR_FUNCTION) {
        expr_free(transposed);
        return NULL;
    }

    Expr* map_args[2] = { expr_new_symbol(func_name), transposed };
    Expr* map_call = expr_new_function(expr_new_symbol(SYM_Map), map_args, 2);
    Expr* result = evaluate(map_call);
    expr_free(map_call);
    return result;
}

Expr* builtin_mean(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Mean[assoc] is the mean of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol == SYM_True) {
        expr_free(is_matrix);
        return apply_columnwise("Mean", data);
    }
    expr_free(is_matrix);

    // Check if it's a vector (list)
    Expr* listq_args[1] = { expr_copy(data) };
    Expr* listq_call = expr_new_function(expr_new_symbol(SYM_ListQ), listq_args, 1);
    Expr* is_list = evaluate(listq_call);
    expr_free(listq_call);
    if (is_list->type == EXPR_SYMBOL && is_list->data.symbol == SYM_False) {
        expr_free(is_list);
        return NULL;
    }
    expr_free(is_list);

    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;

    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (elem->type == EXPR_REAL) {
            has_real = true;
        } else if (elem->type == EXPR_INTEGER) {
            // Keep as rational n/1
        } else if (is_rational(elem, NULL, NULL)) {
            // Keep as rational
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        if (has_real) {
            double sum_val = 0;
            for (size_t i = 0; i < n; i++) {
                double v = 0.0;
                is_numeric(data->data.function.args[i], &v, NULL);
                sum_val += v;
            }
            return expr_new_real(sum_val / (double)n);
        } else {
            // Exact rational arithmetic
            int64_t sum_n = 0;
            int64_t sum_d = 1;
            for (size_t i = 0; i < n; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) {
                    cur_n = elem->data.integer;
                    cur_d = 1;
                } else {
                    is_rational(elem, &cur_n, &cur_d);
                }
                // sum = sum_n/sum_d + cur_n/cur_d = (sum_n*cur_d + cur_n*sum_d) / (sum_d*cur_d)
                int64_t new_n = sum_n * cur_d + cur_n * sum_d;
                int64_t new_d = sum_d * cur_d;
                int64_t common = gcd(new_n < 0 ? -new_n : new_n, new_d);
                sum_n = new_n / common;
                sum_d = new_d / common;
            }
            // mean = sum / n = sum_n / (sum_d * n)
            return make_rational(sum_n, sum_d * (int64_t)n);
        }
    }

    // Fallback to symbolic: Mean = (1/n) * Plus @@ data
    Expr* sum_call = expr_new_function(expr_new_symbol(SYM_Apply), (Expr*[]){expr_new_symbol(SYM_Plus), expr_copy(data)}, 2);
    Expr* sum = evaluate(sum_call);
    expr_free(sum_call);

    Expr* n_inv = make_rational(1, (int64_t)n);
    Expr* mean_args[2] = { n_inv, sum };
    Expr* mean_call = expr_new_function(expr_new_symbol(SYM_Times), mean_args, 2);
    Expr* result = evaluate(mean_call);
    expr_free(mean_call);
    return result;
}

Expr* builtin_rootmeansquare(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol == SYM_True) {
        expr_free(is_matrix);
        return apply_columnwise("RootMeanSquare", data);
    }
    expr_free(is_matrix);

    // Check if it's a vector (list)
    Expr* listq_args[1] = { expr_copy(data) };
    Expr* listq_call = expr_new_function(expr_new_symbol(SYM_ListQ), listq_args, 1);
    Expr* is_list = evaluate(listq_call);
    expr_free(listq_call);
    if (is_list->type == EXPR_SYMBOL && is_list->data.symbol == SYM_False) {
        expr_free(is_list);
        return NULL;
    }
    expr_free(is_list);

    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;

    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (elem->type == EXPR_REAL) {
            has_real = true;
        } else if (elem->type == EXPR_INTEGER) {
            // Keep as rational n/1
        } else if (is_rational(elem, NULL, NULL)) {
            // Keep as rational
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric && has_real) {
        double sum_sq = 0;
        for (size_t i = 0; i < n; i++) {
            double v = 0.0;
            is_numeric(data->data.function.args[i], &v, NULL);
            sum_sq += v * v;
        }
        return expr_new_real(sqrt(sum_sq / (double)n));
    }

    // Fallback to symbolic (also used for exact rationals/integers to ensure correct square root distribution)
    Expr** sq_args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr* power_args[2] = { expr_copy(data->data.function.args[i]), expr_new_integer(2) };
        Expr* power_call = expr_new_function(expr_new_symbol(SYM_Power), power_args, 2);
        sq_args[i] = evaluate(power_call);
        expr_free(power_call);
    }
    
    Expr* sum_call = expr_new_function(expr_new_symbol(SYM_Plus), sq_args, n);
    Expr* sum = evaluate(sum_call);
    expr_free(sum_call);
    free(sq_args);
    
    bool sum_is_numeric = (sum->type == EXPR_INTEGER || sum->type == EXPR_REAL || sum->type == EXPR_BIGINT || is_rational(sum, NULL, NULL));
    
    int64_t root_n = (int64_t)round(sqrt((double)n));
    if (!sum_is_numeric && root_n * root_n == (int64_t)n) {
        Expr* sqrt_sum_args[2] = { sum, make_rational(1, 2) };
        Expr* sqrt_sum_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_sum_args, 2);
        Expr* sqrt_sum = evaluate(sqrt_sum_call);
        expr_free(sqrt_sum_call);

        Expr* n_inv = make_rational(1, root_n);
        Expr* times_args[2] = { n_inv, sqrt_sum };
        Expr* times_call = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
        Expr* result = evaluate(times_call);
        expr_free(times_call);
        return result;
    } else {
        Expr* n_inv = make_rational(1, (int64_t)n);
        Expr* mean_sq_args[2] = { n_inv, sum };
        Expr* mean_sq_call = expr_new_function(expr_new_symbol(SYM_Times), mean_sq_args, 2);
        Expr* mean_sq = evaluate(mean_sq_call);
        expr_free(mean_sq_call);
        
        int64_t num, den;
        if (is_rational(mean_sq, &num, &den)) {
            int64_t root_den = (int64_t)round(sqrt((double)den));
            if (root_den * root_den == den && root_den > 1) {
                Expr* sqrt_num_args[2] = { expr_new_integer(num), make_rational(1, 2) };
                Expr* sqrt_num_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_num_args, 2);
                Expr* sqrt_num = evaluate(sqrt_num_call);
                expr_free(sqrt_num_call);
                
                Expr* times_args[2] = { make_rational(1, root_den), sqrt_num };
                Expr* times_call = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
                Expr* result = evaluate(times_call);
                expr_free(times_call);
                expr_free(mean_sq);
                return result;
            }
        }

        Expr* sqrt_args[2] = { mean_sq, make_rational(1, 2) };
        Expr* sqrt_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_args, 2);
        Expr* result = evaluate(sqrt_call);
        expr_free(sqrt_call);
        return result;
    }
}

Expr* builtin_variance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Variance[assoc] is the variance of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol == SYM_True) {
        expr_free(is_matrix);
        return apply_columnwise("Variance", data);
    }
    expr_free(is_matrix);

    // Check if it's a vector (list)
    Expr* listq_args[1] = { expr_copy(data) };
    Expr* listq_call = expr_new_function(expr_new_symbol(SYM_ListQ), listq_args, 1);
    Expr* is_list = evaluate(listq_call);
    expr_free(listq_call);
    if (is_list->type == EXPR_SYMBOL && is_list->data.symbol == SYM_False) {
        expr_free(is_list);
        return NULL;
    }
    expr_free(is_list);

    size_t n = data->data.function.arg_count;
    if (n <= 1) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (elem->type == EXPR_REAL) {
            has_real = true;
        } else if (elem->type == EXPR_INTEGER || is_rational(elem, NULL, NULL)) {
            // Numeric
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        if (has_real) {
            // Welford's algorithm
            double m = 0;
            double s = 0;
            for (size_t i = 0; i < n; i++) {
                double x = 0.0;
                is_numeric(data->data.function.args[i], &x, NULL);
                double old_m = m;
                m += (x - m) / (double)(i + 1);
                s += (x - m) * (x - old_m);
            }
            return expr_new_real(s / (double)(n - 1));
        } else {
            // Exact calculation for Variance
            // Var = (Sum[x^2] - n*Mean[x]^2) / (n-1)
            // Using Sum[(n*x_i - Sum[x_j])^2] / (n^2 * (n-1))
            int64_t sum_n = 0;
            int64_t sum_d = 1;
            for (size_t i = 0; i < n; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) { cur_n = elem->data.integer; cur_d = 1; }
                else is_rational(elem, &cur_n, &cur_d);
                int64_t new_n = sum_n * cur_d + cur_n * sum_d;
                int64_t new_d = sum_d * cur_d;
                int64_t common = gcd(new_n < 0 ? -new_n : new_n, new_d);
                sum_n = new_n / common;
                sum_d = new_d / common;
            }
            // Now sum_n/sum_d is the sum of elements
            int64_t sq_sum_n = 0;
            int64_t sq_sum_d = 1;
            for (size_t i = 0; i < n; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) { cur_n = elem->data.integer; cur_d = 1; }
                else is_rational(elem, &cur_n, &cur_d);
                
                // (x - mean)^2 = (cur_n/cur_d - sum_n/(n*sum_d))^2 
                // = ( (cur_n * n * sum_d - sum_n * cur_d) / (n * sum_d * cur_d) )^2
                int64_t term_n = cur_n * (int64_t)n * sum_d - sum_n * cur_d;
                int64_t term_d = (int64_t)n * sum_d * cur_d;
                int64_t common = gcd(term_n < 0 ? -term_n : term_n, term_d);
                term_n /= common; term_d /= common;
                
                int64_t term_sq_n = term_n * term_n;
                int64_t term_sq_d = term_d * term_d;
                
                int64_t new_sq_sum_n = sq_sum_n * term_sq_d + term_sq_n * sq_sum_d;
                int64_t new_sq_sum_d = sq_sum_d * term_sq_d;
                common = gcd(new_sq_sum_n < 0 ? -new_sq_sum_n : new_sq_sum_n, new_sq_sum_d);
                sq_sum_n = new_sq_sum_n / common;
                sq_sum_d = new_sq_sum_d / common;
            }
            // Variance = sq_sum / (n-1)
            return make_rational(sq_sum_n, sq_sum_d * ((int64_t)n - 1));
        }
    }

    // Fallback to symbolic: Compute Mean
    Expr* mean_args[1] = { expr_copy(data) };
    Expr* mean_call = expr_new_function(expr_new_symbol(SYM_Mean), mean_args, 1);
    Expr* mu = evaluate(mean_call);
    expr_free(mean_call);

    // Sum[(x - mu) * Conjugate[x - mu]]
    Expr** diffs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* x = data->data.function.args[i];
        Expr* sub_args[2] = { expr_copy(x), expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(mu)}, 2) };
        Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), sub_args, 2));
        
        Expr* conj_args[1] = { expr_copy(diff) };
        Expr* conj_call = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
        Expr* conj_diff = evaluate(conj_call);
        expr_free(conj_call);

        Expr* prod_args[2] = { diff, conj_diff };
        diffs[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), prod_args, 2));
    }
    expr_free(mu);

    Expr* sum_diffs = expr_new_function(expr_new_symbol(SYM_Plus), diffs, n);
    Expr* sum = evaluate(sum_diffs);
    expr_free(sum_diffs);
    free(diffs);

    Expr* n_minus_1_inv = make_rational(1, (int64_t)n - 1);
    Expr* var_call = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ n_minus_1_inv, sum }, 2);
    Expr* result = evaluate(var_call);
    expr_free(var_call);
    return result;
}

Expr* builtin_standard_deviation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* StandardDeviation[assoc] uses the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol == SYM_True) {
        expr_free(is_matrix);
        return apply_columnwise("StandardDeviation", data);
    }
    expr_free(is_matrix);

    // Optimization for real numeric data
    if (data->type == EXPR_FUNCTION) {
        bool all_real = true;
        size_t n = data->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            bool complex;
            if (!is_numeric(data->data.function.args[i], NULL, &complex) || complex) {
                all_real = false;
                break;
            }
        }
        if (all_real && n > 1) {
            Expr* var_call = expr_new_function(expr_new_symbol(SYM_Variance), (Expr*[]){ expr_copy(data) }, 1);
            Expr* var_e = evaluate(var_call);
            expr_free(var_call);
            if (var_e->type == EXPR_REAL) {
                double val = sqrt(var_e->data.real);
                expr_free(var_e);
                return expr_new_real(val);
            }
            expr_free(var_e);
        }
    }

    // StandardDeviation[data] = Variance[data]^(1/2)
    Expr* var_call = expr_new_function(expr_new_symbol(SYM_Variance), (Expr*[]){ expr_copy(data) }, 1);
    Expr* var = evaluate(var_call);
    expr_free(var_call);

    if (!var) return NULL;

    Expr* sqrt_call = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ var, make_rational(1, 2) }, 2);
    Expr* result = evaluate(sqrt_call);
    expr_free(sqrt_call);
    return result;
}



/* ------------------- Median ------------------- */

static bool is_real_numeric(Expr* e) {
    Expr* numq = expr_new_function(expr_new_symbol(SYM_NumericQ), (Expr*[]){expr_copy(e)}, 1);
    Expr* numq_eval = evaluate(numq);
    expr_free(numq);
    if (numq_eval->type != EXPR_SYMBOL || numq_eval->data.symbol != SYM_True) {
        expr_free(numq_eval);
        return false;
    }
    expr_free(numq_eval);
    
    Expr* freeq = expr_new_function(expr_new_symbol(SYM_FreeQ), (Expr*[]){expr_copy(e), expr_new_symbol(SYM_I)}, 2);
    Expr* freeq_eval = evaluate(freeq);
    expr_free(freeq);
    if (freeq_eval->type != EXPR_SYMBOL || freeq_eval->data.symbol != SYM_True) {
        expr_free(freeq_eval);
        return false;
    }
    expr_free(freeq_eval);
    
    return true;
}

Expr* builtin_median(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Median[assoc] is the median of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    // Check if it's a vector or tensor. If it's empty or not a list, return NULL.
    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL || data->data.function.head->data.symbol != SYM_List) {
        return expr_copy(res);
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) return expr_copy(res);

    // Check if it's a matrix/tensor by checking if the first element is a List.
    if (data->data.function.args[0]->type == EXPR_FUNCTION && 
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL && 
        data->data.function.args[0]->data.function.head->data.symbol == SYM_List) {
        // Apply columnwise via Transpose and Map
        return apply_columnwise("Median", data);
    }

    // Now it's a 1D vector. Check if all elements are real numbers.
    bool all_real = true;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (!is_real_numeric(elem)) {
            all_real = false;
            break;
        }
    }

    if (!all_real) {
        char* str = expr_to_string(res);
        printf("Median::rectn: Rectangular array of real numbers is expected at position 1 in %s.\n", str);
        free(str);
        return expr_copy(res);
    }

    // Sort the list
    Expr* sort_expr = expr_new_function(expr_new_symbol(SYM_Sort), (Expr*[]){expr_copy(data)}, 1);
    Expr* sorted = evaluate(sort_expr);
    expr_free(sort_expr);

    if (sorted->type != EXPR_FUNCTION || sorted->data.function.arg_count != n) {
        expr_free(sorted);
        return expr_copy(res);
    }

    Expr* result = NULL;
    if (n % 2 == 1) {
        result = expr_copy(sorted->data.function.args[n / 2]);
    } else {
        Expr* m1 = sorted->data.function.args[n / 2 - 1];
        Expr* m2 = sorted->data.function.args[n / 2];
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(m1), expr_copy(m2)}, 2);
        Expr* sum_eval = evaluate(sum);
        expr_free(sum);
        
        Expr* div = expr_new_function(expr_new_symbol(SYM_Divide), (Expr*[]){sum_eval, expr_new_integer(2)}, 2);
        result = evaluate(div);
        expr_free(div);
    }

    expr_free(sorted);
    return result;
}


/* ------------------- Quartiles ------------------- */

Expr* builtin_quartiles(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* param_expr = NULL;
    if (res->data.function.arg_count == 2) {
        param_expr = res->data.function.args[1];
    }

    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL || data->data.function.head->data.symbol != SYM_List) {
        return expr_copy(res);
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) return expr_copy(res);

    if (data->data.function.args[0]->type == EXPR_FUNCTION && 
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL && 
        data->data.function.args[0]->data.function.head->data.symbol == SYM_List) {
        
        Expr* t_expr = expr_new_function(expr_new_symbol(SYM_Transpose), (Expr*[]){expr_copy(data)}, 1);
        Expr* t_eval = evaluate(t_expr);
        expr_free(t_expr);
        
        if (t_eval->type != EXPR_FUNCTION) {
            expr_free(t_eval);
            return expr_copy(res);
        }
        
        size_t cols = t_eval->data.function.arg_count;
        Expr** q_args = malloc(sizeof(Expr*) * cols);
        for (size_t i = 0; i < cols; i++) {
            Expr* call_args[2];
            call_args[0] = expr_copy(t_eval->data.function.args[i]);
            size_t call_cnt = 1;
            if (param_expr) {
                call_args[1] = expr_copy(param_expr);
                call_cnt = 2;
            }
            Expr* q_call = expr_new_function(expr_new_symbol(SYM_Quartiles), call_args, call_cnt);
            q_args[i] = evaluate(q_call);
            expr_free(q_call);
        }
        Expr* res_list = expr_new_function(expr_new_symbol(SYM_List), q_args, cols);
        free(q_args);
        expr_free(t_eval);
        return res_list;
    }

    bool all_real = true;
    for (size_t i = 0; i < n; i++) {
        if (!is_real_numeric(data->data.function.args[i])) {
            all_real = false;
            break;
        }
    }

    if (!all_real) {
        char* str = expr_to_string(res);
        printf("Quartiles::rectn: Rectangular array of real numbers is expected at position 1 in %s.\n", str);
        free(str);
        return expr_copy(res);
    }

    Expr *a = NULL, *b = NULL, *c = NULL, *d = NULL;
    if (param_expr) {
        if (param_expr->type == EXPR_FUNCTION && param_expr->data.function.arg_count == 2) {
            Expr* row1 = param_expr->data.function.args[0];
            Expr* row2 = param_expr->data.function.args[1];
            if (row1->type == EXPR_FUNCTION && row1->data.function.arg_count == 2 &&
                row2->type == EXPR_FUNCTION && row2->data.function.arg_count == 2) {
                a = expr_copy(row1->data.function.args[0]);
                b = expr_copy(row1->data.function.args[1]);
                c = expr_copy(row2->data.function.args[0]);
                d = expr_copy(row2->data.function.args[1]);
            }
        }
    } 
    
    if (!a || !b || !c || !d) {
        if (a) expr_free(a);
        if (b) expr_free(b);
        if (c) expr_free(c);
        if (d) expr_free(d);
        a = make_rational(1, 2);
        b = expr_new_integer(0);
        c = expr_new_integer(0);
        d = expr_new_integer(1);
    }

    Expr* sort_expr = expr_new_function(expr_new_symbol(SYM_Sort), (Expr*[]){expr_copy(data)}, 1);
    Expr* sorted = evaluate(sort_expr);
    expr_free(sort_expr);

    if (sorted->type != EXPR_FUNCTION || sorted->data.function.arg_count != n) {
        expr_free(sorted);
        expr_free(a); expr_free(b); expr_free(c); expr_free(d);
        return expr_copy(res);
    }

    Expr** sorted_A = sorted->data.function.args;
    Expr* q_vals[3];
    q_vals[0] = make_rational(1, 4);
    q_vals[1] = make_rational(1, 2);
    q_vals[2] = make_rational(3, 4);

    Expr* results[3];
    for (int k = 0; k < 3; k++) {
        Expr* q = q_vals[k];
        Expr* n_expr = expr_new_integer(n);
        Expr* n_plus_b = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){n_expr, expr_copy(b)}, 2));
        Expr* times_q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){n_plus_b, expr_copy(q)}, 2));
        Expr* h = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(a), times_q}, 2));

        double h_val = 0;
        if (!is_numeric(h, &h_val, NULL)) {
            results[k] = expr_new_symbol(SYM_Indeterminate);
            expr_free(h);
            continue;
        }

        if (h_val <= 1.0) {
            results[k] = expr_copy(sorted_A[0]);
            expr_free(h);
            continue;
        }
        if (h_val >= (double)n) {
            results[k] = expr_copy(sorted_A[n - 1]);
            expr_free(h);
            continue;
        }

        Expr* j_expr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Floor), (Expr*[]){expr_copy(h)}, 1));
        int64_t j_idx = 0;
        if (j_expr->type == EXPR_INTEGER) j_idx = j_expr->data.integer;
        else j_idx = (int64_t)floor(h_val);
        expr_free(j_expr);

        if (j_idx < 1) j_idx = 1;
        if (j_idx >= (int64_t)n) j_idx = n - 1;

        Expr* j_expr2 = expr_new_integer(j_idx);
        Expr* neg_j = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), j_expr2}, 2));
        Expr* g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(h), neg_j}, 2));
        expr_free(h);

        Expr* d_times_g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(d), expr_copy(g)}, 2));
        Expr* g_weight = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(c), d_times_g}, 2));
        expr_free(g);

        Expr* neg_Aj1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(sorted_A[j_idx-1])}, 2));
        Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(sorted_A[j_idx]), neg_Aj1}, 2));

        Expr* weight_diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){g_weight, diff}, 2));
        results[k] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(sorted_A[j_idx-1]), weight_diff}, 2));
    }

    expr_free(q_vals[0]);
    expr_free(q_vals[1]);
    expr_free(q_vals[2]);

    expr_free(sorted);
    expr_free(a);
    expr_free(b);
    expr_free(c);
    expr_free(d);

    return expr_new_function(expr_new_symbol(SYM_List), results, 3);
}

/* ------------------- MovingAverage ------------------- */

/*
 * MovingAverage[list, r]                — simple moving average over runs of r elements.
 * MovingAverage[list, {w1, ..., wr}]    — weighted moving average with weights w_k
 *                                         (output uses w_k / Sum[w_k] as effective weights).
 *
 * Output length is Length[list] - r + 1 when 1 <= r <= Length[list]; the function
 * stays unevaluated otherwise.  The unweighted form delegates to Mean so it inherits
 * Mean's exact rational / bigint / real / symbolic handling.  The weighted form builds
 * Plus[Times[w_k/wsum, x_{i+k}], ...] and lets the evaluator simplify.
 */
Expr* builtin_moving_average(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    size_t r = 0;
    Expr** weights = NULL; /* borrowed; non-NULL means weighted form */

    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 1) return NULL;
        if ((uint64_t)spec->data.integer > (uint64_t)n) return NULL;
        r = (size_t)spec->data.integer;
    } else if (spec->type == EXPR_BIGINT) {
        if (mpz_sgn(spec->data.bigint) <= 0) return NULL;
        if (!mpz_fits_ulong_p(spec->data.bigint)) return NULL;
        unsigned long rr = mpz_get_ui(spec->data.bigint);
        if ((size_t)rr > n) return NULL;
        r = (size_t)rr;
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol == SYM_List) {
        r = spec->data.function.arg_count;
        if (r == 0 || r > n) return NULL;
        weights = spec->data.function.args;
    } else {
        return NULL;
    }

    size_t out_n = n - r + 1;
    Expr** out = malloc(sizeof(Expr*) * out_n);
    if (!out) return NULL;

    if (weights == NULL) {
        for (size_t i = 0; i < out_n; i++) {
            Expr** sub = malloc(sizeof(Expr*) * r);
            for (size_t k = 0; k < r; k++) {
                sub[k] = expr_copy(data->data.function.args[i + k]);
            }
            Expr* sublist = expr_new_function(expr_new_symbol(SYM_List), sub, r);
            free(sub);
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Mean), (Expr*[]){ sublist }, 1));
        }
    } else {
        Expr** w_copies = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) w_copies[k] = expr_copy(weights[k]);
        Expr* wsum = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus), w_copies, r));
        free(w_copies);

        Expr* wsum_inv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power),
            (Expr*[]){ wsum, expr_new_integer(-1) }, 2));

        Expr** coefs = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) {
            coefs[k] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(weights[k]), expr_copy(wsum_inv) }, 2));
        }
        expr_free(wsum_inv);

        for (size_t i = 0; i < out_n; i++) {
            Expr** terms = malloc(sizeof(Expr*) * r);
            for (size_t k = 0; k < r; k++) {
                terms[k] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(coefs[k]),
                               expr_copy(data->data.function.args[i + k]) }, 2));
            }
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Plus), terms, r));
            free(terms);
        }

        for (size_t k = 0; k < r; k++) expr_free(coefs[k]);
        free(coefs);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, out_n);
    free(out);
    return result;
}

/* ------------------- MovingMedian ------------------- */

/*
 * MovingMedian[list, r]
 *   gives the moving median of list, computed using spans of r elements.
 *
 * Output length is Length[list] - r + 1 when 1 <= r <= Length[list]; the
 * function stays unevaluated otherwise. Operates on real-valued vectors and
 * matrices; for matrices each row-window is reduced via Median, which yields
 * a column-wise median vector. Non-numeric data triggers the MovingMedian::arg1
 * message and the expression remains unevaluated.
 */
Expr* builtin_moving_median(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    size_t r = 0;
    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 1) return NULL;
        if ((uint64_t)spec->data.integer > (uint64_t)n) return NULL;
        r = (size_t)spec->data.integer;
    } else if (spec->type == EXPR_BIGINT) {
        if (mpz_sgn(spec->data.bigint) <= 0) return NULL;
        if (!mpz_fits_ulong_p(spec->data.bigint)) return NULL;
        unsigned long rr = mpz_get_ui(spec->data.bigint);
        if ((size_t)rr > n) return NULL;
        r = (size_t)rr;
    } else {
        return NULL;
    }

    if (n == 0) return NULL;

    /* Decide vector vs matrix based on whether the first element is a List. */
    bool matrix_mode = (data->data.function.args[0]->type == EXPR_FUNCTION &&
                        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                        data->data.function.args[0]->data.function.head->data.symbol == SYM_List);

    /* Validate that every leaf is a real-valued numeric. Matrices must be rectangular. */
    bool ok = true;
    if (matrix_mode) {
        size_t cols = data->data.function.args[0]->data.function.arg_count;
        for (size_t i = 0; i < n && ok; i++) {
            Expr* row = data->data.function.args[i];
            if (row->type != EXPR_FUNCTION ||
                row->data.function.head->type != EXPR_SYMBOL ||
                row->data.function.head->data.symbol != SYM_List ||
                row->data.function.arg_count != cols) {
                ok = false;
                break;
            }
            for (size_t j = 0; j < cols && ok; j++) {
                if (!is_real_numeric(row->data.function.args[j])) ok = false;
            }
        }
    } else {
        for (size_t i = 0; i < n && ok; i++) {
            if (!is_real_numeric(data->data.function.args[i])) ok = false;
        }
    }

    if (!ok) {
        char* str = expr_to_string(data);
        printf("MovingMedian::arg1: The first argument %s must be a vector or matrix of real values.\n", str);
        free(str);
        return expr_copy(res);
    }

    size_t out_n = n - r + 1;
    Expr** out = malloc(sizeof(Expr*) * out_n);
    if (!out) return NULL;

    for (size_t i = 0; i < out_n; i++) {
        Expr** sub = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) {
            sub[k] = expr_copy(data->data.function.args[i + k]);
        }
        Expr* sublist = expr_new_function(expr_new_symbol(SYM_List), sub, r);
        free(sub);
        out[i] = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Median), (Expr*[]){ sublist }, 1));
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, out_n);
    free(out);
    return result;
}

/* ------------------- ExponentialMovingAverage ------------------- */

/*
 * ExponentialMovingAverage[list, alpha]
 *   gives the exponential moving average of list with smoothing constant alpha.
 *
 * The recurrence is y[1] = x[1], y[i+1] = y[i] + alpha * (x[i+1] - y[i]).
 * Output has the same length as the input list. Stays unevaluated when the
 * first argument is not a List, or when the list is empty.
 *
 * Two evaluation strategies:
 *   1. Fast path (machine precision): if at least one element of the list or
 *      alpha itself is a real number, we evaluate the recurrence in C using
 *      doubles. This is O(n) with no Expr allocations beyond the output and
 *      keeps numeric stability comparable to Mathematica's N[...] form.
 *   2. Symbolic path: build the recurrence as Plus / Times nodes and let the
 *      evaluator handle exact rational arithmetic, bignum promotion, and
 *      symbolic simplification. Works for arbitrary alpha (including symbolic).
 */
Expr* builtin_exponential_moving_average(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* alpha = res->data.function.args[1];

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    /* Decide whether the fast double-precision path applies: at least one
     * element of list or alpha must be EXPR_REAL, and every list element plus
     * alpha must be a real-valued numeric (Integer / Real / Rational, no
     * Complex, no symbolic). Bignums fall through to the symbolic path so we
     * don't lose precision in cases where the exact value is wanted. */
    bool any_real = (alpha->type == EXPR_REAL);
    if (!any_real) {
        for (size_t i = 0; i < n; i++) {
            if (data->data.function.args[i]->type == EXPR_REAL) {
                any_real = true;
                break;
            }
        }
    }
    bool fast_ok = any_real;
    if (fast_ok) {
        bool cplx = false;
        if (!is_numeric(alpha, NULL, &cplx) || cplx) fast_ok = false;
    }
    if (fast_ok) {
        for (size_t i = 0; i < n && fast_ok; i++) {
            bool cplx = false;
            if (!is_numeric(data->data.function.args[i], NULL, &cplx) || cplx) {
                fast_ok = false;
            }
        }
    }

    if (fast_ok) {
        double a = 0.0;
        is_numeric(alpha, &a, NULL);
        Expr** out = malloc(sizeof(Expr*) * n);
        if (!out) return NULL;
        double y = 0.0;
        is_numeric(data->data.function.args[0], &y, NULL);
        out[0] = expr_new_real(y);
        for (size_t i = 1; i < n; i++) {
            double x = 0.0;
            is_numeric(data->data.function.args[i], &x, NULL);
            y = y + a * (x - y);
            out[i] = expr_new_real(y);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        return result;
    }

    /* Symbolic / exact-rational / bignum path. */
    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;
    out[0] = expr_copy(data->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        Expr* x_i = expr_copy(data->data.function.args[i]);
        Expr* neg_y_prev = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), expr_copy(out[i-1]) }, 2));
        Expr* diff = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){ x_i, neg_y_prev }, 2));
        Expr* alpha_diff = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(alpha), diff }, 2));
        out[i] = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(out[i-1]), alpha_diff }, 2));
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);
    return result;
}

void stats_init(void) {
    symtab_add_builtin("Mean", builtin_mean);
    symtab_add_builtin("RootMeanSquare", builtin_rootmeansquare);
    symtab_add_builtin("Median", builtin_median);
    symtab_get_def("Median")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Quartiles", builtin_quartiles);
    symtab_get_def("Quartiles")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Variance", builtin_variance);
    symtab_add_builtin("StandardDeviation", builtin_standard_deviation);
    symtab_add_builtin("MovingAverage", builtin_moving_average);
    symtab_add_builtin("MovingMedian", builtin_moving_median);
    symtab_add_builtin("ExponentialMovingAverage", builtin_exponential_moving_average);

    symtab_get_def("Mean")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RootMeanSquare")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Variance")->attributes |= ATTR_PROTECTED;
    symtab_get_def("StandardDeviation")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MovingAverage")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MovingMedian")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ExponentialMovingAverage")->attributes |= ATTR_PROTECTED;
}
