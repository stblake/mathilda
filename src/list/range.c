#include "list_common.h"
#include "range.h"

Expr* builtin_range(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;
    
    size_t len = res->data.function.arg_count;
    Expr* imin_e = NULL;
    Expr* imax_e = NULL;
    Expr* di_e = NULL;
    
    if (len == 1) {
        imin_e = expr_new_integer(1);
        imax_e = expr_copy(res->data.function.args[0]);
        di_e = expr_new_integer(1);
    } else if (len == 2) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_new_integer(1);
    } else if (len == 3) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_copy(res->data.function.args[2]);
    }
    
    bool is_real = false;
    double min_val = 0, max_val = 0, di_val = 0;
    int64_t n, d;
    
    if (imin_e->type == EXPR_REAL || imax_e->type == EXPR_REAL || di_e->type == EXPR_REAL) is_real = true;
    
    if (imin_e->type == EXPR_INTEGER) min_val = (double)imin_e->data.integer;
    else if (imin_e->type == EXPR_REAL) min_val = imin_e->data.real;
    else if (is_rational(imin_e, &n, &d)) min_val = (double)n / d;
    else goto L_fail_range;
    
    if (imax_e->type == EXPR_INTEGER) max_val = (double)imax_e->data.integer;
    else if (imax_e->type == EXPR_REAL) max_val = imax_e->data.real;
    else if (is_rational(imax_e, &n, &d)) max_val = (double)n / d;
    else goto L_fail_range;
    
    if (di_e->type == EXPR_INTEGER) di_val = (double)di_e->data.integer;
    else if (di_e->type == EXPR_REAL) di_val = di_e->data.real;
    else if (is_rational(di_e, &n, &d)) di_val = (double)n / d;
    else goto L_fail_range;
    
    if (di_val == 0) goto L_fail_range;
    if ((di_val > 0 && min_val > max_val) || (di_val < 0 && min_val < max_val)) {
        expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    
    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);
    
    double val = min_val;
    int steps = 0;
    Expr* curr_e = expr_copy(imin_e);
    
    while ((di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
        Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
        
        if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
        results[results_count++] = i_val;
        
        Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
        Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
        Expr* next_e = evaluate(next_expr);
        expr_free(next_expr);
        expr_free(curr_e);
        curr_e = next_e;
        
        val += di_val;
        steps++;
        if (steps > 1000000) break; 
    }
    
    if (curr_e) expr_free(curr_e);
    expr_free(imin_e);
    expr_free(imax_e);
    expr_free(di_e);
    
    Expr* result_list = expr_new_function(expr_new_symbol(SYM_List), results, results_count);
    free(results);
    return result_list;

L_fail_range:
    if (imin_e) expr_free(imin_e);
    if (imax_e) expr_free(imax_e);
    if (di_e) expr_free(di_e);
    return NULL;
}
