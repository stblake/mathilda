#include "list_common.h"
#include "table.h"

Expr* builtin_table(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    
    if (res->data.function.arg_count > 2) {
        Expr** inner_args = malloc(sizeof(Expr*) * 2);
        inner_args[0] = expr_copy(res->data.function.args[0]);
        inner_args[1] = expr_copy(res->data.function.args[res->data.function.arg_count - 1]);
        Expr* inner_table = expr_new_function(expr_new_symbol(SYM_Table), inner_args, 2);
        free(inner_args);
        
        Expr** outer_args = malloc(sizeof(Expr*) * (res->data.function.arg_count - 1));
        outer_args[0] = inner_table;
        for (size_t i = 1; i < res->data.function.arg_count - 1; i++) {
            outer_args[i] = expr_copy(res->data.function.args[i]);
        }
        Expr* outer_table = expr_new_function(expr_new_symbol(SYM_Table), outer_args, res->data.function.arg_count - 1);
        free(outer_args);
        
        Expr* eval_outer = evaluate(outer_table);
        expr_free(outer_table);
        return eval_outer;
    }
    
    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* ---- Parse the iterator spec (shared helper) ---- */
    IterSpec s;
    if (!iter_spec_parse(spec, &s)) return NULL;

    int is_n_times   = (s.kind == ITER_KIND_COUNT);
    int is_list_iter = (s.kind == ITER_KIND_LIST);
    double min_val = 0, max_val = 0, di_val = 0;
    bool is_real = false, is_inf = false;

    /* Table does not iterate to Infinity; allow_inf = false. */
    if (!is_list_iter) {
        if (!iter_spec_resolve_numeric(&s, /*allow_inf=*/false,
                                       &min_val, &max_val, &di_val,
                                       &is_real, &is_inf)) {
            iter_spec_free(&s);
            return NULL;
        }
    }

    /* Convenience aliases into the owned IterSpec (freed via iter_spec_free). */
    Expr* var_sym = s.var;
    Expr* imin_e  = s.imin;
    Expr* imax_e  = s.imax;
    Expr* di_e    = s.di;
    Expr* list_e  = s.list;

    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);

    Rule* old_own = iter_spec_shadow(var_sym);

    if (is_n_times) {
        int64_t n = imax_e->data.integer;
        for (int64_t i = 0; i < n; i++) {
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else if (is_list_iter) {
        for (size_t i = 0; i < list_e->data.function.arg_count; i++) {
            symtab_add_own_value(var_sym->data.symbol.name, var_sym, list_e->data.function.args[i]);
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else {
        double val = min_val;
        int steps = 0;
        Expr* curr_e = expr_copy(imin_e);
        while ((di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
            Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
            symtab_add_own_value(var_sym->data.symbol.name, var_sym, i_val);
            
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
            
            expr_free(i_val);
            
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
    }

    iter_spec_restore(var_sym, old_own);
    iter_spec_free(&s);

    Expr* result_list = expr_new_function(expr_new_symbol(SYM_List), results, results_count);
    free(results);
    return result_list;
}
