#include "list_common.h"
#include "array.h"
#include "../pack.h"
#include "../checked_alloc.h"

static Expr* array_helper(Expr* f, Expr** n_array, Expr** r_array, size_t dim_count, size_t current_dim, Expr** current_args) {
    if (current_dim == dim_count) {
        Expr** fn_args = malloc(sizeof(Expr*) * dim_count);
        for (size_t i = 0; i < dim_count; i++) fn_args[i] = expr_copy(current_args[i]);
        Expr* fn_expr = expr_new_function(expr_copy(f), fn_args, dim_count);
        free(fn_args);
        Expr* eval_fn = evaluate(fn_expr);
        expr_free(fn_expr);
        return eval_fn;
    }

    Expr* n_expr = n_array[current_dim];
    Expr* r_expr = r_array[current_dim];
    
    int64_t n_val = n_expr->data.integer;
    Expr** results = malloc(sizeof(Expr*) * (size_t)n_val);
    if (n_val == 0) {
        free(results);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    
    int is_range = 0;
    Expr* a_expr = NULL;
    Expr* b_expr = NULL;
    
    if (r_expr && r_expr->type == EXPR_FUNCTION && r_expr->data.function.head->type == EXPR_SYMBOL && r_expr->data.function.head->data.symbol.name == SYM_List && r_expr->data.function.arg_count == 2) {
        is_range = 1;
        a_expr = r_expr->data.function.args[0];
        b_expr = r_expr->data.function.args[1];
    }
    
    Expr* r_base = NULL;
    if (!is_range) {
        r_base = r_expr ? expr_copy(r_expr) : expr_new_integer(1);
    }
    
    for (int64_t i = 0; i < n_val; i++) {
        Expr* arg = NULL;
        if (is_range) {
            if (n_val == 1) {
                arg = expr_copy(a_expr);
            } else {
                Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(b_expr), expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(a_expr)}, 2)}, 2);
                Expr* frac = expr_new_function(expr_new_symbol(SYM_Divide), (Expr*[]){expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(i), diff}, 2), expr_new_integer(n_val - 1)}, 2);
                Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(a_expr), frac}, 2);
                arg = evaluate(sum);
                expr_free(sum);
            }
        } else {
            if (i == 0) {
                arg = expr_copy(r_base);
            } else {
                Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(r_base), expr_new_integer(i)}, 2);
                arg = evaluate(sum);
                expr_free(sum);
            }
        }
        
        current_args[current_dim] = arg;
        results[i] = array_helper(f, n_array, r_array, dim_count, current_dim + 1, current_args);
        expr_free(arg);
    }
    
    if (r_base) expr_free(r_base);
    
    Expr* list_result = expr_new_function(expr_new_symbol(SYM_List), results, (size_t)n_val);
    free(results);
    return list_result;
}

/* Set while the 4-argument form builds its List skeleton: the result is not
 * offered to the packer, since its heads are about to be replaced. */
static bool g_array_no_pack = false;

/* Replace the List heads of the top `levels` levels of the array skeleton
 * `a` (owned) by h: Array[f, dims, origin, h] (audit Array#2). */
static Expr* array_rehead(Expr* a, Expr* h, size_t levels) {
    if (levels == 0 || !head_is(a, SYM_List)) return a;
    size_t n = a->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        items[i] = array_rehead(expr_copy(a->data.function.args[i]), h, levels - 1);
    Expr* out = expr_new_function(expr_copy(h), items, n);
    free(items);
    expr_free(a);
    return out;
}

Expr* builtin_array(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 4) return NULL;
    /* Array[f, dims, origin, h]: the array under head h at every level, as
     * h[f[1], f[2], ...] -- built as the List array, then re-headed. */
    if (res->data.function.arg_count == 4) {
        Expr* a3[3] = { expr_copy(res->data.function.args[0]), expr_copy(res->data.function.args[1]),
                        expr_copy(res->data.function.args[2]) };
        Expr* call = expr_new_function(expr_new_symbol(SYM_Array), a3, 3);
        g_array_no_pack = true;              /* a List skeleton to re-head */
        Expr* arr = builtin_array(call);
        g_array_no_pack = false;
        expr_free(call);
        if (!arr || eval_is_inflight_unwind(arr)) return arr;
        Expr* n_spec = res->data.function.args[1];
        size_t levels = head_is(n_spec, SYM_List) ? n_spec->data.function.arg_count : 1;
        return array_rehead(arr, res->data.function.args[3], levels);
    }
    
    Expr* f = res->data.function.args[0];
    Expr* n_spec = res->data.function.args[1];
    Expr* r_spec = res->data.function.arg_count == 3 ? res->data.function.args[2] : NULL;
    
    size_t dim_count = 1;
    if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->type == EXPR_SYMBOL && n_spec->data.function.head->data.symbol.name == SYM_List) {
        dim_count = n_spec->data.function.arg_count;
        if (dim_count == 0) return NULL;
    }
    
    Expr** n_array = malloc(sizeof(Expr*) * dim_count);
    Expr** r_array = malloc(sizeof(Expr*) * dim_count);
    
    if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->type == EXPR_SYMBOL && n_spec->data.function.head->data.symbol.name == SYM_List) {
        for(size_t i=0; i<dim_count; i++) n_array[i] = n_spec->data.function.args[i];
    } else {
        n_array[0] = n_spec;
    }
    
    if (r_spec) {
        if (r_spec->type == EXPR_FUNCTION && r_spec->data.function.head->type == EXPR_SYMBOL && r_spec->data.function.head->data.symbol.name == SYM_List && r_spec->data.function.arg_count == dim_count) {
            for(size_t i=0; i<dim_count; i++) r_array[i] = r_spec->data.function.args[i];
        } else {
            for(size_t i=0; i<dim_count; i++) r_array[i] = r_spec;
        }
    } else {
        for(size_t i=0; i<dim_count; i++) r_array[i] = NULL;
    }
    
    for (size_t i = 0; i < dim_count; i++) {
        if (n_array[i]->type != EXPR_INTEGER || n_array[i]->data.integer < 0) {
            free(n_array);
            free(r_array);
            return NULL;
        }
    }

    /* The dimensions are the user's: level i allocates d0*...*di result
     * slots in all, so each prefix product must be allocatable.  One past
     * memory (or overflowing) aborts with General::nomem instead of an
     * unchecked malloc written through (audit Array#1, Array[f, 10^15]). */
    {
        int64_t cells = 1;
        for (size_t i = 0; i < dim_count; i++) {
            if (!mth_count_mul(cells, n_array[i]->data.integer, &cells) ||
                !mth_alloc_ok(cells, sizeof(Expr*), NULL)) {
                free(n_array);
                free(r_array);
                return mth_nomem_abort();
            }
        }
    }

    Expr** current_args = malloc(sizeof(Expr*) * dim_count);
    Expr* result = array_helper(f, n_array, r_array, dim_count, 0, current_args);
    free(current_args);
    
    free(n_array);
    free(r_array);
    
    /* Offered rather than direct-built: array_helper runs evaluate(f[i]) per leaf
     * with no advance knowledge of the result type, so a buffer opened up front
     * would be abandoned in the common (symbolic) case. */
    return g_array_no_pack ? result : pack_offer(result);
}
