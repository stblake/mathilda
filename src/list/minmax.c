#include "list_common.h"
#include "minmax.h"
#include "assoc.h"

/* MinMax[list] gives {Min[list], Max[list]} in the natural single pass a caller
 * would otherwise write by hand. Over an association it uses the values (Min and
 * Max already thread there). Delegating to Min / Max keeps the numeric handling
 * — bignums, reals, symbolic extrema, empty-list Infinity — in exactly one place. */
Expr* builtin_minmax(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    bool is_list = arg->type == EXPR_FUNCTION &&
                   arg->data.function.head->type == EXPR_SYMBOL &&
                   arg->data.function.head->data.symbol.name == SYM_List;
    if (!is_list && !is_association(arg)) return NULL;

    Expr* min_arg[1] = { expr_copy(arg) };
    Expr* min_call = expr_new_function(expr_new_symbol(SYM_Min), min_arg, 1);
    Expr* max_arg[1] = { expr_copy(arg) };
    Expr* max_call = expr_new_function(expr_new_symbol(SYM_Max), max_arg, 1);
    Expr* mn = evaluate(min_call); expr_free(min_call);
    Expr* mx = evaluate(max_call); expr_free(max_call);

    Expr* items[2] = { mn, mx };
    return expr_new_function(expr_new_symbol(SYM_List), items, 2);
}

Expr* builtin_min(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_Infinity);

    /* Min[assoc] is the minimum of the association's values. */
    if (n == 1 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
    }

    // Check for List arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_List) {
            has_list = true;
            break;
        }
    }

    if (has_list) {
        size_t new_count = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
                arg->data.function.head->data.symbol.name == SYM_List) {
                new_count += arg->data.function.arg_count;
            } else {
                new_count++;
            }
        }

        Expr** new_args = malloc(sizeof(Expr*) * new_count);
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
                arg->data.function.head->data.symbol.name == SYM_List) {
                for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                    new_args[k++] = expr_copy(arg->data.function.args[j]);
                }
            } else {
                new_args[k++] = expr_copy(arg);
            }
        }
        Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
        free(new_args);
        return ret;
    }

    /* Min of a single element is that element (Min[5] -> 5, Min[x] -> x). Lists
     * were flattened above, so a lone argument here is a genuine scalar. */
    if (n == 1) return expr_copy(res->data.function.args[0]);

    // Check for Overflow[] and -Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        if (is_minus_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* min_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_real_numeric(arg)) {
            if (min_num == NULL) {
                min_num = expr_copy(arg);
            } else {
                if (expr_compare(arg, min_num) < 0) {
                    expr_free(min_num);
                    min_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (min_num ? 1 : 0);
        if (final_count == 0) {
            if (min_num) expr_free(min_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return expr_new_symbol(SYM_Infinity);
        }
        if (final_count == 1) {
            Expr* single = min_num ? min_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (min_num) final_args[k++] = min_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (min_num) expr_free(min_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}

Expr* builtin_max(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return make_minus_infinity();

    /* Max[assoc] is the maximum of the association's values. */
    if (n == 1 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
    }

    // Check for List arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
            arg->data.function.head->data.symbol.name == SYM_List) {
            has_list = true;
            break;
        }
    }
    
    if (has_list) {
        size_t new_count = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol.name == SYM_List) {
                new_count += arg->data.function.arg_count;
            } else {
                new_count++;
            }
        }
        
        Expr** new_args = malloc(sizeof(Expr*) * new_count);
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol.name == SYM_List) {
                for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                    new_args[k++] = expr_copy(arg->data.function.args[j]);
                }
            } else {
                new_args[k++] = expr_copy(arg);
            }
        }
        Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
        free(new_args);
        return ret;
    }

    /* Max of a single element is that element (Max[5] -> 5, Max[x] -> x). Lists
     * were flattened above, so a lone argument here is a genuine scalar. */
    if (n == 1) return expr_copy(res->data.function.args[0]);

    // Check for Overflow[] and Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        if (is_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* max_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_real_numeric(arg)) {
            if (max_num == NULL) {
                max_num = expr_copy(arg);
            } else {
                if (expr_compare(arg, max_num) > 0) {
                    expr_free(max_num);
                    max_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_minus_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (max_num ? 1 : 0);
        if (final_count == 0) {
            if (max_num) expr_free(max_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return make_minus_infinity();
        }
        if (final_count == 1) {
            Expr* single = max_num ? max_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (max_num) final_args[k++] = max_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (max_num) expr_free(max_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}
