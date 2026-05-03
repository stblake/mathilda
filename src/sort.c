#include "sort.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

static Expr* current_sort_p = NULL;

static int custom_compare(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    
    if (current_sort_p == NULL) {
        return expr_compare(ea, eb);
    }
    
    // Call p[ea, eb]
    Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
    Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    
    int cmp = -1; // Default to "in order" (e1 < e2)
    
    if (result->type == EXPR_INTEGER) {
        // 1: e1 before e2 (e1 < e2)
        // 0: same
        // -1: e1 after e2 (e1 > e2)
        if (result->data.integer == 1) cmp = -1;
        else if (result->data.integer == 0) cmp = 0;
        else if (result->data.integer == -1) cmp = 1;
    } else if (result->type == EXPR_SYMBOL) {
        if (result->data.symbol == SYM_True) cmp = -1;
        else if (result->data.symbol == SYM_False) cmp = 1;
    }
    
    expr_free(result);
    return cmp;
}

Expr* builtin_sort(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) {
        return expr_copy(list);
    }
    
    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    size_t count = list->data.function.arg_count;
    if (count == 0) {
        return expr_copy(list);
    }
    
    Expr** sorted_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        sorted_args[i] = expr_copy(list->data.function.args[i]);
    }
    
    // Set global context for comparison
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    qsort(sorted_args, count, sizeof(Expr*), custom_compare);
    
    current_sort_p = old_p;
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), sorted_args, count);
    free(sorted_args);
    
    return result;
}


Expr* builtin_orderedq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) {
        return expr_new_symbol("True");
    }
    
    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    size_t count = list->data.function.arg_count;
    
    if (count <= 1) {
        return expr_new_symbol("True");
    }
    
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    bool ordered = true;
    for (size_t i = 0; i < count - 1; i++) {
        Expr* ea = list->data.function.args[i];
        Expr* eb = list->data.function.args[i+1];
        
        int cmp = 0;
        if (current_sort_p == NULL) {
            cmp = expr_compare(ea, eb);
        } else {
            Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
            Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
            Expr* result = evaluate(call);
            expr_free(call);
            
            cmp = -1; // Default to "in order" (ea <= eb)
            if (result->type == EXPR_INTEGER) {
                if (result->data.integer == 1) cmp = -1;
                else if (result->data.integer == 0) cmp = 0;
                else if (result->data.integer == -1) cmp = 1;
            } else if (result->type == EXPR_SYMBOL) {
                if (result->data.symbol == SYM_True) cmp = -1;
                else if (result->data.symbol == SYM_False) cmp = 1;
            }
            expr_free(result);
        }
        
        if (cmp > 0) {
            ordered = false;
            break;
        }
    }
    
    current_sort_p = old_p;
    
    return expr_new_symbol(ordered ? "True" : "False");
}
