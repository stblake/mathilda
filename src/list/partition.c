#include "list_common.h"
#include "partition.h"

static Expr* partition_rec(Expr* list, Expr* n_spec, Expr* d_spec, size_t level_idx) {
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    // Get n and d for this level
    int64_t n = -1;
    bool n_upto = false;
    Expr* n_e = (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol == SYM_List) ? 
                (level_idx < n_spec->data.function.arg_count ? n_spec->data.function.args[level_idx] : NULL) : 
                (level_idx == 0 ? n_spec : NULL);
    
    if (!n_e) return expr_copy(list);

    if (n_e->type == EXPR_INTEGER) {
        n = n_e->data.integer;
    } else if (n_e->type == EXPR_FUNCTION && n_e->data.function.head->data.symbol == SYM_UpTo && n_e->data.function.arg_count == 1) {
        if (n_e->data.function.args[0]->type == EXPR_INTEGER) {
            n = n_e->data.function.args[0]->data.integer;
            n_upto = true;
        }
    }
    if (n <= 0) return expr_copy(list);

    int64_t d = n;
    if (d_spec) {
        Expr* d_e = (d_spec->type == EXPR_FUNCTION && d_spec->data.function.head->data.symbol == SYM_List) ? 
                    (level_idx < d_spec->data.function.arg_count ? d_spec->data.function.args[level_idx] : NULL) : 
                    (level_idx == 0 ? d_spec : NULL);
        if (d_e && d_e->type == EXPR_INTEGER) {
            d = d_e->data.integer;
        }
    }
    if (d <= 0) return expr_copy(list);

    size_t len = list->data.function.arg_count;
    size_t num_sublists = 0;
    if (n_upto) {
        num_sublists = (len > 0) ? (len - 1) / d + 1 : 0;
    } else {
        if (len >= (size_t)n) {
            num_sublists = (len - (size_t)n) / (size_t)d + 1;
        }
    }

    Expr** sublists = malloc(sizeof(Expr*) * num_sublists);
    for (size_t i = 0; i < num_sublists; i++) {
        size_t start = i * (size_t)d;
        size_t end = start + (size_t)n;
        if (end > len) end = len;
        
        size_t sub_count = end - start;
        Expr** sub_args = malloc(sizeof(Expr*) * sub_count);
        for (size_t j = 0; j < sub_count; j++) {
            sub_args[j] = partition_rec(list->data.function.args[start + j], n_spec, d_spec, level_idx + 1);
        }
        sublists[i] = expr_new_function(expr_copy(list->data.function.head), sub_args, sub_count);
        free(sub_args);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), sublists, num_sublists);
    free(sublists);
    return result;
}

Expr* builtin_partition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* n_spec = res->data.function.args[1];
    Expr* d_spec = (res->data.function.arg_count == 3) ? res->data.function.args[2] : NULL;

    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    return partition_rec(list, n_spec, d_spec, 0);
}
