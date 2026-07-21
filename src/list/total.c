#include "list_common.h"
#include "total.h"
#include "assoc.h"
#include "ndreduce.h"

static int64_t get_depth_for_total(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_depth_for_total(e->data.function.args[i]);
        if (d > max_d) max_d = d;
    }
    return 1 + max_d;
}

static Expr* total_at_exactly_level_k(Expr* e, int64_t k) {
    if (k <= 0) return expr_copy(e);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (k == 1) {
        size_t count = e->data.function.arg_count;
        if (count == 0) return expr_new_integer(0);
        Expr** plus_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) plus_args[i] = expr_copy(e->data.function.args[i]);
        Expr* plus_expr = expr_new_function(expr_new_symbol(SYM_Plus), plus_args, count);
        free(plus_args);
        Expr* res = evaluate(plus_expr);
        expr_free(plus_expr);
        return res;
    } else {
        size_t count = e->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = total_at_exactly_level_k(e->data.function.args[i], k - 1);
        }
        Expr* res = expr_new_function(expr_copy(e->data.function.head), new_args, count);
        free(new_args);
        return res;
    }
}

Expr* builtin_total(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;

    /* Total[assoc] sums the association's values. */
    { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path: sum the flat buffer directly (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_total(res);

    Expr* list = res->data.function.args[0];
    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    int64_t n1 = 1, n2 = 1;
    int64_t depth = get_depth_for_total(list);

    if (level_spec) {
        if (level_spec->type == EXPR_INTEGER) {
            n1 = 1;
            n2 = level_spec->data.integer;
            if (n2 < 0) n2 = depth + n2;
        } else if (is_infinity(level_spec)) {
            n1 = 1;
            n2 = depth - 1;
        } else if (is_listq(level_spec)) {
            if (level_spec->data.function.arg_count == 1) {
                Expr* arg = level_spec->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    n1 = n2 = arg->data.integer;
                    if (n1 < 0) n1 = n2 = depth + n1;
                } else if (is_infinity(arg)) {
                    n1 = n2 = depth - 1;
                } else return NULL;
            } else if (level_spec->data.function.arg_count == 2) {
                Expr* arg1 = level_spec->data.function.args[0];
                Expr* arg2 = level_spec->data.function.args[1];
                if (arg1->type == EXPR_INTEGER) {
                    n1 = arg1->data.integer;
                    if (n1 < 0) n1 = depth + n1;
                } else return NULL;
                
                if (arg2->type == EXPR_INTEGER) {
                    n2 = arg2->data.integer;
                    if (n2 < 0) n2 = depth + n2;
                } else if (is_infinity(arg2)) {
                    n2 = depth - 1;
                } else return NULL;
            } else return NULL;
        } else return NULL;
    }
    
    if (n1 > n2) return expr_copy(list);
    if (n1 < 1) n1 = 1;
    if (n2 >= depth) n2 = depth - 1;

    Expr* current = expr_copy(list);
    for (int64_t k = n2; k >= n1; k--) {
        Expr* next = total_at_exactly_level_k(current, k);
        expr_free(current);
        current = next;
    }

    return current;
}
