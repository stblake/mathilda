#include "list_common.h"
#include "constant_array.h"
#include "../ndarray.h"
#include "../pack.h"

/* ConstantArray[c, n]            -> a flat list of n copies of c.
 * ConstantArray[c, {n1, ..., nk}] -> a nested n1 x ... x nk array of copies of c.
 *
 * Structurally this is Array[] minus the index computation: at each leaf we
 * simply copy the constant element instead of building and evaluating an
 * indexed function call, and there is no index-range machinery. See
 * src/list/array.c for the analog it mirrors. */

/* Build the nested List at recursion depth `current_dim`. The deepest level
 * (current_dim == dim_count) yields a fresh copy of the constant element c.
 * `n_array` holds borrowed pointers to each dimension's EXPR_INTEGER size,
 * already validated non-negative by the caller. */
static Expr* ca_helper(Expr* c, Expr** n_array, size_t dim_count, size_t current_dim) {
    if (current_dim == dim_count) {
        return expr_copy(c);
    }

    int64_t n_val = n_array[current_dim]->data.integer;
    if (n_val == 0) {
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    Expr** results = malloc(sizeof(Expr*) * (size_t)n_val);
    for (int64_t i = 0; i < n_val; i++) {
        results[i] = ca_helper(c, n_array, dim_count, current_dim + 1);
    }

    Expr* list_result = expr_new_function(expr_new_symbol(SYM_List), results, (size_t)n_val);
    free(results);
    return list_result;
}

Expr* builtin_constant_array(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) {
        return builtin_arg_error("ConstantArray", argc, 2, 3);
    }
    /* The optional 3rd (padding) argument is accepted for arity purposes but
     * not yet implemented; leave such calls unevaluated. */
    if (argc != 2) return NULL;

    Expr* c = res->data.function.args[0];
    Expr* dims = res->data.function.args[1];

    /* A dimension spec is either a single integer or a List of integers. */
    size_t dim_count = 1;
    int is_list = (dims->type == EXPR_FUNCTION
                   && dims->data.function.head->type == EXPR_SYMBOL
                   && dims->data.function.head->data.symbol.name == SYM_List);
    if (is_list) {
        dim_count = dims->data.function.arg_count;
        if (dim_count == 0) return NULL; /* ConstantArray[c, {}] stays unevaluated */
    }

    Expr** n_array = malloc(sizeof(Expr*) * dim_count);
    if (is_list) {
        for (size_t i = 0; i < dim_count; i++) n_array[i] = dims->data.function.args[i];
    } else {
        n_array[0] = dims;
    }

    /* Every dimension must be a non-negative machine integer; otherwise leave
     * the call unevaluated (mirrors Array[]). */
    for (size_t i = 0; i < dim_count; i++) {
        if (n_array[i]->type != EXPR_INTEGER || n_array[i]->data.integer < 0) {
            free(n_array);
            return NULL;
        }
    }

    /* A machine-number constant over a rectangular shape is the one case where
     * the whole result is known before anything is built, so fill the buffer
     * directly rather than allocating one Expr per element and packing after.
     * ndbuild_open declines a zero dimension, which routes the ConstantArray[c,
     * {0, ...}] shapes back to ca_helper unchanged. */
    if ((c->type == EXPR_REAL || c->type == EXPR_INTEGER) &&
        dim_count >= 1 && dim_count < NDARRAY_MAX_RANK) {
        int64_t d[NDARRAY_MAX_RANK];
        int64_t total = 1;
        for (size_t i = 0; i < dim_count; i++) {
            d[i] = n_array[i]->data.integer;
            total *= d[i];
        }
        NDType dt = (c->type == EXPR_INTEGER) ? NDT_INT64 : NDT_FLOAT64;
        void* vbuf = NULL;
        Expr* packed = ndbuild_open((int)dim_count, d, dt, &vbuf);
        if (packed) {
            if (dt == NDT_INT64) {
                int64_t v = c->data.integer, *buf = (int64_t*)vbuf;
                for (int64_t i = 0; i < total; i++) buf[i] = v;
            } else {
                double v = c->data.real, *buf = (double*)vbuf;
                for (int64_t i = 0; i < total; i++) buf[i] = v;
            }
            free(n_array);
            return packed;
        }
    }

    Expr* result = ca_helper(c, n_array, dim_count, 0);
    free(n_array);
    return result;
}
