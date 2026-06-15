/* Shared linear-algebra helpers used across the src/linalg/ module and
 * exposed (via linalg.h) to matsol.c, matinv.c, matlstsq.c, mateigen.c.
 *
 * Contains:
 *   get_tensor_dims    -- rank / shape probe for nested-List tensors
 *   flatten_tensor     -- row-major flatten into a pre-allocated array
 *   exact_div_wrapper  -- polynomial exact-division wrapper used by
 *                         RowReduce / LinearSolve to keep intermediate
 *                         entries simplified.
 */

#include "linalg.h"
#include "eval.h"
#include "expand.h"
#include "poly.h"
#include "sym_names.h"
#include <stdlib.h>

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
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){exp_num, t}, 2));
        return r;
    }
}
