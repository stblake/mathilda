/* stats_common.c -- helpers shared across the src/stats/ per-builtin files.
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats_common.h"
#include "arithmetic.h"   /* is_rational */
#include "complex.h"      /* is_complex */
#include "eval.h"         /* evaluate */
#include "sym_names.h"    /* SYM_* */

bool stats_is_numeric(Expr* e, double* val, bool* out_complex) {
    if (e->type == EXPR_INTEGER) {
        if (val) *val = (double)e->data.integer;
        if (out_complex) *out_complex = false;
        return true;
    }
    if (e->type == EXPR_REAL) {
        if (val) *val = e->data.real;
        if (out_complex) *out_complex = false;
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        if (val) *val = (double)n / (double)d;
        if (out_complex) *out_complex = false;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (out_complex) *out_complex = true;
        return true;
    }
    return false;
}

Expr* stats_apply_columnwise(const char* func_name, Expr* matrix) {
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

bool stats_is_real_numeric(Expr* e) {
    Expr* numq = expr_new_function(expr_new_symbol(SYM_NumericQ), (Expr*[]){expr_copy(e)}, 1);
    Expr* numq_eval = evaluate(numq);
    expr_free(numq);
    if (numq_eval->type != EXPR_SYMBOL || numq_eval->data.symbol.name != SYM_True) {
        expr_free(numq_eval);
        return false;
    }
    expr_free(numq_eval);

    Expr* freeq = expr_new_function(expr_new_symbol(SYM_FreeQ), (Expr*[]){expr_copy(e), expr_new_symbol(SYM_I)}, 2);
    Expr* freeq_eval = evaluate(freeq);
    expr_free(freeq);
    if (freeq_eval->type != EXPR_SYMBOL || freeq_eval->data.symbol.name != SYM_True) {
        expr_free(freeq_eval);
        return false;
    }
    expr_free(freeq_eval);

    return true;
}
