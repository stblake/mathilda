/* quantile.c -- Quantile[].
 *
 * Quantile[data, q]                 -- Wolfram default parameters {{0,0},{1,0}}:
 *                                      left-continuous, sorted[[Ceiling[n q]]].
 * Quantile[data, {q1, q2, ...}]     -- list of quantiles, one result per q.
 * Quantile[data, q, {{a,b},{c,d}}]  -- the general parameterized definition
 *                                      (same form Quartiles accepts; Quartiles
 *                                      is this with q = {1/4,1/2,3/4} and
 *                                      parameters {{1/2,0},{0,1}}).
 *
 * The interpolation engine is shared with Quartiles: stats_quantile_point in
 * stats_common.c. Matrix input recurses columnwise; a visible NDArray argument
 * is materialised to the exact List path via pack_unpack (correctness-first --
 * no ndreduce kernel yet; the audit baselines carry the declared reason).
 * See stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "ndreduce.h"
#include "pack.h"      /* pack_eval_plain, pack_unpack */
#include "print.h"     /* expr_to_string */
#include <stdio.h>
#include <stdlib.h>

/* Parse an optional {{a,b},{c,d}} parameter matrix (same shape Quartiles
 * accepts). On success fills the four out-slots with fresh copies and returns
 * true; on a malformed spec returns false with nothing allocated. */
static bool parse_param_matrix(Expr* param_expr, Expr** a, Expr** b, Expr** c, Expr** d) {
    if (param_expr->type != EXPR_FUNCTION || param_expr->data.function.arg_count != 2) return false;
    Expr* row1 = param_expr->data.function.args[0];
    Expr* row2 = param_expr->data.function.args[1];
    if (row1->type != EXPR_FUNCTION || row1->data.function.arg_count != 2 ||
        row2->type != EXPR_FUNCTION || row2->data.function.arg_count != 2) return false;
    *a = expr_copy(row1->data.function.args[0]);
    *b = expr_copy(row1->data.function.args[1]);
    *c = expr_copy(row2->data.function.args[0]);
    *d = expr_copy(row2->data.function.args[1]);
    return true;
}

/* One validated q in [0,1] -> engine result. Returns NULL when q is not real
 * numeric (caller declines) and sets *out_of_range when q is numeric but
 * outside [0,1] (caller emits the message). */
static Expr* quantile_one(Expr** sorted_A, size_t n, Expr* q,
                          Expr* a, Expr* b, Expr* c, Expr* d, bool* out_of_range) {
    *out_of_range = false;
    if (!stats_is_real_numeric(q)) return NULL;
    double qv = 0;
    if (!stats_is_numeric(q, &qv, NULL)) return NULL;
    if (qv < 0.0 || qv > 1.0) { *out_of_range = true; return NULL; }
    return stats_quantile_point(sorted_A, n, q, a, b, c, d);
}

Expr* builtin_quantile(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* data = res->data.function.args[0];
    Expr* qspec = res->data.function.args[1];
    Expr* param_expr = (argc == 3) ? res->data.function.args[2] : NULL;

    /* Visible NDArray / packed data: materialise to the exact List path.
     * pack_unpack returns NULL for a non-ndarray, so this is a no-op for
     * ordinary lists. No kernel yet -- correctness first (plan STATS-1). */
    Expr* unpacked = ndred_call_has_ndarray(res) ? pack_unpack(data) : NULL;
    if (unpacked) data = unpacked;
#define QUANTILE_RET(x) do { Expr* _r = (x); if (unpacked) expr_free(unpacked); return _r; } while (0)

    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List) {
        QUANTILE_RET(expr_copy(res));
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) QUANTILE_RET(expr_copy(res));

    /* Matrix: recurse columnwise, carrying q and the parameter matrix, the way
     * Quartiles' own matrix branch does. */
    if (data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {

        Expr* t_expr = expr_new_function(expr_new_symbol(SYM_Transpose), (Expr*[]){expr_copy(data)}, 1);
        Expr* t_eval = evaluate(t_expr);
        expr_free(t_expr);

        if (t_eval->type != EXPR_FUNCTION) {
            expr_free(t_eval);
            QUANTILE_RET(expr_copy(res));
        }

        size_t cols = t_eval->data.function.arg_count;
        Expr** col_args = malloc(sizeof(Expr*) * cols);
        for (size_t i = 0; i < cols; i++) {
            Expr* call_args[3];
            call_args[0] = expr_copy(t_eval->data.function.args[i]);
            call_args[1] = expr_copy(qspec);
            size_t call_cnt = 2;
            if (param_expr) {
                call_args[2] = expr_copy(param_expr);
                call_cnt = 3;
            }
            Expr* q_call = expr_new_function(expr_new_symbol("Quantile"), call_args, call_cnt);
            col_args[i] = evaluate(q_call);
            expr_free(q_call);
        }
        Expr* res_list = expr_new_function(expr_new_symbol(SYM_List), col_args, cols);
        free(col_args);
        expr_free(t_eval);
        QUANTILE_RET(res_list);
    }

    for (size_t i = 0; i < n; i++) {
        if (!stats_is_real_numeric(data->data.function.args[i])) {
            char* str = expr_to_string(res);
            printf("Quantile::rectn: Rectangular array of real numbers is expected at position 1 in %s.\n", str);
            free(str);
            QUANTILE_RET(expr_copy(res));
        }
    }

    Expr *a = NULL, *b = NULL, *c = NULL, *d = NULL;
    if (param_expr) {
        if (!parse_param_matrix(param_expr, &a, &b, &c, &d)) QUANTILE_RET(expr_copy(res));
    } else {
        /* Wolfram's Quantile defaults {{0,0},{1,0}}: h = n q, result
         * sorted[[Ceiling[h]]] (edge-clamped). */
        a = expr_new_integer(0);
        b = expr_new_integer(0);
        c = expr_new_integer(1);
        d = expr_new_integer(0);
    }
#define QUANTILE_RET_ABCD(x) do { Expr* _r2 = (x); expr_free(a); expr_free(b); expr_free(c); expr_free(d); QUANTILE_RET(_r2); } while (0)

    /* Sort once. pack_eval_plain, not evaluate: Sort returns a PACKED list for
     * a large machine-number input and the args walk below reads
     * .data.function.args (same trap note as median.c / quartiles.c). */
    Expr* sort_expr = expr_new_function(expr_new_symbol(SYM_Sort), (Expr*[]){expr_copy(data)}, 1);
    Expr* sorted = pack_eval_plain(sort_expr);
    expr_free(sort_expr);

    if (sorted->type != EXPR_FUNCTION || sorted->data.function.arg_count != n) {
        expr_free(sorted);
        QUANTILE_RET_ABCD(expr_copy(res));
    }
    Expr** sorted_A = sorted->data.function.args;

    bool qspec_is_list = qspec->type == EXPR_FUNCTION &&
                         qspec->data.function.head->type == EXPR_SYMBOL &&
                         qspec->data.function.head->data.symbol.name == SYM_List;

    if (!qspec_is_list) {
        bool out_of_range = false;
        Expr* r = quantile_one(sorted_A, n, qspec, a, b, c, d, &out_of_range);
        if (!r) {
            if (out_of_range) {
                char* str = expr_to_string(res);
                printf("Quantile::q100: The quantile q is expected to be a number between 0 and 1 in %s.\n", str);
                free(str);
            }
            expr_free(sorted);
            QUANTILE_RET_ABCD(expr_copy(res));
        }
        expr_free(sorted);
        QUANTILE_RET_ABCD(r);
    }

    size_t nq = qspec->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (nq ? nq : 1));
    for (size_t i = 0; i < nq; i++) {
        bool out_of_range = false;
        Expr* r = quantile_one(sorted_A, n, qspec->data.function.args[i], a, b, c, d, &out_of_range);
        if (!r) {
            if (out_of_range) {
                char* str = expr_to_string(res);
                printf("Quantile::q100: The quantile q is expected to be a number between 0 and 1 in %s.\n", str);
                free(str);
            }
            for (size_t k = 0; k < i; k++) expr_free(out[k]);
            free(out);
            expr_free(sorted);
            QUANTILE_RET_ABCD(expr_copy(res));
        }
        out[i] = r;
    }
    Expr* res_list = expr_new_function(expr_new_symbol(SYM_List), out, nq);
    free(out);
    expr_free(sorted);
    QUANTILE_RET_ABCD(res_list);
}
