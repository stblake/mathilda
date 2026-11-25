/* interquartilerange.c -- InterquartileRange[].
 *
 * InterquartileRange[data] = q3 - q1 of Quartiles[data] (Wolfram's IQR uses the
 * Quartiles parameterization {{1/2,0},{0,1}}, so composing over the Quartiles
 * builtin is the definition, not a shortcut). Matrix input recurses per-column
 * BEFORE the vector path: a k-column matrix's Quartiles result is a k-list of
 * triples, and for k == 3 a bare "is it a 3-list" guard would confuse it with a
 * vector's {q1,q2,q3} and compute Quartiles(col3)-Quartiles(col1) -- the silent
 * wrong answer the STATS-1 plan review flagged. The vector-path guard also
 * requires all three quartiles to be SCALARS for the same reason.
 * NDArray input is materialised via pack_unpack (correctness-first, no kernel).
 * See stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "sym_names.h"
#include "ndreduce.h"
#include "pack.h"      /* pack_unpack */
#include "print.h"     /* expr_to_string */
#include "linkedlist.h" /* LLReader */
#include "message.h"
#include <stdio.h>
#include <stdlib.h>

Expr* builtin_interquartilerange(Expr* res) {
    /* InterquartileRange[data] and InterquartileRange[data, {{a, b}, {c, d}}]:
     * q3 - q1 of Quartiles[data(, params)]. */
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 ||
        res->data.function.arg_count > 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* params = res->data.function.arg_count == 2 ? res->data.function.args[1] : NULL;

    Expr* unpacked = ndred_call_has_ndarray(res) ? pack_unpack(data) : NULL;
    if (unpacked) data = unpacked;
#define IQR_RET(x) do { Expr* _r = (x); if (unpacked) expr_free(unpacked); return _r; } while (0)

    /* LinkedList: columnwise (stats_ll.c), or the vector below. */
    { Expr* r; if (!params && stats_ll_arm(res, STATS_IQR, &r)) IQR_RET(r ? r : expr_copy(res)); }

    if (!ll_is_seq(data)) IQR_RET(expr_copy(res));
    if (ll_seq_length(data) == 0) IQR_RET(expr_copy(res));

    if (data->type != EXPR_LINKEDLIST &&
        data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        if (params) IQR_RET(NULL);   /* columnwise with parameters: not handled */
        Expr* r = stats_apply_columnwise("InterquartileRange", data);
        if (!r) stats_msg_rectn("InterquartileRange", res);   /* ragged / mixed depth */
        IQR_RET(r);   /* NULL after the message */
    }

    LLReader e;
    ll_rd_open(data, &e);
    size_t n = e.n;
    for (size_t i = 0; i < n; i++) {
        if (!stats_is_real_numeric(ll_rd_at(&e, i))) {
            ll_rd_close(&e);
            char* str = expr_to_string(res);
            mth_message("InterquartileRange", "rectn", "Rectangular array of real numbers is expected at position 1 in %s.", str);
            free(str);
            IQR_RET(NULL);   /* rest after the message */
        }
    }

    ll_rd_close(&e);

    Expr* q_call = params
        ? expr_new_function(expr_new_symbol(SYM_Quartiles), (Expr*[]){expr_copy(data), expr_copy(params)}, 2)
        : expr_new_function(expr_new_symbol(SYM_Quartiles), (Expr*[]){expr_copy(data)}, 1);
    Expr* q_eval = evaluate(q_call);
    expr_free(q_call);

    bool ok = q_eval->type == EXPR_FUNCTION &&
              q_eval->data.function.head->type == EXPR_SYMBOL &&
              q_eval->data.function.head->data.symbol.name == SYM_List &&
              q_eval->data.function.arg_count == 3 &&
              stats_is_real_numeric(q_eval->data.function.args[0]) &&
              stats_is_real_numeric(q_eval->data.function.args[1]) &&
              stats_is_real_numeric(q_eval->data.function.args[2]);
    if (!ok) {
        expr_free(q_eval);
        IQR_RET(NULL);   /* Quartiles declined (and reported) */
    }

    Expr* neg_q1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(q_eval->data.function.args[0])}, 2));
    Expr* iqr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(q_eval->data.function.args[2]), neg_q1}, 2));
    expr_free(q_eval);
    IQR_RET(iqr);
}
