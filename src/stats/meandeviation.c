/* meandeviation.c -- MeanDeviation[].
 *
 * MeanDeviation[data] = Mean[Abs[data - Mean[data]]], composed through the
 * evaluator so exact input stays exact (the src/stats exactness discipline --
 * see mean.c's overflow note). Matrix input recurses columnwise; NDArray input
 * is materialised via pack_unpack (correctness-first, no kernel yet).
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


/* Infinity / ComplexInfinity / Indeterminate are NumericQ and free of I, so
 * they pass the real-numeric element gate above -- and then the composed
 * Mean/Median tree cannot reduce, handing back a half-evaluated expression.
 * The plan's "center failed to reduce -> decline" rule, enforced. */
static bool has_nonfinite(Expr* e) {
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        return n == SYM_Infinity || n == SYM_ComplexInfinity || n == SYM_Indeterminate;
    }
    if (e->type == EXPR_FUNCTION) {
        if (has_nonfinite(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_nonfinite(e->data.function.args[i])) return true;
    }
    return false;
}

Expr* builtin_meandeviation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    Expr* unpacked = ndred_call_has_ndarray(res) ? pack_unpack(data) : NULL;
    if (unpacked) data = unpacked;
#define MD_RET(x) do { Expr* _r = (x); if (unpacked) expr_free(unpacked); return _r; } while (0)

    /* LinkedList: columnwise (stats_ll.c), or the vector below. */
    { Expr* r; if (stats_ll_arm(res, STATS_MEAN_DEVIATION, &r)) MD_RET(r ? r : expr_copy(res)); }

    if (!ll_is_seq(data)) MD_RET(expr_copy(res));
    if (ll_seq_length(data) == 0) MD_RET(expr_copy(res));

    if (data->type != EXPR_LINKEDLIST &&
        data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        Expr* r = stats_apply_columnwise("MeanDeviation", data);
        MD_RET(r ? r : expr_copy(res));
    }

    LLReader e;
    ll_rd_open(data, &e);
    size_t n = e.n;
    for (size_t i = 0; i < n; i++) {
        if (!stats_is_real_numeric(ll_rd_at(&e, i)) || has_nonfinite(ll_rd_at(&e, i))) {
            ll_rd_close(&e);
            char* str = expr_to_string(res);
            mth_message("MeanDeviation", "rectn",
                        "Rectangular array of real numbers is expected at position 1 in %s.", str);
            free(str);
            MD_RET(expr_copy(res));
        }
    }

    /* Mean[Abs[data - Mean[data]]] over the elements: the vector Mean core
     * (stats_mean_vec), the centring and Abs one element at a time -- what the
     * Listable Plus and Abs do over the List, so the List and the LinkedList
     * (and an unpacked packed array) agree to the bit. All-numeric input was
     * just verified, so every layer reduces. */
    Expr* mu = stats_mean_vec(&e);
    Expr* neg_mu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), mu }, 2));
    Expr** dev = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(ll_rd_at(&e, i)), expr_copy(neg_mu) }, 2));
        dev[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){ c }, 1));
    }
    expr_free(neg_mu);
    ll_rd_close(&e);
    Expr* r = stats_mean_of(dev, n);
    for (size_t i = 0; i < n; i++) expr_free(dev[i]);
    free(dev);
    MD_RET(r ? r : expr_copy(res));
}
