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
#include <stdio.h>
#include <stdlib.h>

Expr* builtin_meandeviation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    Expr* unpacked = ndred_call_has_ndarray(res) ? pack_unpack(data) : NULL;
    if (unpacked) data = unpacked;
#define MD_RET(x) do { Expr* _r = (x); if (unpacked) expr_free(unpacked); return _r; } while (0)

    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List) {
        MD_RET(expr_copy(res));
    }
    size_t n = data->data.function.arg_count;
    if (n == 0) MD_RET(expr_copy(res));

    if (data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        Expr* r = stats_apply_columnwise("MeanDeviation", data);
        MD_RET(r ? r : expr_copy(res));
    }

    for (size_t i = 0; i < n; i++) {
        if (!stats_is_real_numeric(data->data.function.args[i])) {
            char* str = expr_to_string(res);
            printf("MeanDeviation::rectn: Rectangular array of real numbers is expected at position 1 in %s.\n", str);
            free(str);
            MD_RET(expr_copy(res));
        }
    }

    /* Mean[Abs[data - Mean[data]]], one tree, one evaluation. All-numeric input
     * was just verified, so every layer reduces. */
    Expr* inner_mean = expr_new_function(expr_new_symbol(SYM_Mean), (Expr*[]){expr_copy(data)}, 1);
    Expr* neg_mean = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), inner_mean}, 2);
    Expr* centered = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(data), neg_mean}, 2);
    Expr* absdev = expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){centered}, 1);
    Expr* outer = expr_new_function(expr_new_symbol(SYM_Mean), (Expr*[]){absdev}, 1);
    MD_RET(eval_and_free(outer));
}
