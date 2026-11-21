/* movingaverage.c -- MovingAverage[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "eval.h"
#include "sym_names.h"
#include "ndarray.h"
#include "ndreduce.h"
#include "stats_common.h"
#include "linkedlist.h" /* LLElems */
#include <stdlib.h>
#include <gmp.h>

/*
 * MovingAverage[list, r]                — simple moving average over runs of r elements.
 * MovingAverage[list, {w1, ..., wr}]    — weighted moving average with weights w_k
 *                                         (output uses w_k / Sum[w_k] as effective weights).
 *
 * Output length is Length[list] - r + 1 when 1 <= r <= Length[list]; the function
 * stays unevaluated otherwise.  The unweighted form delegates to Mean so it inherits
 * Mean's exact rational / bigint / real / symbolic handling.  The weighted form builds
 * Plus[Times[w_k/wsum, x_{i+k}], ...] and lets the evaluator simplify.
 */
/* The window el[i, i + r) as a sequence of the data's kind: a List of the
 * elements, or a LinkedList -- a sub-sequence of a LinkedList is a
 * LinkedList, its slots copied (typed stays typed). */
Expr* stats_moving_window(const Expr* data, Expr** el, size_t i, size_t r) {
    const LinkedList* ll = expr_linkedlist(data);
    if (ll) {
        LinkedList* wl = ll_new();
        ll_append_range(wl, ll, (int64_t)i, (int64_t)(i + r));
        return expr_new_linkedlist(wl);
    }
    Expr** sub = malloc(sizeof(Expr*) * r);
    for (size_t k = 0; k < r; k++) sub[k] = expr_copy(el[i + k]);
    Expr* sublist = expr_new_function(expr_new_symbol(SYM_List), sub, r);
    free(sub);
    return sublist;
}

/* MovingAverage over the n elements el (of `data`, a List or LinkedList) with
 * window spec `spec`; a weight sequence is read into *w (the caller closes
 * it). */
static Expr* moving_average_vec(Expr* data, Expr** el, size_t n, Expr* spec, LLElems* w) {
    size_t r = 0;
    Expr** weights = NULL; /* borrowed; non-NULL means weighted form */

    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 1) return NULL;
        if ((uint64_t)spec->data.integer > (uint64_t)n) return NULL;
        r = (size_t)spec->data.integer;
    } else if (spec->type == EXPR_BIGINT) {
        if (mpz_sgn(spec->data.bigint) <= 0) return NULL;
        if (!mpz_fits_ulong_p(spec->data.bigint)) return NULL;
        unsigned long rr = mpz_get_ui(spec->data.bigint);
        if ((size_t)rr > n) return NULL;
        r = (size_t)rr;
    } else if (ll_elems_open(spec, w)) {
        r = w->n;
        if (r == 0 || r > n) return NULL;
        weights = w->v;
    } else {
        return NULL;
    }

    size_t out_n = n - r + 1;
    Expr** out = malloc(sizeof(Expr*) * out_n);
    if (!out) return NULL;

    if (weights == NULL) {
        for (size_t i = 0; i < out_n; i++) {
            Expr* sublist = stats_moving_window(data, el, i, r);
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Mean), (Expr*[]){ sublist }, 1));
        }
    } else {
        Expr** w_copies = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) w_copies[k] = expr_copy(weights[k]);
        Expr* wsum = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus), w_copies, r));
        free(w_copies);

        Expr* wsum_inv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power),
            (Expr*[]){ wsum, expr_new_integer(-1) }, 2));

        Expr** coefs = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) {
            coefs[k] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(weights[k]), expr_copy(wsum_inv) }, 2));
        }
        expr_free(wsum_inv);

        for (size_t i = 0; i < out_n; i++) {
            Expr** terms = malloc(sizeof(Expr*) * r);
            for (size_t k = 0; k < r; k++) {
                terms[k] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(coefs[k]),
                               expr_copy(el[i + k]) }, 2));
            }
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Plus), terms, r));
            free(terms);
        }

        for (size_t k = 0; k < r; k++) expr_free(coefs[k]);
        free(coefs);
    }

    Expr* result = ll_seq_like(data, out, out_n);
    free(out);
    return result;
}

Expr* builtin_moving_average(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* NDArray fast path: sliding-window mean on the buffer (see ndreduce.c). */
    if (is_ndarray(data)) return ndred_moving_average(res);

    /* LinkedList: the packed kernel on the chunks (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_MOVING_AVERAGE, &r)) return r; }

    /* The data and a weight list may each be a List or a LinkedList; the
     * answer is in the image of the data. */
    LLElems d, w = { NULL, 0, false };
    if (!ll_elems_open(data, &d)) return NULL;
    Expr* out_seq = moving_average_vec(data, d.v, d.n, spec, &w);
    ll_elems_close(&d);
    ll_elems_close(&w);
    return out_seq;
}
