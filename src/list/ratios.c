#include "list_common.h"
#include "message.h"   /* mth_message */
#include "ndarray.h"    /* is_ndarray, ndarray_delist_and_reeval */
#include "ndstruct.h"   /* NDArray buffer fast paths */
#include "ratios.h"
#include "assoc.h"
#include "linkedlist_seq.h"  /* LLSeq: a LinkedList's ratios from its slots */
#include <math.h>

/* ---- Ratios ----------------------------------------------------------------
 *
 * Ratios[list]                   successive (first) ratios of list.
 * Ratios[list, n]                n-fold iterated ratios (length l - n).
 * Ratios[list, {n1, n2, ...}]    successive n_k-th ratios at level k of a
 *                                nested list / array.
 *
 * Ratios is the multiplicative analog of Differences: the element-wise
 * subtraction list[[i+1]] - list[[i]] becomes the division list[[i+1]] /
 * list[[i]], built as Divide[hi, lo] and reduced by the evaluator. Integers,
 * bignums (-> exact Rationals), machine doubles, symbolic terms, and whole
 * sublists (matrix rows, via the Listable Power/Times that Divide rewrites to)
 * all combine correctly. Because row division threads element-wise,
 * Ratios[m, n] on a matrix m takes ratios of successive rows within each
 * column, matching Ratios[m, {n, 0}]. FoldList[Times, x, Ratios[list]] inverts
 * Ratios.
 */

/* Returns the evaluated expression (numerator / denominator). The operands are
 * copied; ownership of the result transfers to the caller. */
static Expr* ratio_divide(Expr* numerator, Expr* denominator) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(numerator);
    args[1] = expr_copy(denominator);
    Expr* div = expr_new_function(expr_new_symbol(SYM_Divide), args, 2);
    free(args);
    return eval_and_free(div);
}

/* One ratio pass over the top level of lst (an EXPR_FUNCTION). The result keeps
 * lst's head; the i-th element is elem[i+1] / elem[i]. A list of length <= 1
 * yields the empty list. */
/* Ratios of a machine-real LinkedList on its slots: s[i+1] / s[i] per pair,
 * as ndstruct_ratios computes the packed List's. NULL for any other kind, a
 * zero divisor (ComplexInfinity) or a non-finite quotient: the element path. */
static Expr* ratio_once_ll_machine(const LinkedList* ll) {
    int64_t n = ll_length(ll);
    if (n < 2 || ll_kind(ll) != LL_REAL) return NULL;
    double* v = malloc(sizeof(double) * (size_t)n);
    if (!v) return NULL;
    ll_gather(ll, 0, n, v);
    bool ok = true;
    for (int64_t i = 0; i < n - 1 && ok; i++) {
        ok = v[i] != 0.0;
        if (ok) { v[i] = v[i + 1] / v[i]; ok = isfinite(v[i]); }
    }
    Expr* r = ok ? expr_new_linkedlist(ll_from_buffer(LL_REAL, v, n - 1)) : NULL;
    free(v);
    return r;
}

static Expr* ratio_once(Expr* lst) {
    if (lst->type == EXPR_LINKEDLIST) {
        Expr* fast = ratio_once_ll_machine(lst->data.llist);
        if (fast) return fast;
    }
    LLSeq lv;
    llq_open(lst, &lv);
    size_t n = lv.n;

    if (n <= 1) {
        Expr* r = llq_result(&lv, NULL, 0);
        llq_close(&lv);
        return r;
    }

    size_t outn = n - 1;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = llq_peek(&lv, i);
        Expr* hi = llq_peek(&lv, i + 1);
        out[i] = ratio_divide(hi, lo);
    }
    Expr* result = llq_result(&lv, out, outn);
    free(out);
    llq_close(&lv);
    return result;
}

/* Apply ratio_once a total of n times. Returns a fresh expression (a copy of
 * lst when n == 0). */
static Expr* ratio_n(Expr* lst, int64_t n) {
    /* Zero passes still answer a NEW LinkedList, never the argument itself. */
    Expr* cur = lst->type == EXPR_LINKEDLIST ? expr_new_linkedlist(ll_copy(lst->data.llist))
                                             : expr_copy(lst);
    for (int64_t k = 0; k < n; k++) {
        Expr* nxt = ratio_once(cur);
        expr_free(cur);
        cur = nxt;
    }
    return cur;
}

/* True iff every level spec entry is 0 (nothing left to do below here). */
static bool levels_all_zero(const int64_t* levels, size_t num) {
    for (size_t i = 0; i < num; i++) if (levels[i] != 0) return false;
    return true;
}

/* Multidimensional ratios: apply levels[0] first ratios at the top level, then
 * recurse into each element with the remaining levels.
 * A level that does not exist (an atom where the spec still asks for
 * ratios) answers NULL: the level spec is deeper than the array. */
static Expr* ratio_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0) return expr_copy(lst);
    if (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST) {
        /* An atom (or a buffer) at this level: nothing to recurse into. */
        if (!is_ndarray(lst) && !levels_all_zero(levels, num))
            return NULL;                       /* level spec deeper than the array */
        return expr_copy(lst);
    }

    Expr* cur = ratio_n(lst, levels[0]);
    if (num == 1) {
        return cur;
    }

    LLSeq cv;
    llq_open(cur, &cv);
    size_t m = cv.n;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        Expr* sub = llq_peek(&cv, i);
        if (sub->type == EXPR_FUNCTION && !head_is(sub, SYM_List)) {
            /* Below the top level only a List is a level of the array:
             * recursing into the arguments of 3/2 or b/a built
             * Rational[Rational[2, 3]] and a b (Ratios#1 of the 2026-11
             * audit).  An element that is not a sublist fits only an
             * all-zero tail of the level spec. */
            out[i] = levels_all_zero(levels + 1, num - 1) ? expr_copy(sub) : NULL;
        } else {
            out[i] = ratio_levels(sub, levels + 1, num - 1);
        }
        if (!out[i]) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            llq_close(&cv);
            expr_free(cur);
            return NULL;
        }
    }
    Expr* result = llq_result(&cv, out, m);
    free(out);
    llq_close(&cv);
    expr_free(cur);
    return result;
}

Expr* builtin_ratios(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Successive quotients straight off the buffer. Declines an integer dtype
     * (the ratios are exact Rationals) and a zero element (ComplexInfinity),
     * both of which the List path answers correctly. */
    if (is_ndarray(lst)) {
        Expr* nd = ndstruct_ratios(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }

    /* Ratios[assoc] gives the successive value ratios, keyed by the trailing key
     * of each pair (the leading key drops, as in Wolfram and like Differences). */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    if (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST) return NULL;

    /* Ratios[list, {n1, n2, ...}] — per-level ratios. */
    if (argc == 2 && (is_listq(res->data.function.args[1]) ||
                      res->data.function.args[1]->type == EXPR_LINKEDLIST)) {
        Expr* spec = res->data.function.args[1];      /* a List or a LinkedList */
        size_t num = (size_t)ll_seq_length(spec);
        if (num == 0) return ratio_n(lst, 0);

        int64_t* levels = malloc(sizeof(int64_t) * num);
        for (size_t i = 0; i < num; i++) {
            Expr* e = ll_seq_element(spec, (int64_t)i);
            bool ok = e->type == EXPR_INTEGER && e->data.integer >= 0;
            if (ok) levels[i] = e->data.integer;
            expr_free(e);
            if (!ok) {
                free(levels);
                return NULL;
            }
        }
        Expr* result = ratio_levels(lst, levels, num);
        free(levels);
        if (!result) {
            char* sspec = expr_to_string(spec);
            mth_message("Ratios", "lvldep", "The level specification %s is deeper than the array.",
                        sspec ? sspec : "?");
            free(sspec);
        }
        return result;
    }

    /* Ratios[list], Ratios[list, n]. */
    int64_t n = 1;
    if (argc == 2) {
        Expr* an = res->data.function.args[1];
        if (an->type != EXPR_INTEGER || an->data.integer < 0) return NULL;
        n = an->data.integer;
    }

    return ratio_n(lst, n);
}
