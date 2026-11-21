#include "list_common.h"
#include "ndarray.h"    /* is_ndarray, ndarray_delist_and_reeval */
#include "ndstruct.h"   /* NDArray buffer fast paths */
#include "ratios.h"
#include "assoc.h"
#include "linkedlist_seq.h"  /* LLView: a LinkedList's ratios from its slots */
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
    LLView lv;
    ll_view(lst, &lv);
    size_t n = lv.n;

    if (n <= 1) {
        Expr* r = ll_view_result(&lv, NULL, 0);
        ll_view_done(&lv);
        return r;
    }

    size_t outn = n - 1;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = lv.a[i];
        Expr* hi = lv.a[i + 1];
        out[i] = ratio_divide(hi, lo);
    }
    Expr* result = ll_view_result(&lv, out, outn);
    free(out);
    ll_view_done(&lv);
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

/* Multidimensional ratios: apply levels[0] first ratios at the top level, then
 * recurse into each element with the remaining levels. */
static Expr* ratio_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0 || (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST)) {
        return expr_copy(lst);
    }

    Expr* cur = ratio_n(lst, levels[0]);
    if (num == 1) {
        return cur;
    }

    LLView cv;
    ll_view(cur, &cv);
    size_t m = cv.n;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = ratio_levels(cv.a[i], levels + 1, num - 1);
    }
    Expr* result = ll_view_result(&cv, out, m);
    free(out);
    ll_view_done(&cv);
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
