#include "list_common.h"
#include "ratios.h"
#include "assoc.h"

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
static Expr* ratio_once(Expr* lst) {
    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;

    if (n <= 1) {
        return expr_new_function(expr_copy(head), NULL, 0);
    }

    size_t outn = n - 1;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = lst->data.function.args[i];
        Expr* hi = lst->data.function.args[i + 1];
        out[i] = ratio_divide(hi, lo);
    }
    Expr* result = expr_new_function(expr_copy(head), out, outn);
    free(out);
    return result;
}

/* Apply ratio_once a total of n times. Returns a fresh expression (a copy of
 * lst when n == 0). */
static Expr* ratio_n(Expr* lst, int64_t n) {
    Expr* cur = expr_copy(lst);
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
    if (num == 0 || lst->type != EXPR_FUNCTION) {
        return expr_copy(lst);
    }

    Expr* cur = ratio_n(lst, levels[0]);
    if (num == 1) {
        return cur;
    }

    size_t m = cur->data.function.arg_count;
    Expr* head = cur->data.function.head;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = ratio_levels(cur->data.function.args[i], levels + 1, num - 1);
    }
    Expr* result = expr_new_function(expr_copy(head), out, m);
    free(out);
    expr_free(cur);
    return result;
}

Expr* builtin_ratios(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Ratios[assoc] gives the successive value ratios, keyed by the trailing key
     * of each pair (the leading key drops, as in Wolfram and like Differences). */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    if (lst->type != EXPR_FUNCTION) return NULL;

    /* Ratios[list, {n1, n2, ...}] — per-level ratios. */
    if (argc == 2 && is_listq(res->data.function.args[1])) {
        Expr* spec = res->data.function.args[1];
        size_t num = spec->data.function.arg_count;
        if (num == 0) return expr_copy(lst);

        int64_t* levels = malloc(sizeof(int64_t) * num);
        for (size_t i = 0; i < num; i++) {
            Expr* e = spec->data.function.args[i];
            if (e->type != EXPR_INTEGER || e->data.integer < 0) {
                free(levels);
                return NULL;
            }
            levels[i] = e->data.integer;
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
