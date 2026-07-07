#include "list_common.h"
#include "differences.h"
#include "assoc.h"

/* ---- Differences ----------------------------------------------------------
 *
 * Differences[list]              successive (first) differences of list.
 * Differences[list, n]           n-fold first differences (length l - n).
 * Differences[list, n, s]        n-fold differences of elements step s apart
 *                                (length l - n|s|).
 * Differences[list, {n1, n2,..}] successive n_k-th differences at level k of a
 *                                nested list / array.
 *
 * The element-wise subtraction list[[i+s]] - list[[i]] is built as
 * Subtract[hi, lo] and reduced by the evaluator, so integers, bignums, machine
 * doubles, symbolic terms, and whole sublists (matrix rows, via the Listable
 * Plus/Times that Subtract rewrites to) all combine correctly. Because row
 * subtraction threads element-wise, Differences[m, n] on a matrix m takes
 * differences of successive rows, matching Differences[m, {n, 0}].
 */

/* Returns the evaluated expression (minuend - subtrahend). The operands are
 * copied; ownership of the result transfers to the caller. */
static Expr* diff_minus(Expr* minuend, Expr* subtrahend) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(minuend);
    args[1] = expr_copy(subtrahend);
    Expr* sub = expr_new_function(expr_new_symbol(SYM_Subtract), args, 2);
    free(args);
    return eval_and_free(sub);
}

/* One difference pass with nonzero integer step s over the top level of lst
 * (an EXPR_FUNCTION). The result keeps lst's head. For s > 0 the i-th element
 * is elem[i+s] - elem[i]; for s < 0 it is elem[i] - elem[i+|s|]. A list whose
 * length does not exceed |s| yields the empty list. */
static Expr* diff_once(Expr* lst, int64_t s) {
    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;
    int64_t a = (s < 0) ? -s : s;

    if ((int64_t)n <= a) {
        return expr_new_function(expr_copy(head), NULL, 0);
    }

    size_t outn = n - (size_t)a;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = lst->data.function.args[i];
        Expr* hi = lst->data.function.args[i + (size_t)a];
        out[i] = (s > 0) ? diff_minus(hi, lo) : diff_minus(lo, hi);
    }
    Expr* result = expr_new_function(expr_copy(head), out, outn);
    free(out);
    return result;
}

/* Apply diff_once with step s a total of n times. Returns a fresh expression
 * (a copy of lst when n == 0). */
static Expr* diff_n_step(Expr* lst, int64_t n, int64_t s) {
    Expr* cur = expr_copy(lst);
    for (int64_t k = 0; k < n; k++) {
        Expr* nxt = diff_once(cur, s);
        expr_free(cur);
        cur = nxt;
    }
    return cur;
}

/* Multidimensional differences: apply levels[0] step-1 first differences at the
 * top level, then recurse into each element with the remaining levels. */
static Expr* diff_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0 || lst->type != EXPR_FUNCTION) {
        return expr_copy(lst);
    }

    Expr* cur = diff_n_step(lst, levels[0], 1);
    if (num == 1) {
        return cur;
    }

    size_t m = cur->data.function.arg_count;
    Expr* head = cur->data.function.head;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = diff_levels(cur->data.function.args[i], levels + 1, num - 1);
    }
    Expr* result = expr_new_function(expr_copy(head), out, m);
    free(out);
    expr_free(cur);
    return result;
}

Expr* builtin_differences(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Differences[assoc] gives the successive value differences, keyed by the
     * trailing key of each pair (so the leading key drops, as in Wolfram). */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    if (lst->type != EXPR_FUNCTION) return NULL;

    /* Differences[list, {n1, n2, ...}] — per-level differences. */
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
        Expr* result = diff_levels(lst, levels, num);
        free(levels);
        return result;
    }

    /* Differences[list], Differences[list, n], Differences[list, n, s]. */
    int64_t n = 1, s = 1;
    if (argc >= 2) {
        Expr* an = res->data.function.args[1];
        if (an->type != EXPR_INTEGER || an->data.integer < 0) return NULL;
        n = an->data.integer;
    }
    if (argc == 3) {
        Expr* as = res->data.function.args[2];
        if (as->type != EXPR_INTEGER || as->data.integer == 0) return NULL;
        s = as->data.integer;
    }

    return diff_n_step(lst, n, s);
}
