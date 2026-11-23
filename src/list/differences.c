#include "list_common.h"
#include "message.h"   /* mth_message */
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "differences.h"
#include "assoc.h"
#include "checked_int.h"     /* ci_sub_i64 */
#include "linkedlist_seq.h"  /* LLSeq: a LinkedList differences from its slots */
#include <math.h>

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
/* One differencing pass over a typed LinkedList's slots, gathered into a
 * scratch buffer: one int64 or one double subtraction per element, as the
 * packed Differences computes it. NULL for a boxed or complex list, an int64
 * overflow, or a non-finite double (the element path answers). */
static Expr* diff_once_ll_machine(const LinkedList* ll, int64_t s) {
    LLKind k = ll_kind(ll);
    int64_t n = ll_length(ll);
    int64_t a = (s < 0) ? (s == INT64_MIN ? INT64_MAX : -s) : s;  /* |INT64_MIN| > any length (audit Differences#1) */
    if (n <= a || (k != LL_I64 && k != LL_REAL)) return NULL;
    void* buf = malloc(sizeof(double) * (size_t)n);
    if (!buf) return NULL;
    ll_gather(ll, 0, n, buf);
    int64_t outn = n - a;
    bool ok = true;
    if (k == LL_I64) {
        int64_t* v = buf;
        for (int64_t i = 0; i < outn && ok; i++)
            ok = !(s > 0 ? ci_sub_i64(v[i + a], v[i], &v[i]) : ci_sub_i64(v[i], v[i + a], &v[i]));
    } else {
        double* v = buf;
        for (int64_t i = 0; i < outn && ok; i++) {
            v[i] = s > 0 ? v[i + a] - v[i] : v[i] - v[i + a];
            ok = isfinite(v[i]);
        }
    }
    Expr* r = ok ? expr_new_linkedlist(ll_from_buffer(k, buf, outn)) : NULL;
    free(buf);
    return r;
}

static Expr* diff_once(Expr* lst, int64_t s) {
    if (lst->type == EXPR_LINKEDLIST) {
        Expr* fast = diff_once_ll_machine(lst->data.llist, s);
        if (fast) return fast;
    }
    LLSeq lv;
    llq_open(lst, &lv);
    size_t n = lv.n;
    int64_t a = (s < 0) ? (s == INT64_MIN ? INT64_MAX : -s) : s;  /* |INT64_MIN| > any length (audit Differences#1) */

    if ((int64_t)n <= a) {
        Expr* r = llq_result(&lv, NULL, 0);
        llq_close(&lv);
        return r;
    }

    size_t outn = n - (size_t)a;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = llq_peek(&lv, i);
        Expr* hi = llq_peek(&lv, i + (size_t)a);
        out[i] = (s > 0) ? diff_minus(hi, lo) : diff_minus(lo, hi);
    }
    Expr* result = llq_result(&lv, out, outn);
    free(out);
    llq_close(&lv);
    return result;
}

/* Apply diff_once with step s a total of n times. Returns a fresh expression
 * (a copy of lst when n == 0). */
static Expr* diff_n_step(Expr* lst, int64_t n, int64_t s) {
    /* Zero passes still answer a NEW LinkedList, never the argument itself. */
    Expr* cur = lst->type == EXPR_LINKEDLIST ? expr_new_linkedlist(ll_copy(lst->data.llist))
                                             : expr_copy(lst);
    for (int64_t k = 0; k < n; k++) {
        Expr* nxt = diff_once(cur, s);
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

/* Multidimensional differences: apply levels[0] step-1 first differences at the
 * top level, then recurse into each element with the remaining levels.
 * A level that does not exist (an atom where the spec still asks for
 * differences) answers NULL: the level spec is deeper than the array. */
static Expr* diff_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0) return expr_copy(lst);
    if (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST) {
        /* An atom (or a buffer) at this level: nothing to recurse into. */
        if (!is_ndarray(lst) && !levels_all_zero(levels, num))
            return NULL;                       /* level spec deeper than the array */
        return expr_copy(lst);
    }

    Expr* cur = diff_n_step(lst, levels[0], 1);
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
            out[i] = diff_levels(sub, levels + 1, num - 1);
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

Expr* builtin_differences(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Native buffer path first; ndstruct_delist_repack below stays the fallback
     * for every form it declines. See ndstruct.h. */
    if (is_ndarray(lst)) {
        Expr* fast = ndstruct_differences(res);
        if (fast) return fast;
        return ndstruct_delist_repack(res, lst);
    }

    /* Differences[assoc] gives the successive value differences, keyed by the
     * trailing key of each pair (so the leading key drops, as in Wolfram). */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    if (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST) return NULL;

    /* Differences[list, {n1, n2, ...}] — per-level differences. */
    if (argc == 2 && (is_listq(res->data.function.args[1]) ||
                      res->data.function.args[1]->type == EXPR_LINKEDLIST)) {
        Expr* spec = res->data.function.args[1];      /* a List or a LinkedList */
        size_t num = (size_t)ll_seq_length(spec);
        if (num == 0) return diff_n_step(lst, 0, 1);

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
        Expr* result = diff_levels(lst, levels, num);
        free(levels);
        if (!result) {
            char* sspec = expr_to_string(spec);
            mth_message("Differences", "lvldep", "The level specification %s is deeper than the array.",
                        sspec ? sspec : "?");
            free(sspec);
        }
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
