#include "list_common.h"
#include "accumulate.h"
#include "assoc.h"
#include "ndarray.h"
#include "ndreduce.h"
#include "../numloop.h"
#include "numeric.h"
#include "checked_int.h"     /* ci_add_i64 */
#include "linkedlist_seq.h"  /* LLSeq: a LinkedList accumulates from its slots */
#include <math.h>

/* Accumulate of a typed LinkedList straight on its slots: gathered into one
 * scratch buffer and prefix-summed exactly as the packed Accumulate does
 * (ndred_accumulate: one running int64 sum, or one running double add per
 * element), so the answers agree to the bit. NULL when the list is boxed,
 * complex or empty, or the sum leaves the machine range (int64 overflow, a
 * non-finite double) -- the element path then answers as the List one does. */
static Expr* accumulate_ll_machine(const LinkedList* ll) {
    LLKind k = ll_kind(ll);
    int64_t n = ll_length(ll);
    if (n == 0 || (k != LL_I64 && k != LL_REAL)) return NULL;
    void* buf = malloc(sizeof(double) * (size_t)n);
    if (!buf) return NULL;
    ll_gather(ll, 0, n, buf);
    bool ok = true;
    if (k == LL_I64) {
        int64_t* v = buf;
        for (int64_t i = 1; i < n && ok; i++) ok = !ci_add_i64(v[i - 1], v[i], &v[i]);
    } else {
        double* v = buf;
        double acc = v[0];
        for (int64_t i = 1; i < n; i++) { acc += v[i]; v[i] = acc; }
        ok = isfinite(acc);
        for (int64_t i = 0; i < n && ok; i++) ok = isfinite(v[i]);
    }
    Expr* r = ok ? expr_new_linkedlist(ll_from_buffer(k, buf, n)) : NULL;
    free(buf);
    return r;
}

/* Accumulate[list] returns the running cumulative totals of list, with the
 * same head and the same length as the input. The intermediate sums are
 * built as Plus[acc, next] and reduced by the evaluator so that integers,
 * arbitrary-precision bignums, machine doubles, lists (matrix columns via
 * Listable Plus), and symbolic expressions all combine correctly.
 *
 * An optional second argument of the form Method -> "CompensatedSummation"
 * triggers Kahan compensated summation when every element reduces to a
 * machine double (EXPR_REAL, EXPR_INTEGER, or EXPR_BIGINT). For other
 * inputs the option is silently ignored and the symbolic accumulation
 * is used. */
static bool accumulate_is_compensated_method(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol.name;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method) return false;
    if (rhs->type != EXPR_STRING) return false;
    return strcmp(rhs->data.string, "CompensatedSummation") == 0;
}

static bool accumulate_to_double(Expr* e, double* out) {
    if (e->type == EXPR_REAL)    { *out = e->data.real;             return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;  return true; }
    if (e->type == EXPR_BIGINT)  { *out = numeric_mpz_to_double(e->data.bigint); return true; }
    return false;
}

Expr* builtin_accumulate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Accumulate[assoc] gives the running totals of the values, keeping every
     * key aligned with its cumulative sum. */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    /* Accumulate[ndarray] prefix-sums the leading axis on the buffer (ndreduce.c). */
    if (argc == 1 && is_ndarray(lst)) return ndred_accumulate(res);

    if (lst->type != EXPR_FUNCTION && lst->type != EXPR_LINKEDLIST) return NULL;

    bool kahan_requested = false;
    if (argc == 2) {
        if (!accumulate_is_compensated_method(res->data.function.args[1])) return NULL;
        kahan_requested = true;
    }

    if (!kahan_requested && lst->type == EXPR_LINKEDLIST) {
        Expr* fast = accumulate_ll_machine(lst->data.llist);
        if (fast) return fast;
    }

    /* The elements: a List's own, or a LinkedList's read from its slots. */
    LLSeq lv;
    llq_open(lst, &lv);
    size_t n = lv.n;

    if (n == 0) {
        Expr* r = llq_linked(&lv) ? llq_result(&lv, NULL, 0) : expr_copy(lst);
        llq_close(&lv);
        return r;
    }

    /* Kahan compensated summation in double precision when every element
     * is a machine number. */
    if (kahan_requested) {
        /* Only a list holding a machine Real is summed in doubles: an
         * all-exact list (integers, bignums) is accumulated exactly below,
         * as without the option -- compensated summation is a rounding
         * remedy, and exact sums have no rounding (it turned
         * {10^20, 1, -10^20} into {1.*^20, 1.*^20, 0.}). */
        bool all_numeric = true, any_real = false;
        for (size_t i = 0; i < n; i++) {
            double tmp;
            Expr* el = llq_peek(&lv, i);
            if (!accumulate_to_double(el, &tmp)) {
                all_numeric = false;
                break;
            }
            if (el->type == EXPR_REAL) any_real = true;
        }
        if (all_numeric && any_real) {
            Expr** out = malloc(sizeof(Expr*) * n);
            if (!out) { llq_close(&lv); return NULL; }
            double sum = 0.0;
            double c = 0.0;
            for (size_t i = 0; i < n; i++) {
                double x = 0.0;
                accumulate_to_double(llq_peek(&lv, i), &x);
                double y = x - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
                out[i] = expr_new_real(sum);
            }
            Expr* result = llq_result(&lv, out, n);
            free(out);
            llq_close(&lv);
            return result;
        }
        /* Mixed/symbolic input: fall through to the standard accumulator. */
    }

    /* Automatic numeric fast-path: an all-machine-Real list accumulates as a
     * running double sum, in place of the evaluate(Plus[..]) per element below.
     * NULL means the list is not uniformly Real -- fall through. */
    if (lst->type == EXPR_FUNCTION) {
        Expr* fast = numloop_accumulate(lst);
        if (fast) { llq_close(&lv); return fast; }
    }

    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) { llq_close(&lv); return NULL; }

    out[0] = expr_copy(llq_peek(&lv, 0));
    for (size_t i = 1; i < n; i++) {
        Expr** plus_args = malloc(sizeof(Expr*) * 2);
        if (!plus_args) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            llq_close(&lv);
            return NULL;
        }
        plus_args[0] = expr_copy(out[i - 1]);
        plus_args[1] = expr_copy(llq_peek(&lv, i));
        Expr* plus_expr = expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
        free(plus_args);
        out[i] = evaluate(plus_expr);
        expr_free(plus_expr);
    }

    Expr* result = llq_result(&lv, out, n);
    free(out);
    llq_close(&lv);
    return result;
}
