#include "../ds_ops.h"
#include "linkedlist.h"
#include "list_common.h"
#include "total.h"
#include "assoc.h"
#include "ndreduce.h"
#include "plus.h"
#include "funcprog.h"   /* mth_level_spec_parse */
#include "sym_names.h"
#include "ndarray.h"
#include <string.h>

/* The depth that negative and Infinity level specs count from: the nesting
 * of Lists only. An element with any other head is a leaf, so
 * Total[{z^2, 1}, Infinity] is 1 + z^2 (descending into Power[z, 2] made it
 * 3 + z) and Total[{{a + 1, 2}}, {-1}] totals the rows. */
static int64_t get_depth_for_total(Expr* e);
static bool total_depth_fn(Expr* el, void* ud) {
    int64_t d = get_depth_for_total(el);
    if (d > *(int64_t*)ud) *(int64_t*)ud = d;
    return true;
}
static int64_t get_depth_for_total(Expr* e) {
    if (e->type == EXPR_NDARRAY) return 1 + e->data.ndarray.rank;
    if (ds_has_each(e->type)) {
        /* A payload sequence (LinkedList, PVector, ...) nests like the List it
         * stands for: a typed one holds machine atoms, so its depth is 2 without
         * a scan; a boxed one is one more than its deepest element. */
        if (ds_ops(e->type)->length(e) == 0) return 1;
        if (ds_atoms_only(e)) return 2;
        int64_t max_d = 0;
        ds_each(e, total_depth_fn, &max_d);
        return 1 + max_d;
    }
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List) return 1;
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_depth_for_total(e->data.function.args[i]);
        if (d > max_d) max_d = d;
    }
    return 1 + max_d;
}

/* Plus of `count` owned summands (the array is the caller's; its entries are
 * adopted), collected the way Total's level 1 always has. */
static Expr* total_plus_of(Expr** plus_args, size_t count) {
    if (count == 0) return expr_new_integer(0);
    {
        Expr* plus_expr = expr_new_function(expr_new_symbol(SYM_Plus), plus_args, count);
        /* The list elements are already evaluated (Total is not HoldAll), so
         * collect like terms directly via builtin_plus — skipping the O(n)
         * re-evaluation of every already-normal term and the redundant
         * ORDERLESS sort of the raw input. builtin_plus does NOT free its
         * argument (only the evaluator frees a builtin's input), so we own and
         * free plus_expr here. eval_and_free then finalizes the small collapsed
         * result: Flat-flattens any nested Plus element and sorts. */
        Expr* collected = builtin_plus(plus_expr);
        if (collected) {
            expr_free(plus_expr);
            return eval_and_free(collected);
        }
        return eval_and_free(plus_expr);  /* builtin declined: evaluate as usual */
    }
}

static Expr* total_at_exactly_level_k(Expr* e, int64_t k);

/* Total of a LinkedList at exactly level k, natively.
 *   k == 1: the sum of the elements. A typed list sums its slots in place
 *           (ll_sum_machine: exact int64 promoting to GMP, or the pairwise
 *           float tree of the packed Total, bit for bit); a boxed list -- or
 *           a machine sum that left the machine range -- collects one Plus
 *           straight from the slots, exactly as a List's level 1 does.
 *   k >= 2: a new LinkedList of each element totalled at level k - 1, the
 *           closure rule (a List gives a List).
 * No List is ever built. */
Expr* total_ll_level_k(Expr* e, int64_t k) {
    LinkedList* ll = e->data.llist;
    int64_t n = ll_length(ll);
    if (k == 1) {
        if (n == 0) return expr_new_integer(0);
        Expr* fast = ll_sum_machine(ll);
        if (fast) return fast;
        Expr** a = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t i = 0; i < n; i++) a[i] = ll_get(ll, i);
        Expr* r = total_plus_of(a, (size_t)n);
        free(a);
        return r;
    }
    LinkedList* out = ll_new();
    for (int64_t i = 0; i < n; i++) {
        Expr* el = ll_get(ll, i);
        ll_push_back(out, total_at_exactly_level_k(el, k - 1));
        expr_free(el);
    }
    return expr_new_linkedlist(out);
}

static Expr* total_at_exactly_level_k(Expr* e, int64_t k) {
    if (k <= 0) return expr_copy(e);
    { Expr* tl = ds_total_level(e, k); if (tl) return tl; }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (k == 1) {
        size_t count = e->data.function.arg_count;
        if (count == 0) return expr_new_integer(0);
        Expr** plus_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) plus_args[i] = expr_copy(e->data.function.args[i]);
        Expr* r = total_plus_of(plus_args, count);
        free(plus_args);
        return r;
    } else {
        size_t count = e->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = total_at_exactly_level_k(e->data.function.args[i], k - 1);
        }
        Expr* res = expr_new_function(expr_copy(e->data.function.head), new_args, count);
        free(new_args);
        return res;
    }
}

/* A trailing Method -> m option of Total.  Returns 1 for
 * "CompensatedSummation", 0 for any other Method value, -1 when e is not a
 * Method rule (so it is the level spec). */
static int total_method_option(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL) return -1;
    const char* hd = e->data.function.head->data.symbol.name;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return -1;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method) return -1;
    const Expr* rhs = e->data.function.args[1];
    return rhs->type == EXPR_STRING &&
           strcmp(rhs->data.string, "CompensatedSummation") == 0;
}

/* Kahan-compensated total of a flat list of machine numbers, or NULL when an
 * element is not a machine real (then the ordinary Plus fold is used: exact
 * and symbolic sums have no rounding to compensate). */
static Expr* total_kahan(const Expr* list) {
    if (is_ndarray(list)) {
        /* A packed / NDArray rank-1 float64 buffer (Total is packed-aware). */
        if (list->data.ndarray.rank != 1 || list->data.ndarray.dtype != NDT_FLOAT64)
            return NULL;
        const double* d = (const double*)list->data.ndarray.data;
        int64_t n = list->data.ndarray.dims[0];
        double sum = 0.0, c = 0.0;
        for (int64_t i = 0; i < n; i++) {
            double y = d[i] - c, t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return expr_new_real(sum);
    }
    if (list->type != EXPR_FUNCTION) return NULL;
    size_t n = list->data.function.arg_count;
    bool any_real = false;
    for (size_t i = 0; i < n; i++) {
        const Expr* x = list->data.function.args[i];
        if (x->type == EXPR_REAL) any_real = true;
        else if (x->type != EXPR_INTEGER) return NULL;
    }
    if (!any_real) return NULL;
    double sum = 0.0, c = 0.0;
    for (size_t i = 0; i < n; i++) {
        const Expr* x = list->data.function.args[i];
        double v = x->type == EXPR_REAL ? x->data.real : (double)x->data.integer;
        double y = v - c, t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return expr_new_real(sum);
}

Expr* builtin_total(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    /* Total[list, (level,) Method -> m].  Read the option off and decide the
     * rest on the positional arguments, so the shared level parser never sees
     * a Method rule (it said Total::level for one, v0.6727). */
    size_t npos = res->data.function.arg_count;
    int method = -1;
    if (npos >= 2) {
        method = total_method_option(res->data.function.args[npos - 1]);
        if (method >= 0) npos--;
    }
    if (npos > 2) return NULL;          /* Total[l, a, b]: arity says nonopt */
    if (method == 1) {
        const Expr* lv = npos == 2 ? res->data.function.args[1] : NULL;
        bool level1 = !lv || (lv->type == EXPR_INTEGER && lv->data.integer == 1);
        if (level1) {
            Expr* k = total_kahan(res->data.function.args[0]);
            if (k) return k;
        }
    }
    if (method >= 0) {
        /* Any other Method, or input Kahan cannot improve: the ordinary sum. */
        Expr** a = malloc(sizeof(Expr*) * npos);
        if (!a) return NULL;
        for (size_t i = 0; i < npos; i++) a[i] = expr_copy(res->data.function.args[i]);
        Expr* call = expr_new_function(expr_new_symbol("Total"), a, npos);
        free(a);
        Expr* r = builtin_total(call);
        expr_free(call);
        return r;                  /* NULL: stays unevaluated */
    }

    /* Total[assoc] sums the association's values. */
    { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path: sum the flat buffer directly (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_total(res);

    Expr* list = res->data.function.args[0];
    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* Total of an empty list is the additive identity, at any level spec:
     * Total[{}] is 0, as Plus @@ {} already was. The level walk below has no
     * element to fold and handed the empty List straight back, so Total[{}]
     * answered {} -- an empty LIST where a number is expected, which then
     * propagates as a non-number through whatever consumes it.
     *
     * Found by the coverage sweep (experiment 20): a Select that matched
     * nothing fed Total, and the checksum came back {} instead of 0. Nothing in
     * the suite covered the empty case because nothing in the suite produced an
     * empty result by accident. */
    if (is_listq(list) && list->data.function.arg_count == 0)
        return expr_new_integer(0);
    /* A LinkedList takes the same level walk below; total_at_exactly_level_k
     * dispatches it to the native kernels (total_ll_level_k). */
    bool is_ll = ds_has_each(list->type) && ds_ops(list->type)->total_level != NULL;
    if (is_ll && ds_ops(list->type)->length(list) == 0) return expr_new_integer(0);

    int64_t n1 = 1, n2 = 1;
    int64_t depth = get_depth_for_total(list);

    /* The shared level-spec parser (funcprog.h); a bad spec is Total::level
     * (audit Total#3: it stayed unevaluated without a word).  Infinity counts
     * as the deepest level that can be totalled, depth - 1, and a negative
     * bound -k as depth - k. */
    if (level_spec) {
        int64_t lo, hi;
        if (!mth_level_spec_parse(level_spec, &lo, &hi)) {
            mth_level_spec_message("Total", level_spec);
            return NULL;
        }
        n1 = lo >= MTH_LEVEL_INF ? depth - 1 : lo < 0 ? depth + lo : lo;
        n2 = hi >= MTH_LEVEL_INF ? depth - 1 : hi < 0 ? depth + hi : hi;
    }

    if (n1 > n2) return expr_copy(list);
    if (n1 < 1) n1 = 1;
    if (n2 >= depth) n2 = depth - 1;

    Expr* current = expr_copy(list);
    for (int64_t k = n2; k >= n1; k--) {
        Expr* next = total_at_exactly_level_k(current, k);
        expr_free(current);
        current = next;
    }

    return current;
}
