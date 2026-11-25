#include "sym_intern.h"   /* intern_symbol: Scaled */
#include "list_common.h"
#include "minmax.h"
#include "assoc.h"
#include "ndarray.h"
#include "ndreduce.h"
#include "interval.h"
#include "numeric.h"
#include "linkedlist_seq.h"  /* LLSeq: LinkedList arguments splice like Lists */

/* Max / Min of a typed LinkedList straight on its slots (the extreme of the
 * raw int64 or double values, as ndred_max / ndred_min scan a packed
 * buffer). NULL for a boxed, complex or empty list. */
static Expr* minmax_ll_machine(const LinkedList* ll, bool want_max) {
    LLKind k = ll_kind(ll);
    int64_t n = ll_length(ll);
    if (n == 0 || (k != LL_I64 && k != LL_REAL)) return NULL;
    LLCursor c;
    const void* p;
    int dir;
    int64_t base, m;
    ll_cursor(&c, ll, 0, n);
    if (k == LL_I64) {
        int64_t best = 0;
        bool first = true;
        while ((m = ll_cursor_next(&c, &p, &dir, &base)) > 0)
            for (int64_t t = 0; t < m; t++) {
                int64_t v = LL_SLOT(int64_t, p, dir, t);
                if (first || (want_max ? v > best : v < best)) best = v;
                first = false;
            }
        return expr_new_integer(best);
    }
    double best = 0.0;
    bool first = true;
    while ((m = ll_cursor_next(&c, &p, &dir, &base)) > 0)
        for (int64_t t = 0; t < m; t++) {
            double v = LL_SLOT(double, p, dir, t);
            if (first || (want_max ? v > best : v < best)) best = v;
            first = false;
        }
    return expr_new_real(best);
}

/* A List argument, or a LinkedList (which stands for one): the arguments
 * Max and Min splice into themselves. */
static bool minmax_spliced(const Expr* arg) {
    return arg->type == EXPR_LINKEDLIST ||
           (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_List);
}

/* Max / Min (`head`) of a LinkedList's elements, read one at a time from a
 * snapshot: head[acc, e] evaluated element by element (Max and Min are Flat,
 * so an unevaluated head[acc, e] splices into the next step and the answer is
 * head[all elements]'s) -- or the typed raw scan. NULL for an empty list,
 * which contributes nothing. */
static Expr* minmax_ll_reduce(const Expr* head, const Expr* arg) {
    const LinkedList* ll = arg->data.llist;
    if (ll_length(ll) == 0) return NULL;
    bool want_max = head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Max;
    Expr* fast = minmax_ll_machine(ll, want_max);
    if (fast) return fast;
    LLSeq q;
    llq_open(arg, &q);
    Expr* acc = expr_copy(llq_peek(&q, 0));
    for (size_t i = 1; i < q.n; i++) {
        Expr* two[2] = { acc, expr_copy(llq_peek(&q, i)) };
        acc = eval_and_free(expr_new_function(expr_copy((Expr*)head), two, 2));
    }
    size_t nq = q.n;
    llq_close(&q);
    if (nq == 1) {                 /* head[e] alone: the evaluator's answer */
        Expr* one[1] = { acc };
        acc = eval_and_free(expr_new_function(expr_copy((Expr*)head), one, 1));
    }
    return acc;
}

/* head[args...] with every List argument replaced by its elements and every
 * LinkedList argument by the reduction of its elements (minmax_ll_reduce). */
static Expr* minmax_splice(Expr* res) {
    size_t n = res->data.function.arg_count;
    size_t new_count = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        new_count += (minmax_spliced(arg) && arg->type != EXPR_LINKEDLIST) ? arg->data.function.arg_count : 1;
    }
    Expr** new_args = malloc(sizeof(Expr*) * (new_count ? new_count : 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_LINKEDLIST) {
            Expr* r = minmax_ll_reduce(res->data.function.head, arg);
            if (r) new_args[k++] = r;
        } else if (minmax_spliced(arg)) {
            for (size_t j = 0; j < arg->data.function.arg_count; j++)
                new_args[k++] = expr_copy(arg->data.function.args[j]);
        } else {
            new_args[k++] = expr_copy(arg);
        }
    }
    Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, k);
    free(new_args);
    return ret;
}

/* A real numeric quantity that is not a number: Pi, Sqrt[2], Log[3]. Max and
 * Min order these with the numbers, keeping the exact form of the winner
 * (Max[1., Pi] is Pi, Min[1/2, Sqrt[2]] is 1/2). A complex quantity (I Pi)
 * has no order and stays symbolic. */
static bool minmax_real_quantity(Expr* e) {
    if (e->type != EXPR_FUNCTION && e->type != EXPR_SYMBOL) return false;
    if (!expr_is_numeric_quantity(e)) return false;
    Expr* v = numericalize(e, numeric_machine_spec());
    bool ok = v && (v->type == EXPR_REAL
#ifdef USE_MPFR
                    || v->type == EXPR_MPFR
#endif
                    );
    if (v) expr_free(v);
    return ok;
}

static bool minmax_orderable(Expr* e) {
    return is_real_numeric(e) || minmax_real_quantity(e);
}

/* Sign of a - b for two orderable arguments. Two numbers keep the canonical
 * order the heads always used; a numeric quantity is decided by Less /
 * Greater, which settle exact constants to whatever precision they need.
 * An undecided pair (equal values) reports 0, so the earlier argument stays. */
static int minmax_cmp(Expr* a, Expr* b) {
    if (is_real_numeric(a) && is_real_numeric(b)) return expr_compare(a, b);
    Expr* args[2] = { expr_copy(a), expr_copy(b) };
    Expr* lt = eval_and_free(expr_new_function(expr_new_symbol(SYM_Less), args, 2));
    bool is_lt = lt->type == EXPR_SYMBOL && lt->data.symbol.name == SYM_True;
    expr_free(lt);
    if (is_lt) return -1;
    Expr* args2[2] = { expr_copy(a), expr_copy(b) };
    Expr* gt = eval_and_free(expr_new_function(expr_new_symbol(SYM_Greater), args2, 2));
    bool is_gt = gt->type == EXPR_SYMBOL && gt->data.symbol.name == SYM_True;
    expr_free(gt);
    return is_gt ? 1 : 0;
}

/* MinMax[list] gives {Min[list], Max[list]} in the natural single pass a caller
 * would otherwise write by hand. Over an association it uses the values (Min and
 * Max already thread there). Delegating to Min / Max keeps the numeric handling
 * — bignums, reals, symbolic extrema, empty-list Infinity — in exactly one place. */
/* MinMax[list, d]: the range {min, max} widened to {min - dmin, max + dmax}
 * for d (both sides), {dmin, dmax}, or Scaled[s] (s times the range on each
 * side) -- audit MinMax#1. */
static Expr* minmax_padded(Expr* res) {
    Expr* d = res->data.function.args[1];
    Expr* one[1] = { expr_copy(res->data.function.args[0]) };
    Expr* call = expr_new_function(expr_new_symbol("MinMax"), one, 1);
    Expr* mm = builtin_minmax(call);
    expr_free(call);
    if (!mm || !head_is(mm, SYM_List) || mm->data.function.arg_count != 2) return mm ? mm : NULL;
    Expr* mn = mm->data.function.args[0];
    Expr* mx = mm->data.function.args[1];
    Expr *dlo, *dhi;
    if (head_is(d, SYM_List) && d->data.function.arg_count == 2) {
        dlo = expr_copy(d->data.function.args[0]);
        dhi = expr_copy(d->data.function.args[1]);
    } else if (head_is(d, intern_symbol("Scaled")) && d->data.function.arg_count == 1) {
        Expr* neg[2] = { expr_new_integer(-1), expr_copy(mn) };
        Expr* span[2] = { expr_copy(mx), expr_new_function(expr_new_symbol(SYM_Times), neg, 2) };
        Expr* t[2] = { expr_copy(d->data.function.args[0]),
                       expr_new_function(expr_new_symbol(SYM_Plus), span, 2) };
        dlo = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t, 2));
        dhi = expr_copy(dlo);
    } else if (!head_is(d, SYM_List)) {
        dlo = expr_copy(d);
        dhi = expr_copy(d);
    } else {
        expr_free(mm);
        return NULL;
    }
    Expr* nlo[2] = { expr_new_integer(-1), dlo };
    Expr* lo[2] = { expr_copy(mn), expr_new_function(expr_new_symbol(SYM_Times), nlo, 2) };
    Expr* hi[2] = { expr_copy(mx), dhi };
    Expr* items[2] = { eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), lo, 2)),
                       eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), hi, 2)) };
    expr_free(mm);
    return expr_new_function(expr_new_symbol(SYM_List), items, 2);
}

Expr* builtin_minmax(Expr* res) {
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count == 2) return minmax_padded(res);
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    bool is_list = arg->type == EXPR_FUNCTION &&
                   arg->data.function.head->type == EXPR_SYMBOL &&
                   arg->data.function.head->data.symbol.name == SYM_List;
    /* A packed buffer is accepted here so that MinMax can go on pack.c's AWARE
     * list. Without both halves the head is the recurring defect of
     * docs/experiments/README.md wearing another costume: the transparency gate
     * materialised 10^6 Expr nodes, and the Min and Max calls built below --
     * each of which HAS a buffer path (ndred_min / ndred_max) -- then ran the
     * generic List code on the boxed copy. Measured 313 ms on a 10^6 vector
     * against 426 us for np.array([v.min(), v.max()]), 735x, while the same
     * Min[v] and Max[v] called directly cost 255 us and 232 us.
     *
     * Nothing below needs changing: expr_copy is a refcount bump, so the buffer
     * is handed to Min and Max intact and each takes its own fast path. */
    if (!is_list && !is_ndarray(arg) && !is_association(arg) && arg->type != EXPR_LINKEDLIST)
        return NULL;

    Expr* min_arg[1] = { expr_copy(arg) };
    Expr* min_call = expr_new_function(expr_new_symbol(SYM_Min), min_arg, 1);
    Expr* max_arg[1] = { expr_copy(arg) };
    Expr* max_call = expr_new_function(expr_new_symbol(SYM_Max), max_arg, 1);
    Expr* mn = evaluate(min_call); expr_free(min_call);
    Expr* mx = evaluate(max_call); expr_free(max_call);

    Expr* items[2] = { mn, mx };
    return expr_new_function(expr_new_symbol(SYM_List), items, 2);
}

/* Min / Max over several arguments, some of them arrays (packed or visible).
 * The heads flatten Lists, so Max[{a, b}, c] is Max[Max[{a, b}], c]: reduce
 * each array argument through its own buffer path (ndred_max / ndred_min,
 * which also answers -Infinity / Infinity for an empty array), then evaluate
 * the head over the reduced values and the other arguments. Without this a
 * visible NDArray stayed an atom and the call was left unevaluated, and a
 * packed one was materialised element by element (Max[v, 2] over 10^6
 * elements took ~200 ms). NULL when no argument is an array. */
static Expr* minmax_reduce_array_args(Expr* res, const char* head) {
    size_t n = res->data.function.arg_count;
    bool any = false;
    for (size_t i = 0; i < n && !any; i++) any = is_ndarray(res->data.function.args[i]);
    if (!any) return NULL;
    Expr** args = malloc(sizeof(Expr*) * n);
    if (!args) return NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        if (is_ndarray(a)) {
            Expr* one[1] = { expr_copy(a) };
            Expr* call = expr_new_function(expr_new_symbol(head), one, 1);
            args[i] = evaluate(call);
            expr_free(call);
        } else {
            args[i] = expr_copy(a);
        }
    }
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    free(args);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;
}

/* Min/Max with Interval arguments (audit Max#1).  An interval enters an
 * extremum only through its extreme endpoint: Min sees the lowest lower
 * endpoint (the first pair's, after canonical sorting), Max the highest upper
 * endpoint (the last pair's).  Each Interval argument is replaced by that
 * endpoint and the call rebuilt, so Max[Interval[{1,3}], 5] -> Max[3, 5] -> 5
 * and Min[Interval[{1,3}]] -> Min[1] -> 1.  Returns NULL when no argument is an
 * Interval. */
static Expr* minmax_interval_endpoints(Expr* res, bool is_max) {
    size_t n = res->data.function.arg_count;
    bool any = false;
    for (size_t i = 0; i < n && !any; i++) {
        Expr* a = res->data.function.args[i];
        any = is_interval(a) && interval_pair_count(a) > 0;
    }
    if (!any) return NULL;
    Expr** args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        if (is_interval(a) && interval_pair_count(a) > 0)
            args[i] = expr_copy(is_max ? interval_pair_hi(a, interval_pair_count(a) - 1)
                                       : interval_pair_lo(a, 0));
        else
            args[i] = expr_copy(a);
    }
    Expr* out = expr_new_function(expr_copy(res->data.function.head), args, n);
    free(args);
    return out;
}

Expr* builtin_min(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_Infinity);

    /* Min[assoc] is the minimum of the association's values. */
    if (n == 1 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
    }

    /* Min[ndarray] scans the flat buffer (see ndreduce.c). */
    if (n == 1 && is_ndarray(res->data.function.args[0])) return ndred_min(res);
    /* Min[ll] of a typed LinkedList scans its slots. */
    if (n == 1 && res->data.function.args[0]->type == EXPR_LINKEDLIST) {
        Expr* fast = minmax_ll_machine(res->data.function.args[0]->data.llist, false);
        if (fast) return fast;
    }
    /* Min[nd, x, ...]: see minmax_reduce_array_args. */
    if (n > 1) { Expr* r = minmax_reduce_array_args(res, SYM_Min); if (r) return r; }

    /* Min[Interval[{a,b}, ...]] is the lowest lower endpoint (the first pair's
     * lower endpoint after canonical sorting). */
    { Expr* r = minmax_interval_endpoints(res, false); if (r) return r; }

    // Check for List (or LinkedList) arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n && !has_list; i++) has_list = minmax_spliced(res->data.function.args[i]);
    if (has_list) return minmax_splice(res);

    /* An undetermined element leaves the extremum undetermined (audit
     * MinMax#2): Min/Max are numeric functions of every element. */
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Indeterminate)
            return expr_new_symbol(SYM_Indeterminate);
    }

    /* Min of a single element is that element (Min[5] -> 5, Min[x] -> x). Lists
     * were flattened above, so a lone argument here is a genuine scalar. */
    if (n == 1) return expr_copy(res->data.function.args[0]);

    // Check for Overflow[] and -Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        if (is_minus_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* min_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (minmax_orderable(arg)) {
            if (min_num == NULL) {
                min_num = expr_copy(arg);
            } else {
                if (minmax_cmp(arg, min_num) < 0) {
                    expr_free(min_num);
                    min_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (min_num ? 1 : 0);
        if (final_count == 0) {
            if (min_num) expr_free(min_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return expr_new_symbol(SYM_Infinity);
        }
        if (final_count == 1) {
            Expr* single = min_num ? min_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (min_num) final_args[k++] = min_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (min_num) expr_free(min_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}

Expr* builtin_max(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return make_minus_infinity();

    /* Max[assoc] is the maximum of the association's values. */
    if (n == 1 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
    }

    /* Max[ndarray] scans the flat buffer (see ndreduce.c). */
    if (n == 1 && is_ndarray(res->data.function.args[0])) return ndred_max(res);
    /* Max[ll] of a typed LinkedList scans its slots. */
    if (n == 1 && res->data.function.args[0]->type == EXPR_LINKEDLIST) {
        Expr* fast = minmax_ll_machine(res->data.function.args[0]->data.llist, true);
        if (fast) return fast;
    }
    /* Max[nd, x, ...]: see minmax_reduce_array_args. */
    if (n > 1) { Expr* r = minmax_reduce_array_args(res, SYM_Max); if (r) return r; }

    /* Max[Interval[..., {a,b}]] is the highest upper endpoint (the last pair's
     * upper endpoint after canonical sorting). */
    { Expr* r = minmax_interval_endpoints(res, true); if (r) return r; }

    // Check for List (or LinkedList) arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n && !has_list; i++) has_list = minmax_spliced(res->data.function.args[i]);
    if (has_list) return minmax_splice(res);

    /* An undetermined element leaves the extremum undetermined (audit
     * MinMax#2): Min/Max are numeric functions of every element. */
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Indeterminate)
            return expr_new_symbol(SYM_Indeterminate);
    }


    /* Max of a single element is that element (Max[5] -> 5, Max[x] -> x). Lists
     * were flattened above, so a lone argument here is a genuine scalar. */
    if (n == 1) return expr_copy(res->data.function.args[0]);

    // Check for Overflow[] and Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        if (is_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* max_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (minmax_orderable(arg)) {
            if (max_num == NULL) {
                max_num = expr_copy(arg);
            } else {
                if (minmax_cmp(arg, max_num) > 0) {
                    expr_free(max_num);
                    max_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_minus_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (max_num ? 1 : 0);
        if (final_count == 0) {
            if (max_num) expr_free(max_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return make_minus_infinity();
        }
        if (final_count == 1) {
            Expr* single = max_num ? max_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (max_num) final_args[k++] = max_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (max_num) expr_free(max_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}
