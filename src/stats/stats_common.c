/* stats_common.c -- helpers shared across the src/stats/ per-builtin files.
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats_common.h"
#include "numeric.h"   /* expr_is_numeric_quantity */
#include <float.h>
#include <math.h>
#include "arithmetic.h"   /* is_rational, make_rational */
#include "complex.h"      /* is_complex */
#include "eval.h"         /* evaluate, eval_and_free */
#include "sym_names.h"    /* SYM_* */
#include "assoc.h"        /* is_association, assoc_apply_over_values */
#include "ndreduce.h"     /* ndred_call_has_ndarray, ndred_skewness, ndred_kurtosis */
#include "numeric_order.h" /* numord_key / numord_cmp: stats_sort_order */
#include "linkedlist.h"   /* ll_seq_length: stats_is_matrix */
#include "message.h"      /* mth_message: stats_msg_rectn */
#include "print.h"        /* expr_to_string */
#include <stdlib.h>
#include <math.h>         /* floor -- stats_quantile_point index fallback */

/* True when e is a function call with the given symbol head. */
static bool stats_head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

bool stats_is_numeric(Expr* e, double* val, bool* out_complex) {
    if (e->type == EXPR_INTEGER) {
        if (val) *val = (double)e->data.integer;
        if (out_complex) *out_complex = false;
        return true;
    }
    if (e->type == EXPR_REAL) {
        if (val) *val = e->data.real;
        if (out_complex) *out_complex = false;
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        if (val) *val = (double)n / (double)d;
        if (out_complex) *out_complex = false;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (out_complex) *out_complex = true;
        return true;
    }
    return false;
}

bool stats_is_matrix(const Expr* data) {
    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List ||
        data->data.function.arg_count == 0) return false;
    int64_t cols = ll_seq_length(data->data.function.args[0]);
    if (cols < 0) return false;
    for (size_t i = 1; i < data->data.function.arg_count; i++)
        if (ll_seq_length(data->data.function.args[i]) != cols) return false;
    return true;
}

void stats_msg_rectn(const char* head, Expr* res) {
    char* str = expr_to_string(res);
    mth_message(head, "rectn", "Rectangular array of real numbers is expected at position 1 in %s.", str);
    free(str);
}

bool stats_reject_ragged(const char* head, Expr* res, const Expr* data) {
    if (!stats_head_is(data, SYM_List)) return false;
    for (size_t i = 0; i < data->data.function.arg_count; i++) {
        if (ll_seq_length(data->data.function.args[i]) >= 0) {
            /* A full array of higher rank (a list of matrices) is fine: the
             * statistic is taken elementwise over the first level. */
            Expr* aq = eval_and_free(expr_new_function(expr_new_symbol(SYM_ArrayQ),
                                                       (Expr*[]){ expr_copy((Expr*)data) }, 1));
            bool full = aq->type == EXPR_SYMBOL && aq->data.symbol.name == SYM_True;
            expr_free(aq);
            if (full) return false;
            stats_msg_rectn(head, res);
            return true;
        }
    }
    return false;
}

bool stats_reject_short(const char* head, const Expr* data) {
    if (!stats_head_is(data, SYM_List) || data->data.function.arg_count >= 2) return false;
    char* str = expr_to_string((Expr*)data);
    mth_message(head, "shlen", "The argument %s should have at least two elements.", str);
    free(str);
    return true;
}

Expr* stats_apply_columnwise(const char* func_name, Expr* matrix) {
    if (!stats_is_matrix(matrix)) return NULL;     /* ragged: no Transpose */
    // Result is Map[func_name, Transpose[matrix]]
    Expr* transpose_args[1] = { expr_copy(matrix) };
    Expr* transpose_call = expr_new_function(expr_new_symbol(SYM_Transpose), transpose_args, 1);
    Expr* transposed = evaluate(transpose_call);
    expr_free(transpose_call);

    /* A List, not merely a call: an unevaluated Transpose[...] mapped over
     * re-enters func_name with the same data, without end. */
    if (!stats_head_is(transposed, SYM_List)) {
        expr_free(transposed);
        return NULL;
    }

    Expr* map_args[2] = { expr_new_symbol(func_name), transposed };
    Expr* map_call = expr_new_function(expr_new_symbol(SYM_Map), map_args, 2);
    Expr* result = evaluate(map_call);
    expr_free(map_call);
    return result;
}

Expr* stats_standardized_moment(Expr* res, int p) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Skewness[assoc] / Kurtosis[assoc] work on the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray / packed fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res))
        return (p == 3) ? ndred_skewness(res) : ndred_kurtosis(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c);
     * otherwise CentralMoment below reads the LinkedList itself. */
    { Expr* r; if (stats_ll_arm(res, p == 3 ? STATS_SKEWNESS : STATS_KURTOSIS, &r)) return r; }

    /* A machine vector (a Real among the numbers): the buffer kernel, so a
     * short boxed list answers what its packed / NDArray form does, and data
     * whose moments leave the double range (1.*^-200 cubes to 0.) is
     * rescaled -- Skewness[{1.*^-200, 2.*^-200, 4.*^-200}] is
     * Skewness[{1, 2, 4}] // N, not Indeterminate. */
    if (stats_head_is(data, SYM_List) && data->data.function.arg_count >= 2) {
        size_t n = data->data.function.arg_count;
        bool machine = false, ok = true;
        for (size_t i = 0; i < n && ok; i++) {
            Expr* x = data->data.function.args[i];
            if (x->type == EXPR_REAL) machine = true;
            else if (x->type != EXPR_INTEGER) ok = false;
        }
        if (ok && machine) {
            double* xs = malloc(sizeof(double) * n);
            if (xs) {
                for (size_t i = 0; i < n; i++) stats_is_numeric(data->data.function.args[i], &xs[i], NULL);
                double k = ndred_stdmoment_f64(xs, n, p);
                free(xs);
                if (isfinite(k)) return expr_new_real(k);
            }
        }
    }
    /* m_p = CentralMoment[data, p] and m_2 = CentralMoment[data, 2]. Both reduce
     * for any list/vector/matrix/array (numeric, exact, or symbolic); if either
     * stays a CentralMoment[...] the data was not reducible, so leave the caller's
     * head unevaluated too. */
    Expr* mp = eval_and_free(expr_new_function(expr_new_symbol(SYM_CentralMoment),
                   (Expr*[]){ expr_copy(data), expr_new_integer(p) }, 2));
    Expr* m2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_CentralMoment),
                   (Expr*[]){ expr_copy(data), expr_new_integer(2) }, 2));
    if (stats_head_is(mp, SYM_CentralMoment) || stats_head_is(m2, SYM_CentralMoment)) {
        expr_free(mp); expr_free(m2);
        return NULL;
    }

    /* result = m_p / m_2^(p/2). Power and Divide thread, so a columnwise vector
     * of moments yields a columnwise vector of standardized moments. */
    Expr* denom = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ m2, make_rational(p, 2) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Divide),
               (Expr*[]){ mp, denom }, 2));
}

bool stats_is_real_numeric(Expr* e) {
    /* Leaf fast path. Every head in this subsystem calls this ONCE PER LIST
     * ELEMENT to gate its input, and the general path below costs two full
     * evaluator round-trips (NumericQ, then FreeQ) per call -- 400,000 of them
     * for a 200,000-element list. A machine integer, a machine real, a bignum,
     * an MPFR real and an exact Rational are all NumericQ and all free of I by
     * construction, so they need no evaluator at all. Anything else (symbols,
     * Pi, Sqrt[2], I, unevaluated heads) still goes the long way, so the
     * accepted set is unchanged. */
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (is_rational_like(e)) return true;   /* Integer, BigInt, Rational[p, q] */

    /* An ALREADY-EVALUATED complex number is a hole in the FreeQ[e, I] test
     * below: 2 + I evaluates to the structural node Complex[2, 1], which is
     * NumericQ and contains no literal symbol I, so FreeQ answered True and
     * this gate passed it. Every head in this subsystem then produced a silent
     * answer for data it had just promised was a rectangular array of REAL
     * numbers -- Median[{1, 2 + I, 3}] returned 3 with no message, and
     * Quartiles[{1, 2 + I, 3, 4}] returned {2, 7/2, 3 + I/2}, a complex
     * quartile.
     *
     * Decide it on the imaginary part rather than on the head: Complex[x, 0] is
     * a real number wearing a Complex head, and it DOES reach here. The
     * evaluator normalises an int64 or double zero away (builtin_complex), but
     * not an MPFR one, so N[2+I, 30] + N[Conjugate[2+I], 30] arrives as
     * Complex[4.0, 0.0] at 30 digits and is real. N[] is taken for exactly the
     * types stats_is_numeric cannot read directly, the same fallback the
     * quantile engine uses for h and q.
     *
     * KNOWN REMAINING GAP, stated rather than smoothed over: this closes the
     * bare-Complex case only. A complex value nested under a numeric head --
     * Sqrt[2 + I] is Power[Complex[2,1], 1/2] -- is still NumericQ, still free
     * of literal I, and still passes. Closing that needs a real-valuedness test
     * (Im[e] == 0) rather than a structural one, which is a larger change than
     * this one. */
    Expr* im_part = NULL;
    if (is_complex(e, NULL, &im_part)) {
        double imv = 0.0;
        bool known = stats_is_numeric(im_part, &imv, NULL);
        if (!known) {
            Expr* nim = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
                            (Expr*[]){expr_copy(im_part)}, 1));
            known = stats_is_numeric(nim, &imv, NULL);
            expr_free(nim);
        }
        return known && imv == 0.0;
    }

    /* The numeric-quantity test itself, not the NumericQ builtin: an infinity
     * is a legitimate real datum here (Median[{-Infinity, 5, 1}] is 1), while
     * NumericQ[Infinity] is False (audit NumericQ#1).  It also saves an
     * evaluator round trip per element. */
    if (!expr_is_numeric_quantity(e)) return false;

    Expr* freeq = expr_new_function(expr_new_symbol(SYM_FreeQ), (Expr*[]){expr_copy(e), expr_new_symbol(SYM_I)}, 2);
    Expr* freeq_eval = evaluate(freeq);
    expr_free(freeq);
    if (freeq_eval->type != EXPR_SYMBOL || freeq_eval->data.symbol.name != SYM_True) {
        expr_free(freeq_eval);
        return false;
    }
    expr_free(freeq_eval);

    return true;
}

Expr* stats_quantile_point(const StatsOrder* sorted, size_t n, Expr* q,
                           Expr* a, Expr* b, Expr* c, Expr* d) {
    Expr* n_expr = expr_new_integer((int64_t)n);
    Expr* n_plus_b = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){n_expr, expr_copy(b)}, 2));
    Expr* times_q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){n_plus_b, expr_copy(q)}, 2));
    Expr* h = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(a), times_q}, 2));

    double h_val = 0;
    if (!stats_is_numeric(h, &h_val, NULL)) {
        /* h is NumericQ but not machine-representable (an exact irrational such
         * as Pi/4 * n): take its numeric value via N[h] for the clamp and index
         * decisions only. The result itself is still built from exact h. */
        Expr* nh = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){expr_copy(h)}, 1));
        bool ok = stats_is_numeric(nh, &h_val, NULL);
        expr_free(nh);
        if (!ok) {
            expr_free(h);
            return expr_new_symbol(SYM_Indeterminate);
        }
    }

    if (h_val <= 1.0) {
        expr_free(h);
        return expr_copy(stats_order_at(sorted, 0));
    }
    if (h_val >= (double)n) {
        expr_free(h);
        return expr_copy(stats_order_at(sorted, n - 1));
    }

    Expr* j_expr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Floor), (Expr*[]){expr_copy(h)}, 1));
    int64_t j_idx = 0;
    if (j_expr->type == EXPR_INTEGER) j_idx = j_expr->data.integer;
    else j_idx = (int64_t)floor(h_val);
    expr_free(j_expr);

    if (j_idx < 1) j_idx = 1;
    if (j_idx >= (int64_t)n) j_idx = n - 1;

    Expr* j_expr2 = expr_new_integer(j_idx);
    Expr* neg_j = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), j_expr2}, 2));
    Expr* g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(h), neg_j}, 2));
    expr_free(h);

    /* Wolfram's definition interpolates between x_(Floor[h]) and x_(Ceiling[h]).
     * At integer h the two neighbours COINCIDE, so the upper index is j, not
     * j+1 -- consulting j+1 there is what let a Real neighbour leak into an
     * otherwise exact result. */
    bool g_is_zero = (g->type == EXPR_INTEGER && g->data.integer == 0) ||
                     (g->type == EXPR_REAL && g->data.real == 0.0);
    int64_t upper_idx = g_is_zero ? j_idx : j_idx + 1;

    Expr* d_times_g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(d), expr_copy(g)}, 2));
    Expr* g_weight = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(c), d_times_g}, 2));
    expr_free(g);

    /* SELECT rather than recompute at the two weights that name an element
     * outright. weight == 1 is Quantile's own default (c=1, d=0): computing
     * A[j] + 1*(A[j+1] - A[j]) is an identity in exact arithmetic but NOT in
     * floating point -- with a huge negative first element the subtraction
     * rounds and the sum returns 0. rather than the element. It also drags an
     * exact list to Real whenever q is inexact. weight == 0 is the same
     * argument for the lower neighbour (Quartiles' c=0 at integer h). */
    bool w_is_zero = (g_weight->type == EXPR_INTEGER && g_weight->data.integer == 0) ||
                     (g_weight->type == EXPR_REAL && g_weight->data.real == 0.0);
    bool w_is_one = (g_weight->type == EXPR_INTEGER && g_weight->data.integer == 1) ||
                    (g_weight->type == EXPR_REAL && g_weight->data.real == 1.0);
    if (w_is_zero || w_is_one) {
        int64_t pick = w_is_zero ? j_idx : upper_idx;
        expr_free(g_weight);
        return expr_copy(stats_order_at(sorted, pick - 1));
    }

    /* Two algebraically identical forms, each with a floating-point failure the
     * other does not have. INSIDE the unit interval use the convex combination
     * (1-w) A[j] + w A[j+1]: it never forms a quantity larger in magnitude than
     * the two neighbours, whereas A[j] + w (A[j+1] - A[j]) makes the difference
     * first, and that difference overflows to Infinity for two neighbours of
     * opposite sign near the double range -- Infinity then survives the multiply
     * and the add, so Quantile[{-1.0*10^308, 1.0*10^308}, 1/2, {{1/2,0},{0,1}}]
     * returned Infinity where the answer is 0. Same class as the w == 1 case
     * handled just above: an identity in the reals, not in floating point.
     *
     * OUTSIDE it, keep the difference form. w = c + d g is user-controlled
     * through {{a,b},{c,d}} and nothing constrains it to [0,1]; every standard
     * Hyndman-Fan type lands inside, but for w outside, (1-w) and w have
     * opposite signs and the two products can overflow independently, giving
     * inf + -inf = NaN on input the difference form handles exactly (equal huge
     * neighbours, where A + w*0 is simply A). Neither form is safe everywhere;
     * each is used where it is safe. A weight that does not reduce to a machine
     * number keeps the historical form. */
    double w_val = 0.0;
    bool w_known = stats_is_numeric(g_weight, &w_val, NULL);
    if (!w_known) {
        Expr* nw = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
                       (Expr*[]){expr_copy(g_weight)}, 1));
        w_known = stats_is_numeric(nw, &w_val, NULL);
        expr_free(nw);
    }

    if (w_known && w_val >= 0.0 && w_val <= 1.0) {
        Expr* one_minus_w = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_new_integer(1),
                       eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){expr_new_integer(-1), expr_copy(g_weight)}, 2)) }, 2));

        Expr* lo_term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){one_minus_w, expr_copy(stats_order_at(sorted, j_idx-1))}, 2));
        Expr* hi_term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){g_weight, expr_copy(stats_order_at(sorted, upper_idx-1))}, 2));

        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){lo_term, hi_term}, 2));
    }

    Expr* neg_Aj1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(stats_order_at(sorted, j_idx-1))}, 2));
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(stats_order_at(sorted, upper_idx-1)), neg_Aj1}, 2));

    Expr* weight_diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){g_weight, diff}, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(stats_order_at(sorted, j_idx-1)), weight_diff}, 2));
}

bool stats_machine_ok(double v) {
    return isfinite(v) && !(v != 0.0 && fabs(v) < DBL_MIN);
}

/* ------------------------------------------------------------------------
 * stats_sort_order: an index sort (bottom-up merge, stable) into NUMERIC order
 * (numeric_order.h) -- the order statistics are defined on the data sorted by
 * value, and canonical order is not that order once an exact non-rational is
 * present (Median[{Sqrt[2], 1, 3}] was 3). Numerically equal elements (1 and
 * 1.) fall back to canonical order, so the permutation is deterministic.
 * A typed LinkedList is compared on its raw slots, gathered into a machine
 * buffer; otherwise each element's key is computed once (O(1) for a number
 * atom) and the comparator reads the elements in place.
 * ---------------------------------------------------------------------- */
typedef struct { LLReader* rd; const int64_t* iv; const double* dv; const NumOrdKey* keys; } SortCtx;

static Expr* sort_elem(const SortCtx* c, size_t i, Expr** owned) {
    *owned = NULL;
    if (!c->rd->snap) return c->rd->v[i];
    Expr* b = (Expr*)ll_peek_boxed(c->rd->snap, (int64_t)i);
    if (b) return b;
    *owned = ll_get(c->rd->snap, (int64_t)i);       /* complex slots */
    return *owned;
}

static int sort_cmp(const SortCtx* c, size_t i, size_t j) {
    if (c->iv) return (c->iv[i] > c->iv[j]) - (c->iv[i] < c->iv[j]);
    if (c->dv) return (c->dv[i] > c->dv[j]) - (c->dv[i] < c->dv[j]);
    const NumOrdKey* ki = &c->keys[i];
    const NumOrdKey* kj = &c->keys[j];
    /* Two machine values: decided here, with no element read. */
    if (ki->exact && kj->exact && ki->v != kj->v) return ki->v < kj->v ? -1 : 1;
    Expr *oi, *oj;
    Expr* a = sort_elem(c, i, &oi);
    Expr* b = sort_elem(c, j, &oj);
    int r = numord_cmp(a, ki, b, kj);
    if (r == 0) r = expr_compare(a, b);
    if (oi) expr_free(oi);
    if (oj) expr_free(oj);
    return r;
}

size_t* stats_sort_order(LLReader* rd) {
    size_t n = rd->n;
    size_t* p = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* t = malloc(sizeof(size_t) * (n ? n : 1));
    if (!p || !t) { free(p); free(t); return NULL; }
    SortCtx c = { rd, NULL, NULL, NULL };
    void* raw = NULL;
    NumOrdKey* keys = NULL;
    if (rd->snap && n > 0) {
        LLKind k = ll_kind(rd->snap);
        if (k == LL_I64 || k == LL_REAL) {
            raw = malloc(8 * n);
            if (raw) {
                ll_gather(rd->snap, 0, (int64_t)n, raw);
                if (k == LL_I64) c.iv = raw; else c.dv = raw;
            }
        }
    }
    if (!c.iv && !c.dv) {
        keys = malloc(sizeof(NumOrdKey) * (n ? n : 1));
        if (!keys) { free(p); free(t); free(raw); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* o;
            numord_key(sort_elem(&c, i, &o), &keys[i]);
            if (o) expr_free(o);
        }
        c.keys = keys;
    }
    for (size_t i = 0; i < n; i++) p[i] = i;
    for (size_t w = 1; w < n; w *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * w) {
            size_t mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
            size_t a = lo, b = mid, o = lo;
            while (a < mid && b < hi) t[o++] = sort_cmp(&c, p[b], p[a]) < 0 ? p[b++] : p[a++];
            while (a < mid) t[o++] = p[a++];
            while (b < hi) t[o++] = p[b++];
        }
        size_t* x = p; p = t; t = x;
    }
    free(t);
    free(raw);
    free(keys);
    return p;
}
