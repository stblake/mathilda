#include "list_common.h"
#include "range.h"
#include "../pack.h"
#include "../checked_int.h"
#include "../ndarray.h"   /* ndarray_warn_once */

/* Range[] is the system's most-used list producer, so it is the first place
 * automatic packing pays for itself. Both branches below now decide the element
 * COUNT before producing anything, then fill either a dense buffer (packed) or
 * an Expr array (plain) from that one count -- so the two representations agree
 * by construction, which is the property the differential suite checks.
 *
 * Deciding the count up front also fixes a bug on the exact branch: it used to
 * drive the loop off the `double` shadow of the bounds, and one ulp at 10^18 is
 * 128, so `val += 1` never advanced and `Range[10^18, 10^18 + 3]` ran to the
 * safety cap instead of giving four elements. The exact branch now counts in
 * int64. */

/* How much a Range result may occupy before the call is refused.
 *
 * This replaces a 10^6-ELEMENT cap that TRUNCATED SILENTLY: `Range[2000000]`
 * answered with 1000001 elements and said nothing, so every downstream length,
 * total and plot was quietly computed on half the data. A resource limit is
 * legitimate; a resource limit that changes the answer instead of declining is
 * not, and one that does so without a message cannot even be noticed. Found on
 * 2026-08-02 while sizing a benchmark vector, which is exactly how a silent
 * truncation gets found — by accident, downstream of the damage.
 *
 * Two things changed. The ceiling is now on BYTES rather than elements, which
 * is what actually runs out, and it is ~268 million int64 elements instead of
 * 10^6. And reaching it is REPORTED and leaves the call unevaluated, so the
 * caller sees a Range[...] it can act on rather than a plausible short list. */
#define RANGE_MAX_BYTES ((int64_t)1 << 31)          /* 2 GiB of result */
#define RANGE_MAX_ELEMENTS (RANGE_MAX_BYTES / (int64_t)sizeof(int64_t))

/* The ceiling is deliberately the SAME for the packed and the boxed paths,
 * even though a boxed element costs an order of magnitude more than the eight
 * bytes it is sized from. Giving the boxed path a tighter one would make
 * MATHILDA_NO_PACK=1 refuse a Range that packing accepts, and the two
 * representations agreeing on what is admissible is worth more than a
 * better-calibrated number. The boxed paths check their allocation instead and
 * report the same refusal if it fails. */

/* Report the refusal once and let the caller return the call unevaluated. */
static void range_too_large(double count) {
    char msg[192];
    snprintf(msg, sizeof msg,
             "Range::toobig: Range of about %.0f elements exceeds the "
             "%lld-element limit; the expression is left unevaluated.\n",
             count, (long long)RANGE_MAX_ELEMENTS);
    /* Once per evaluation: returning NULL leaves the node unevaluated and the
     * fixed-point loop comes back to it, so a bare fprintf printed this five
     * times for one Range. */
    ndarray_warn_once(msg);
}

/* The inexact recurrence, in one place so the counting pass and the filling
 * pass cannot drift apart. Accumulating `val += di` rather than computing
 * `min + k*di` is deliberate: it is what this branch has always emitted, and
 * the two differ in the last bits (Range[0., 1., 0.1] gives
 * 0.30000000000000004 for the third element either way, but not for every
 * step count). */
#define RANGE_MORE(val, max_val, di_val) \
    (((di_val) > 0 && (val) <= (max_val) + 1e-14) || \
     ((di_val) < 0 && (val) >= (max_val) - 1e-14))

Expr* builtin_range(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    size_t len = res->data.function.arg_count;
    Expr* imin_e = NULL;
    Expr* imax_e = NULL;
    Expr* di_e = NULL;

    if (len == 1) {
        imin_e = expr_new_integer(1);
        imax_e = expr_copy(res->data.function.args[0]);
        di_e = expr_new_integer(1);
    } else if (len == 2) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_new_integer(1);
    } else if (len == 3) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_copy(res->data.function.args[2]);
    }

    bool is_real = false;
    double min_val = 0, max_val = 0, di_val = 0;
    int64_t n, d;

    if (imin_e->type == EXPR_REAL || imax_e->type == EXPR_REAL || di_e->type == EXPR_REAL) is_real = true;

    /* All three bounds machine integers: the count, and every element, are
     * computable exactly. Taken before the double coercion below so values past
     * 2^53 are handled at full precision. */
    bool all_int = (imin_e->type == EXPR_INTEGER &&
                    imax_e->type == EXPR_INTEGER &&
                    di_e->type == EXPR_INTEGER);

    if (imin_e->type == EXPR_INTEGER) min_val = (double)imin_e->data.integer;
    else if (imin_e->type == EXPR_REAL) min_val = imin_e->data.real;
    else if (is_rational(imin_e, &n, &d)) min_val = (double)n / d;
    else goto L_fail_range;

    if (imax_e->type == EXPR_INTEGER) max_val = (double)imax_e->data.integer;
    else if (imax_e->type == EXPR_REAL) max_val = imax_e->data.real;
    else if (is_rational(imax_e, &n, &d)) max_val = (double)n / d;
    else goto L_fail_range;

    if (di_e->type == EXPR_INTEGER) di_val = (double)di_e->data.integer;
    else if (di_e->type == EXPR_REAL) di_val = di_e->data.real;
    else if (is_rational(di_e, &n, &d)) di_val = (double)n / d;
    else goto L_fail_range;

    if (all_int) {
        int64_t a = imin_e->data.integer;
        int64_t b = imax_e->data.integer;
        int64_t s = di_e->data.integer;
        int64_t span;
        /* A span that does not fit (a near INT64_MIN, b near INT64_MAX) means a
         * count far past the cap anyway; drop through to the generic path, which
         * refuses it the same way it always has. */
        if (s != 0 && !ci_sub_i64(b, a, &span)) {
            if ((s > 0 && span < 0) || (s < 0 && span > 0)) {
                expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
                return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            }
            /* span and s share a sign, so C99's truncation toward zero is
             * floor() here. */
            int64_t count = span / s + 1;
            if (count > RANGE_MAX_ELEMENTS) {
                range_too_large((double)count);
                goto L_fail_range;          /* unevaluated, not truncated */
            }

            expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
            imin_e = imax_e = di_e = NULL;

            /* Every element lies between a and b, so `cur += s` cannot leave
             * int64 -- no overflow check needed inside either loop. */
            int64_t* buf = NULL;
            Expr* packed = ndbuild_open_i64(count, &buf);
            if (packed) {
                int64_t cur = a;
                for (int64_t i = 0; i < count; i++) { buf[i] = cur; cur += s; }
                return packed;
            }
            Expr** results = malloc(sizeof(Expr*) * (size_t)count);
            if (!results) { range_too_large((double)count); return NULL; }
            int64_t cur = a;
            for (int64_t i = 0; i < count; i++) {
                results[i] = expr_new_integer(cur);
                cur += s;
            }
            Expr* out = expr_new_function(expr_new_symbol(SYM_List), results, (size_t)count);
            free(results);
            return out;
        }
    }

    if (di_val == 0) goto L_fail_range;
    if ((di_val > 0 && min_val > max_val) || (di_val < 0 && min_val < max_val)) {
        expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    /* How many elements this can possibly be, arithmetically, BEFORE the
     * counting loop runs. The loop below accumulates `val += di` one step at a
     * time -- deliberately, so the count agrees with the fill bit for bit --
     * and `Range[0., 1., 10^-300]` would spin in it for the rest of the
     * universe. The estimate bounds that work without deciding the count: a
     * range that cannot fit is refused here, and everything else is counted
     * exactly the way it always was. */
    {
        double est = (max_val - min_val) / di_val + 1.0;
        if (!(est >= 0.0) || est > (double)RANGE_MAX_ELEMENTS) {
            range_too_large(est);
            goto L_fail_range;
        }
    }

    /* Counting pass. Doubles only -- no allocation, no evaluation. */
    size_t count = 0;
    {
        double val = min_val;
        while (RANGE_MORE(val, max_val, di_val)) {
            count++;
            val += di_val;
            /* The estimate above is arithmetic and this loop accumulates, so
             * rounding can carry it one or two past. Stopping a little beyond
             * the ceiling keeps the two consistent without ever truncating a
             * range the ceiling admitted. */
            if ((int64_t)count > RANGE_MAX_ELEMENTS + 2) {
                range_too_large((double)count);
                goto L_fail_range;
            }
        }
    }

    if (is_real) {
        double* buf = NULL;
        Expr* packed = ndbuild_open_f64((int64_t)count, &buf);
        if (packed) {
            double val = min_val;
            for (size_t i = 0; i < count; i++) { buf[i] = val; val += di_val; }
            expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
            return packed;
        }
    }

    Expr** results = malloc(sizeof(Expr*) * (count ? count : 1));
    if (!results) { range_too_large((double)count); goto L_fail_range; }
    size_t results_count = 0;

    double val = min_val;
    /* `curr_e` is the EXACT running value, advanced by a real Plus so that a
     * Rational range stays exact and an Integer one promotes to a BigInt when it
     * has to. The inexact branch does not need it at all -- it emits from `val`
     * -- and advancing it there would cost a full evaluate() per element for a
     * result that is never read. src/list/table.c makes the same distinction
     * ("advance the exact running value (unused when compiling)"). */
    Expr* curr_e = is_real ? NULL : expr_copy(imin_e);

    while (results_count < count) {
        results[results_count++] = is_real ? expr_new_real(val) : expr_copy(curr_e);

        if (!is_real) {
            Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
            Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
            Expr* next_e = evaluate(next_expr);
            expr_free(next_expr);
            expr_free(curr_e);
            curr_e = next_e;
        }

        val += di_val;
    }

    if (curr_e) expr_free(curr_e);
    expr_free(imin_e);
    expr_free(imax_e);
    expr_free(di_e);

    Expr* result_list = expr_new_function(expr_new_symbol(SYM_List), results, results_count);
    free(results);
    return result_list;

L_fail_range:
    if (imin_e) expr_free(imin_e);
    if (imax_e) expr_free(imax_e);
    if (di_e) expr_free(di_e);
    return NULL;
}
