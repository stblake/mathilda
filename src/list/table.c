#include "list_common.h"
#include "table.h"
#include "compile/autocompile.h"
#include "iter.h"
#include "../pack.h"

/* Element backstop for a finite Table range, matching Sum's SUM_MAX_FINITE_TERMS
 * and Product's PRODUCT_MAX_FINITE_TERMS. Reaching it returns Table[...]
 * unevaluated rather than silently truncating (the old 1000000 cap did the
 * latter, and at a hundredth of this size — so Table[i,{i,1,2000000}] came back
 * with 1000001 elements). With the exact termination test in iter_range_continue
 * a well-formed finite range terminates on its own; this only guards a range
 * whose element count genuinely exceeds the limit. */
#define TABLE_MAX_FINITE_ELEMENTS 100000000LL

Expr* builtin_table(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    
    if (res->data.function.arg_count > 2) {
        Expr** inner_args = malloc(sizeof(Expr*) * 2);
        inner_args[0] = expr_copy(res->data.function.args[0]);
        inner_args[1] = expr_copy(res->data.function.args[res->data.function.arg_count - 1]);
        Expr* inner_table = expr_new_function(expr_new_symbol(SYM_Table), inner_args, 2);
        free(inner_args);
        
        Expr** outer_args = malloc(sizeof(Expr*) * (res->data.function.arg_count - 1));
        outer_args[0] = inner_table;
        for (size_t i = 1; i < res->data.function.arg_count - 1; i++) {
            outer_args[i] = expr_copy(res->data.function.args[i]);
        }
        Expr* outer_table = expr_new_function(expr_new_symbol(SYM_Table), outer_args, res->data.function.arg_count - 1);
        free(outer_args);
        
        Expr* eval_outer = evaluate(outer_table);
        expr_free(outer_table);
        return eval_outer;
    }
    
    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* ---- Parse the iterator spec (shared helper) ---- */
    IterSpec s;
    if (!iter_spec_parse(spec, &s)) return NULL;

    int is_n_times   = (s.kind == ITER_KIND_COUNT);
    int is_list_iter = (s.kind == ITER_KIND_LIST);
    double min_val = 0, max_val = 0, di_val = 0;
    bool is_real = false, is_inf = false;

    /* Table does not iterate to Infinity; allow_inf = false. */
    if (!is_list_iter) {
        if (!iter_spec_resolve_numeric(&s, /*allow_inf=*/false,
                                       &min_val, &max_val, &di_val,
                                       &is_real, &is_inf)) {
            iter_spec_free(&s);
            return NULL;
        }
    }

    /* Convenience aliases into the owned IterSpec (freed via iter_spec_free). */
    Expr* var_sym = s.var;
    Expr* imin_e  = s.imin;
    Expr* imax_e  = s.imax;
    Expr* di_e    = s.di;
    Expr* list_e  = s.list;

    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);
    bool overflow = false;   /* element count exceeded TABLE_MAX_FINITE_ELEMENTS */

    /* Exact-integer range longer than the backstop: decline up front, before
     * allocating a single element. The in-loop guard below (and Sum/Product)
     * only discover the overflow mid-fill; for an integer range the length is
     * O(1) to compute, so Table[i, {i, 1, 10^9}] returns unevaluated at once
     * instead of allocating a hundred million Exprs and then giving up. Real and
     * rational ranges fall through to the in-loop backstop. */
    if (!is_n_times && !is_list_iter && !is_real
        && expr_is_integer_like(imin_e) && expr_is_integer_like(imax_e)
        && di_e && expr_is_integer_like(di_e)) {
        mpz_t lo, hi, step, span;
        mpz_inits(lo, hi, step, span, NULL);
        expr_to_mpz(imin_e, lo);
        expr_to_mpz(imax_e, hi);
        expr_to_mpz(di_e, step);
        mpz_sub(span, hi, lo);                    /* imax - imin */
        bool decline = false;
        if (mpz_sgn(span) != 0 && mpz_sgn(span) == mpz_sgn(step)) {
            mpz_abs(span, span);
            mpz_abs(step, step);
            mpz_fdiv_q(span, span, step);         /* q = floor(|span|/|step|); count = q + 1 */
            decline = mpz_cmp_si(span, (long)TABLE_MAX_FINITE_ELEMENTS) >= 0;  /* count > cap */
        }
        mpz_clears(lo, hi, step, span, NULL);
        if (decline) { free(results); iter_spec_free(&s); return NULL; }
    }

    Rule* old_own = iter_spec_shadow(var_sym);

    if (is_n_times) {
        int64_t n = imax_e->data.integer;
        for (int64_t i = 0; i < n; i++) {
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else if (is_list_iter) {
        for (size_t i = 0; i < list_e->data.function.arg_count; i++) {
            symtab_add_own_value(var_sym->data.symbol.name, var_sym, list_e->data.function.args[i]);
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else {
        /* Auto-compile fast path — ONLY for an inexact (machine-real) iterator,
         * where the interpreter already produces machine reals, so exactness is
         * never at stake.  Exact (Integer/BigInt/Rational) iterators keep the
         * pure-interpreter path below and their results are bit-for-bit
         * unchanged.  A non-finite / complex compiled result falls back to the
         * interpreter for that element (which yields the complex/singular value).
         * When compiling, the exact running value `curr_e` is not needed, so its
         * per-element evaluate(Plus[...]) advance is skipped too. */
        AutoCompiled* ac = is_real ? autocompile_new(expr, (const Expr* const*)&var_sym, 1) : NULL;
        double val = min_val;
        Expr* curr_e = ac ? NULL : expr_copy(imin_e);

        /* Packed fast path: a compiled body whose result type is CT_REAL emits
         * nothing but machine reals, so the whole list can be written straight
         * into a buffer. Counted first with the loop's own recurrence -- doubles
         * only, no evaluation -- so the packed and interpreted paths produce the
         * same number of elements from the same sequence of `val`s.
         *
         * compiled_eval_real still self-guards finiteness, so an element where
         * the interpreter would give ComplexInfinity or a complex value returns
         * false. That abandons the buffer, hands back the elements already
         * written, and the ordinary loop below picks up from exactly where it
         * stopped -- same fallback the boxed path has always had, per element. */
        if (ac && autocompiled_result_is_real(ac)) {
            size_t total = 0;
            double v = min_val;
            /* Real packed path: is_real is true here, so the count uses the same
             * double recurrence as the fill below (no curr_e). */
            while (iter_range_continue(is_real, /*is_inf=*/false, NULL, NULL,
                                       v, max_val, di_val)) {
                if ((int64_t)total >= TABLE_MAX_FINITE_ELEMENTS) { overflow = true; break; }
                total++;
                v += di_val;
            }
            double* buf = NULL;
            Expr* packed = overflow ? NULL : ndbuild_open_f64((int64_t)total, &buf);
            if (packed) {
                size_t i = 0;
                while (i < total && autocompiled_eval_real(ac, &val, &buf[i])) {
                    i++;
                    val += di_val;
                }
                if (i == total) {
                    autocompiled_free(ac);
                    iter_spec_restore(var_sym, old_own);
                    iter_spec_free(&s);
                    free(results);
                    return packed;
                }
                /* Element i is not machine-real: finish on the List path. */
                if (total > results_cap) {
                    results_cap = total;
                    results = realloc(results, sizeof(Expr*) * results_cap);
                }
                ndbuild_abandon(packed, i, results);
                results_count = i;
            }
        }

        while (iter_range_continue(is_real, /*is_inf=*/false, curr_e, imax_e,
                                   val, max_val, di_val)) {
            if ((int64_t)results_count >= TABLE_MAX_FINITE_ELEMENTS) { overflow = true; break; }
            /* Boxed (not _eval_real) so an integer-valued body stays an Integer
             * — the element type is user-visible in the returned list. */
            Expr* eval_expr = ac ? autocompiled_eval_boxed(ac, &val) : NULL;
            if (!eval_expr) {   /* interpreter path (and per-element fallback) */
                Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
                symtab_add_own_value(var_sym->data.symbol.name, var_sym, i_val);
                eval_expr = evaluate(expr);
                expr_free(i_val);
            }
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;

            if (!ac) {   /* advance the exact running value (unused when compiling) */
                Expr* next_e = iter_step_add(curr_e, di_e);
                if (!next_e) {
                    Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
                    Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
                    next_e = evaluate(next_expr);
                    expr_free(next_expr);
                }
                expr_free(curr_e);
                curr_e = next_e;
            }

            val += di_val;
        }
        if (curr_e) expr_free(curr_e);
        autocompiled_free(ac);
    }

    /* A finite range whose element count exceeded the backstop: free what was
     * collected and leave Table[...] unevaluated, as Sum/Product do on overflow,
     * rather than returning a silently-truncated list. */
    if (overflow) {
        for (size_t i = 0; i < results_count; i++) expr_free(results[i]);
        free(results);
        iter_spec_restore(var_sym, old_own);
        iter_spec_free(&s);
        return NULL;
    }

    iter_spec_restore(var_sym, old_own);
    iter_spec_free(&s);

    Expr* result_list = expr_new_function(expr_new_symbol(SYM_List), results, results_count);
    free(results);
    /* Catches every branch the direct path above does not: the exact-iterator
     * loop, Table[expr, {n}], Table[expr, {i, list}], and a body that compiled to
     * an integer or complex result. Declines in O(1) on anything symbolic. */
    return pack_offer(result_list);
}
