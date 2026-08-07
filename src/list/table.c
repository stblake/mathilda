#include "list_common.h"
#include "table.h"
#include "compile/autocompile.h"
#include "iter.h"
#include "../pack.h"

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
        int steps = 0;
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
            while ((di_val > 0 && v <= max_val + 1e-14) || (di_val < 0 && v >= max_val - 1e-14)) {
                total++;
                v += di_val;
                if (total > 1000000) break;
            }
            double* buf = NULL;
            Expr* packed = ndbuild_open_f64((int64_t)total, &buf);
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
                steps = (int)i;
            }
        }

        while ((di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
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
            steps++;
            if (steps > 1000000) break;
        }
        if (curr_e) expr_free(curr_e);
        autocompiled_free(ac);
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
