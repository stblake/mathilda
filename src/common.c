/* Mathilda — common preprocessing helpers.  See common.h. */

#include "common.h"

#include "numeric.h"
#include "rationalize.h"
#include "sym_names.h"   /* SYM_Rule / SYM_RuleDelayed / SYM_Method */
#include "message.h"
#include "eval.h"      /* eval_and_free: alg_thread_call */
#include "print.h"     /* expr_to_string for the ::normal message */     /* mth_message: route arg-count errors through Quiet/Check */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

bool head_is(const Expr* e, const char* sym) {
    return e && sym &&
           e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

Expr* builtin_arg_error(const char* head, size_t got, size_t min, size_t max) {
    const char* got_s = (got == 1) ? "" : "s"; /* "argument" vs "arguments" */

    if (min == max) {
        /* Fixed arity. The reference tags an expected count of one `argx`, and
         * a call with ONE argument against a larger fixed count `argr`. */
        const char* tag = (min == 1) ? "argx" : (got == 1) ? "argr" : "argrx";
        mth_message(head, tag,
                "%s called with %zu argument%s; %zu argument%s %s expected.",
                head, got, got_s,
                min, (min == 1) ? "" : "s", (min == 1) ? "is" : "are");
    } else if (max == SIZE_MAX) {
        mth_message(head, "argm",
                "%s called with %zu argument%s; %zu or more arguments are expected.",
                head, got, got_s, min);
    } else if (max == min + 1) {
        mth_message(head, (got == 1) ? "argtu" : "argt",
                "%s called with %zu argument%s; %zu or %zu arguments are expected.",
                head, got, got_s, min, max);
    } else {
        mth_message(head, "argb",
                "%s called with %zu argument%s; between %zu and %zu arguments are expected.",
                head, got, got_s, min, max);
    }
    return NULL;
}

bool alg_thread_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL)
        return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_Rule || h == SYM_Equal || h == SYM_Unequal || h == SYM_Less ||
           h == SYM_Greater || h == SYM_LessEqual || h == SYM_GreaterEqual ||
           h == SYM_Inequality || h == SYM_And || h == SYM_Or || h == SYM_Not ||
           h == SYM_Xor || h == SYM_Implies;
}

Expr* alg_thread_call(const Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    const Expr* e = res->data.function.args[0];
    if (!alg_thread_head(e)) return NULL;
    size_t n = e->data.function.arg_count, argc = res->data.function.arg_count;
    bool ineq = e->data.function.head->data.symbol.name == SYM_Inequality;
    /* Threading is transparent, as Listable threading is: the inner calls run
     * at the caller's builtin depth, so a user-level Together[x -> ...] still
     * answers in its user-level form.  (The counter is registered with the
     * TimeConstrained unwind, which restores it on an abort.) */
    int* depth = eval_builtin_depth_counter();
    int saved = *depth;
    if (*depth > 0) (*depth)--;
    Expr** na = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        if (ineq && (i & 1u)) { na[i] = expr_copy(e->data.function.args[i]); continue; }
        Expr** fa = malloc(sizeof(Expr*) * argc);
        fa[0] = expr_copy(e->data.function.args[i]);
        for (size_t k = 1; k < argc; k++) fa[k] = expr_copy(res->data.function.args[k]);
        na[i] = eval_and_free(expr_new_function(expr_copy(res->data.function.head), fa, argc));
        free(fa);
    }
    Expr* out = expr_new_function(expr_copy(e->data.function.head), na, n);
    free(na);
    out = eval_and_free(out);
    *depth = saved;
    return out;
}

Expr* builtin_normal_error(const char* head, const Expr* call, int pos) {
    char* shown = expr_to_string((Expr*)call);
    mth_message(head, "normal", "Nonatomic expression expected at position %d in %s.",
                pos, shown ? shown : "");
    free(shown);
    return NULL;
}

/* Bit precision of a single inexact leaf.  Real → 53 (IEEE 754 double).
 * MPFR → mpfr_get_prec().  Caller is expected to only invoke this on
 * leaves where is_inexact_leaf would have returned true. */
static long leaf_precision_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_REAL) return 53;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
#endif
    return 0;
}

static bool is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Recursive walk: update *info in place with min precision and has_inexact
 * flag.  We track the running minimum so a single pass captures both. */
static void scan_inexact_walk(const Expr* e, CommonInexactInfo* info) {
    if (!e) return;
    if (is_inexact_leaf(e)) {
        long bits = leaf_precision_bits(e);
        if (!info->has_inexact || bits < info->min_bits) {
            info->min_bits = bits;
        }
        info->has_inexact = true;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    scan_inexact_walk(e->data.function.head, info);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        scan_inexact_walk(e->data.function.args[i], info);
    }
}

CommonInexactInfo common_scan_inexact(const Expr* e) {
    CommonInexactInfo info = { false, 0 };
    scan_inexact_walk(e, &info);
    return info;
}

Expr* common_rationalize_input(const Expr* e, long bits) {
    return internal_force_rationalize_bits(e, bits);
}

bool common_has_machine_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_NDARRAY)
        return e->data.ndarray.dtype != NDT_INT64;
    if (e->type == EXPR_FUNCTION) {
        if (common_has_machine_real(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (common_has_machine_real(e->data.function.args[i])) return true;
    }
    return false;
}

bool common_machine_real_value(const Expr* e, double* out) {
    if (!e || !out) return false;
    switch (e->type) {
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
        case EXPR_FUNCTION: {
            /* Rational[p, q]. Deliberately NOT a general numeric reader: MPFR
             * and Complex are excluded because the contagion rule this serves
             * is about MACHINE reals, and neither is one. An MPFR diagonal
             * keeps its exact zeros in Mathematica, so it must decline here. */
            const Expr* h = e->data.function.head;
            if (h && h->type == EXPR_SYMBOL &&
                h->data.symbol.name == SYM_Rational &&
                e->data.function.arg_count == 2) {
                double p, q;
                if (common_machine_real_value(e->data.function.args[0], &p) &&
                    common_machine_real_value(e->data.function.args[1], &q) &&
                    q != 0.0) {
                    *out = p / q;
                    return true;
                }
            }
            return false;
        }
        default: return false;
    }
}

Expr* common_numericalize_result(const Expr* e, long bits) {
    if (!e) return NULL;
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (bits > 53) {
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = bits;
    }
#else
    (void)bits;
#endif
    return numericalize(e, spec);
}

/* Is `e` a `Method -> ...` (or `Method :> ...`) option rule? */
static bool is_method_option(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!head_is(e, SYM_Rule) && !head_is(e, SYM_RuleDelayed)) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    return lhs && lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Method;
}

Expr* common_method_alias(Expr* res, const char* head, size_t n_positional,
                          const char* method, Expr* (*impl)(Expr*)) {
    if (!res || res->type != EXPR_FUNCTION || !head || !method || !impl) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < n_positional) return NULL;

    /* One extra slot for the Method option we append; dropping the caller's
     * Method options only shortens the list. */
    Expr** args = calloc(argc + 1, sizeof(Expr*));
    if (!args) return NULL;

    size_t n = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (i >= n_positional && is_method_option(a)) continue;
        args[n++] = expr_copy(a);
    }

    /* expr_new_function copies the argument array into its own storage, so a
     * stack array is enough here (the Expr* elements themselves are adopted). */
    Expr* opt_args[2];
    opt_args[0] = expr_new_symbol("Method");
    opt_args[1] = expr_new_string(method);
    args[n++] = expr_new_function(expr_new_symbol("Rule"), opt_args, 2);

    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    free(args);

    /* Direct call, so the ownership contract makes us the evaluator here. */
    Expr* out = impl(call);
    expr_free(call);
    return out;
}
