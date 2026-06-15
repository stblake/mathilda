/* Mathilda -- the logarithmic integral li.
 *
 *   LogIntegral[z]   li(z) = PV Int_0^z dt / ln t
 *
 * li has a branch cut along (-Infinity, +1); the principal value is taken on
 * the cut. The implementation rests on the identity
 *
 *   li(z) = Ei(Log z),
 *
 * where Ei is ExpIntegralEi and Log is the principal logarithm. This lets us
 * reuse ExpIntegralEi's fully-tested numeric stack (mpfr_eint / the real and
 * complex convergent series with cancellation guard bits) without duplicating
 * any of it, and the principal Log automatically supplies the +-i Pi jump that
 * places the branch cut on (-Infinity, +1).
 *
 * Evaluation is layered so each kind of argument takes the cheapest route:
 *
 *   exact special values     ->  0, -Infinity, Infinity, Indeterminate
 *   numeric (inexact) z       ->  evaluate ExpIntegralEi[Log[z]]
 *   everything else           ->  stays symbolic (return NULL)
 *
 * Exact non-special numbers (e.g. LogIntegral[2], LogIntegral[1/2]) stay
 * symbolic, matching the Wolfram Language; only inexact input or an explicit
 * N[...] (which rewrites the argument to an MPFR number) evaluates numerically.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "logintegral.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "arithmetic.h"   /* is_complex, make_complex */
#include "attr.h"
#include "eval.h"         /* eval_and_free */
#include "expr.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool li_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` is a function with head symbol `name`. */
static bool li_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool li_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries inexact numeric content: a Real/MPFR leaf, or a
 * Complex[..] with at least one inexact part. These are the arguments for
 * which li(z) is evaluated numerically. */
static bool li_is_numeric_inexact(Expr* e) {
    if (li_is_inexact(e)) return true;
    Expr *re, *im;
    if (is_complex(e, &re, &im) && (li_is_inexact(re) || li_is_inexact(im)))
        return true;
    return false;
}

/* Build -Infinity = Times[-1, Infinity]. */
static Expr* li_make_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* ------------------------------------------------------------------ */
/* LogIntegral[z]                                                     */
/* ------------------------------------------------------------------ */

static Expr* li_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(0);                          /* li(0) = 0      */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 1)
        return li_make_neg_infinity();                       /* li(1) = -Inf   */
    if (li_is_symbol(arg, "Infinity"))         return expr_new_symbol(SYM_Infinity);
    if (li_is_symbol(arg, "ComplexInfinity"))  return expr_new_symbol(SYM_Indeterminate);
    if (li_is_symbol(arg, "Indeterminate"))    return expr_new_symbol(SYM_Indeterminate);

    /* 2. Numeric (inexact) arguments: li(z) = Ei(Log z). Build
     * ExpIntegralEi[Log[z]] and evaluate it through the existing kernels. */
    if (li_is_numeric_inexact(arg)) {
        Expr* inner = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){ expr_copy(arg) }, 1);
        Expr* call = expr_new_function(expr_new_symbol(SYM_ExpIntegralEi),
            (Expr*[]){ inner }, 1);
        Expr* out = eval_and_free(call);
        /* Defensive: only accept a reduced (numeric) result. If the
         * composition failed to evaluate (head still ExpIntegralEi/Log), drop
         * it and stay symbolic rather than leak the wrong head. */
        if (out && !li_head_is(out, "ExpIntegralEi") && !li_head_is(out, "Log"))
            return out;
        expr_free(out);
        return NULL;
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* li_emit_argx(size_t argc) {
    fprintf(stderr,
            "LogIntegral::argx: LogIntegral called with %zu arguments; "
            "1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_logintegral(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return li_one_arg(args[0]);
    return li_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void logintegral_init(void) {
    symtab_add_builtin("LogIntegral", builtin_logintegral);
    symtab_get_def("LogIntegral")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
