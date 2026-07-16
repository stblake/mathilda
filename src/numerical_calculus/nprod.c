/*
 * nprod.c — NProduct[f, {i, imin, imax (, di)}, opts]   (see nprod.h)
 *
 * Strategy
 * --------
 * Per Keiper 1992 ("The N functions of Mathematica", G.15), a numerical product
 * is evaluated as the exponential of a numerical sum of logarithms:
 *
 *     Prod_{i=imin}^{imax} f(i)  =  Exp[ NSum[ Log[f(i)], {i, imin, imax} ] ].
 *
 * This is exact for finite ranges (Exp inverts Log on any branch and the
 * principal-log phases that wind by multiples of 2*pi are unwrapped by Exp) and
 * correct for convergent infinite products (factors -> 1 => Log f -> 0 is smooth
 * on the principal branch, and Euler-Maclaurin uses only branch-independent
 * derivatives f'/f).  We therefore delegate every hard part — method selection,
 * Euler-Maclaurin, Wynn epsilon, Cohen-Villegas-Zagier, MPFR working precision,
 * large-finite tail differences, and divergence detection — to the existing,
 * tested NSum engine, and only:
 *   - parse NProduct's own option names and map them onto NSum's;
 *   - add guard digits on the arbitrary-precision path (Exp amplifies the
 *     absolute error of the exponent into relative error of the product);
 *   - handle multidimensional products by recursion (inner NProduct as body);
 *   - special-case divergence (NSum -> ComplexInfinity).
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.
 */

#include "nprod.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* Extra digits carried on the MPFR path before the final Exp/round. */
#define NP_GUARD_DIGITS 10.0

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */

static void np_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NProduct::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

/* Iterator spec: List with 2..4 args whose first is a symbol (mirrors NSum). */
static bool np_is_spec(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    return n >= 2 && n <= 4 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

static bool np_is_known_option(const char* s) {
    return s == SYM_Method || s == SYM_WorkingPrecision
        || s == SYM_NProductFactors || s == SYM_NProductExtraFactors
        || s == SYM_WynnDegree || s == SYM_VerifyConvergence
        || s == SYM_AccuracyGoal || s == SYM_PrecisionGoal
        || s == SYM_Compiled || s == SYM_EvaluationMonitor;
}

static bool np_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && np_is_known_option(lhs->data.symbol.name);
}

/* WorkingPrecision values are small integer/real literals in practice. */
static bool np_to_double(Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    return false;
}

static bool np_is_numberlike(Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: {
            Expr* h = e->data.function.head;
            return h->type == EXPR_SYMBOL
                && (h->data.symbol.name == SYM_Complex || h->data.symbol.name == SYM_Rational);
        }
        default:
            return false;
    }
}

static bool np_is_infinity(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_ComplexInfinity
            || e->data.symbol.name == SYM_Indeterminate
            || e->data.symbol.name == SYM_Infinity;
    if (e->type == EXPR_FUNCTION) {
        Expr* h = e->data.function.head;
        return h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_DirectedInfinity;
    }
    return false;
}

/* Build Rule[Symbol(name), rhs] adopting rhs. */
static Expr* np_rule(const char* name, Expr* rhs) {
    Expr* a[2];
    a[0] = expr_new_symbol(name);
    a[1] = rhs;
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* Wrap head[child] adopting child. */
static Expr* np_apply1(const char* head, Expr* child) {
    Expr* a[1];
    a[0] = child;
    return expr_new_function(expr_new_symbol(head), a, 1);
}

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr* factors;       /* NProductFactors rhs (borrowed) or NULL      */
    Expr* extra_factors; /* NProductExtraFactors rhs or NULL            */
    Expr* method;        /* Method rhs or NULL                          */
    Expr* wynn;          /* WynnDegree rhs or NULL                      */
    Expr* verify;        /* VerifyConvergence rhs or NULL               */
    Expr* wprec;         /* WorkingPrecision rhs or NULL                */
    Expr* accgoal;       /* AccuracyGoal rhs or NULL                    */
    Expr* precgoal;      /* PrecisionGoal rhs or NULL                   */
    bool   mpfr;         /* WorkingPrecision selects arbitrary precision */
    double wdigits;      /* requested working-precision digits (mpfr)   */
} NpOpts;

static bool np_apply_option(Expr* rule, NpOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_NProductFactors)      { o->factors = rhs;       return true; }
    if (name == SYM_NProductExtraFactors) { o->extra_factors = rhs; return true; }
    if (name == SYM_Method)               { o->method = rhs;        return true; }
    if (name == SYM_WynnDegree)           { o->wynn = rhs;          return true; }
    if (name == SYM_VerifyConvergence)    { o->verify = rhs;        return true; }
    if (name == SYM_AccuracyGoal)         { o->accgoal = rhs;       return true; }
    if (name == SYM_PrecisionGoal)        { o->precgoal = rhs;      return true; }
    if (name == SYM_WorkingPrecision) {
        o->wprec = rhs;
        o->mpfr = false; o->wdigits = 0.0;
        if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_MachinePrecision)) {
            double d;
            if (np_to_double(rhs, &d) && d > NUMERIC_MACHINE_PRECISION_DIGITS) {
                o->mpfr = true; o->wdigits = d;
            }
        }
        return true;
    }
    /* Compiled / EvaluationMonitor: accepted, ignored. */
    if (name == SYM_Compiled || name == SYM_EvaluationMonitor) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nprod(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    Expr** args = res->data.function.args;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && np_is_option_arg(args[pos_end - 1])) pos_end--;
    if (pos_end < 2) return NULL;                       /* need f + >=1 spec */

    Expr* body = args[0];
    for (size_t i = 1; i < pos_end; i++)
        if (!np_is_spec(args[i])) return NULL;

    NpOpts o; memset(&o, 0, sizeof o);
    for (size_t i = pos_end; i < argc; i++)
        if (!np_apply_option(args[i], &o)) return NULL;

    size_t nspecs = pos_end - 1;
    Expr* spec0 = args[1];

    /* Effective body of the outer product: for >1 spec it is an inner NProduct
     * over the remaining specs (HoldAll + index localisation lets a dependent
     * inner bound such as {k,1,n} see the outer index). The inner NProduct
     * receives the ORIGINAL options (it maps them itself). */
    Expr* eff_body;
    if (nspecs >= 2) {
        size_t n = 1 /*body*/ + (nspecs - 1) /*specs[2..]*/ + (argc - pos_end) /*opts*/;
        Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!v) return NULL;
        size_t w = 0;
        v[w++] = expr_copy(body);
        for (size_t i = 2; i < pos_end; i++) v[w++] = expr_copy(args[i]);
        for (size_t i = pos_end; i < argc; i++) v[w++] = expr_copy(args[i]);
        eff_body = expr_new_function(expr_new_symbol(SYM_NProduct), v, w);
        free(v);
    } else {
        eff_body = expr_copy(body);
    }

    /* Build the inner NSum[Log[eff_body], spec0, <mapped options>]. */
    Expr* logbody = np_apply1("Log", eff_body);

    Expr** nv = malloc(sizeof(Expr*) * 12);
    if (!nv) { expr_free(logbody); return NULL; }
    size_t nw = 0;
    nv[nw++] = logbody;                 /* adopts eff_body                  */
    nv[nw++] = expr_copy(spec0);
    if (o.factors)       nv[nw++] = np_rule("NSumTerms",      expr_copy(o.factors));
    if (o.extra_factors) nv[nw++] = np_rule("NSumExtraTerms", expr_copy(o.extra_factors));
    if (o.method)        nv[nw++] = np_rule("Method",         expr_copy(o.method));
    if (o.wynn)          nv[nw++] = np_rule("WynnDegree",     expr_copy(o.wynn));
    if (o.verify)        nv[nw++] = np_rule("VerifyConvergence", expr_copy(o.verify));
    if (o.accgoal)       nv[nw++] = np_rule("AccuracyGoal",   expr_copy(o.accgoal));
    if (o.precgoal)      nv[nw++] = np_rule("PrecisionGoal",  expr_copy(o.precgoal));
    if (o.mpfr) {
        long inner = (long)(o.wdigits + NP_GUARD_DIGITS + 0.5);
        nv[nw++] = np_rule("WorkingPrecision", expr_new_integer(inner));
    } else if (o.wprec) {
        nv[nw++] = np_rule("WorkingPrecision", expr_copy(o.wprec));
    }

    Expr* nsum = expr_new_function(expr_new_symbol(SYM_NSum), nv, nw);
    free(nv);

    Expr* s = eval_and_free(nsum);      /* the numeric Log-sum, or special  */

    /* Divergent product: NSum reports the divergent log-sum as ComplexInfinity. */
    if (s && s->type == EXPR_SYMBOL && s->data.symbol.name == SYM_ComplexInfinity) {
        expr_free(s);
        return expr_new_symbol(SYM_ComplexInfinity);
    }
    /* Could not reduce to a number (symbolic body / bound): stay unevaluated. */
    if (!s || (!np_is_numberlike(s) && !np_is_infinity(s))) {
        expr_free(s);
        return NULL;
    }

    Expr* result = eval_and_free(np_apply1("Exp", s));   /* s -> Exp node    */

    /* Round back to the requested working precision on the MPFR path. */
    if (o.mpfr && np_is_numberlike(result)) {
        Expr* a[2];
        a[0] = result;
        a[1] = expr_new_integer((long)(o.wdigits + 0.5));
        result = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), a, 2));
    }

    (void)np_warn;   /* reserved for future divergence/convergence messages  */
    return result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nprod_init(void) {
    symtab_add_builtin("NProduct", builtin_nprod);
    /* HoldAll: the factor and the iterator specs must not be pre-evaluated; the
     * index is Block-localised inside the delegated NSum. Not Listable. */
    symtab_get_def("NProduct")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
