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
 *   - on an infinite range, sum a cancellation-free form of Log[f]
 *     (2 ArcTanh[g/(2 + g)], g = f - 1; see np_log_summand), since machine
 *     Log[1 + g] loses g's digits and vanishes in a far tail;
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
#include "iter.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "arithmetic.h" /* is_rational, is_complex */
#include "nsum.h"       /* nsum_eval_for_nproduct */
#include "pack.h"       /* pack_eval_plain */

/* Extra digits carried on the MPFR path before the final Exp/round. */
#define NP_GUARD_DIGITS 10.0

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */

static void np_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(0, "NProduct", tag, fmt, ap);
    va_end(ap);
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
 *  The log-summand of an infinite product                             *
 *                                                                     *
 *  A convergent infinite product has factors f = 1 + g with g -> 0,    *
 *  and Log[f] evaluated in floating point is Log of a number already   *
 *  rounded to 1 + g: it carries only the digits of g that survive the  *
 *  addition, |log f - g| / |g| ~ eps / |g|.  For 1 + 1/k^2 that is a    *
 *  relative error of 1e-4 at k = 1e6, and from k ~ 1e8 on Log[f] is an  *
 *  exact 0.  NSum's tail integral samples the summand far beyond the   *
 *  head (its power-law ladder starts eight decades out), so it saw a   *
 *  summand that stops dead: the tail model failed, Wynn took over and  *
 *  NProduct[1 + 1/k^2, {k, 10^6 + 1, Infinity}] came back 1 + 8e-11    *
 *  (true value 1 + 1e-6) with NSum::ncvg; even from k = 1 the answer   *
 *  was wrong in the ninth digit, with an NSum::accgl.                  *
 *                                                                     *
 *  The summand is therefore written without the cancellation,          *
 *      Log[1 + g] = 2 ArcTanh[g / (2 + g)],                            *
 *  with g = f - 1 formed symbolically (and put over a common           *
 *  denominator when the 1 does not cancel by itself, so               *
 *  (k^2 + 1)/k^2 - 1 becomes 1/k^2).  The identity is exact on the     *
 *  principal branch for every w = 1 + g off the cut (-oo, 0]: with     *
 *  z = (w - 1)/(w + 1), 1 + z = 2w/(w + 1) and 1 - z = 2/(w + 1), and   *
 *  for Im w != 0 or w > 0 the arguments arg w - arg(w + 1) and          *
 *  -arg(w + 1) lie in (-Pi, Pi), so Log(1 + z) - Log(1 - z) = Log w.    *
 *  On the cut the two differ by 2 Pi I, which Exp removes.  The        *
 *  derivatives the Euler-Maclaurin corrections need are taken of the   *
 *  equal Log[f] (NSum's derivative twin, nsum_eval_for_nproduct):      *
 *  f'/f is compact and cancellation-free, while those of the ArcTanh    *
 *  form balloon.  A factor E^h needs no Log at all: Exp[Sum h] is the   *
 *  product, whatever branch Log[E^h] would pick.                        *
 *                                                                     *
 *  The ArcTanh form is singular where Log[f] is (f = 0) and also at     *
 *  f = -1 (g = -2), which can only occur where |g| is not small, i.e.   *
 *  among the leading factors.  Those are probed until NP_PROBE_SETTLED *
 *  consecutive |g| < 1/2; when a -1 is among them, every factor before *
 *  the settled run is multiplied out directly and the sum starts after *
 *  it.  A product whose factors are not seen to settle within           *
 *  NP_PROBE_MAX keeps Log[f] (the former summand, exact at the head).   *
 * ------------------------------------------------------------------ */

#define NP_PROBE_CHUNK   32     /* factors evaluated per Table call          */
#define NP_PROBE_MAX     1024   /* leading factors probed for f = -1         */
#define NP_PROBE_SETTLED 16     /* consecutive |g| < 1/2 that end the probe  */

/* Does `e` mention the (interned) symbol `name`? */
static bool np_mentions(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type != EXPR_FUNCTION) return false;
    if (np_mentions(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (np_mentions(e->data.function.args[i], name)) return true;
    return false;
}

/* A sum with an addend free of the index: the signature of a 1 that did not
 * cancel in f - 1 (-1 + (1 + k^2)/k^2). */
static bool np_has_constant_addend(const Expr* g, const char* var) {
    if (!g || g->type != EXPR_FUNCTION || g->data.function.head->type != EXPR_SYMBOL
        || g->data.function.head->data.symbol.name != SYM_Plus) return false;
    for (size_t i = 0; i < g->data.function.arg_count; i++)
        if (!np_mentions(g->data.function.args[i], var)) return true;
    return false;
}

/* A machine number (integer, real, rational, or complex of those). */
static bool np_to_complex(Expr* e, double* re, double* im) {
    int64_t rn, rd;
    *im = 0.0;
    if (np_to_double(e, re)) return true;
    if (e && e->type == EXPR_FUNCTION && is_rational(e, &rn, &rd)) {
        *re = (double)rn / (double)rd; return true;
    }
    Expr *a, *b;
    if (e && is_complex(e, &a, &b)) {
        double x, y;
        if (np_to_complex(a, &x, &y) && np_to_complex(b, im, &y)) { *re = x; return true; }
    }
    return false;
}

/* head[a, b], adopting both. */
static Expr* np_fn2(const char* head, Expr* a, Expr* b) {
    Expr* v[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), v, 2);
}

/* The iterator pieces of `spec` (borrowed; imin/di NULL when defaulted). */
static void np_spec_parts(Expr* spec, Expr** var, Expr** imin, Expr** di) {
    size_t n = spec->data.function.arg_count;
    *var = spec->data.function.args[0];
    *imin = n >= 3 ? spec->data.function.args[1] : NULL;
    *di = n >= 4 ? spec->data.function.args[3] : NULL;
}

/* imin + j di, unevaluated (fresh). */
static Expr* np_index(Expr* imin, Expr* di, long j) {
    return np_fn2("Plus", imin ? expr_copy(imin) : expr_new_integer(1),
                  np_fn2("Times", expr_new_integer(j), di ? expr_copy(di) : expr_new_integer(1)));
}

/* Table[body, {var, imin + a di, imin + b di, di}] (adopts body). */
static Expr* np_table(Expr* body, Expr* spec, long a, long b) {
    Expr *var, *imin, *di;
    np_spec_parts(spec, &var, &imin, &di);
    Expr* it[4] = { expr_copy(var), np_index(imin, di, a), np_index(imin, di, b),
                    di ? expr_copy(di) : expr_new_integer(1) };
    return np_fn2("Table", body, expr_new_function(expr_new_symbol(SYM_List), it, 4));
}

/* Block[{var}, Hold[Evaluate[body]]], unwrapped (adopts body): evaluate a
 * factor-derived expression symbolically without a global value of the index
 * leaking in -- neither inside the Block nor when its result is re-evaluated
 * outside it, which the Hold prevents.  NULL if the shape is lost. */
static Expr* np_block_eval(Expr* body, Expr* spec) {
    Expr* v[1] = { expr_copy(spec->data.function.args[0]) };
    Expr* held = np_apply1("Hold", np_apply1("Evaluate", body));
    Expr* r = eval_and_free(np_fn2("Block", expr_new_function(expr_new_symbol(SYM_List), v, 1),
                                   held));
    Expr* out = NULL;
    if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
        && r->data.function.head->data.symbol.name == SYM_Hold
        && r->data.function.arg_count == 1) {
        out = r->data.function.args[0];
        r->data.function.args[0] = NULL;
    }
    expr_free(r);
    return out;
}

/* Probe the leading factors f(imin + j di) through g = f - 1.  True when they
 * settle (NP_PROBE_SETTLED consecutive |g| < 1/2) within NP_PROBE_MAX factors,
 * each numeric; *peel is then the number of factors before the settled run
 * when a factor -1 was among them, else 0.  Messages from the probe are the
 * evaluation's own business (NSum meets the same terms), so it runs Quiet. */
static bool np_probe_head(Expr* g, Expr* spec, long* peel) {
    int settled = 0;
    bool minus_one = false;
    *peel = 0;
    for (long a = 0; a < NP_PROBE_MAX; a += NP_PROBE_CHUNK) {
        Expr* table = np_table(np_apply1("N", expr_copy(g)), spec, a, a + NP_PROBE_CHUNK - 1);
        /* Table packs a numeric result: materialise it to walk the elements. */
        Expr* qt = np_apply1("Quiet", table);
        Expr* vals = pack_eval_plain(qt);
        expr_free(qt);
        bool ok = vals && vals->type == EXPR_FUNCTION
               && vals->data.function.head->type == EXPR_SYMBOL
               && vals->data.function.head->data.symbol.name == SYM_List;
        for (size_t i = 0; ok && i < vals->data.function.arg_count; i++) {
            double re, im;
            if (!np_to_complex(vals->data.function.args[i], &re, &im)
                || !isfinite(re) || !isfinite(im)) { ok = false; break; }
            if (re == -2.0 && im == 0.0) minus_one = true;           /* f = -1 */
            settled = (hypot(re, im) < 0.5) ? settled + 1 : 0;
            if (settled >= NP_PROBE_SETTLED) {
                if (minus_one) *peel = a + (long)i + 1 - NP_PROBE_SETTLED;
                expr_free(vals);
                return true;
            }
        }
        expr_free(vals);
        if (!ok) return false;
    }
    return false;
}

/* The log-sum of an infinite product: Prod f = H * Exp[NSum[summand, spec']],
 * spec' starting `peel` factors later and H their direct product. */
typedef struct {
    Expr* summand;       /* owned                                         */
    Expr* twin;          /* owned derivative twin (Log[f]) or NULL        */
    long  peel;          /* leading factors multiplied out directly       */
} NpLogSum;

/* Build the summand for an infinite product of f over `spec` (see above).
 * Borrows f. */
static void np_log_summand(Expr* f, Expr* spec, NpLogSum* out) {
    const char* var = spec->data.function.args[0]->data.symbol.name;
    out->twin = NULL; out->peel = 0;
    Expr* fe = np_block_eval(expr_copy(f), spec);
    if (fe && fe->type == EXPR_FUNCTION && fe->data.function.head->type == EXPR_SYMBOL
        && fe->data.function.head->data.symbol.name == SYM_Power
        && fe->data.function.arg_count == 2
        && fe->data.function.args[0]->type == EXPR_SYMBOL
        && fe->data.function.args[0]->data.symbol.name == SYM_E) {
        out->summand = expr_copy(fe->data.function.args[1]);        /* E^h -> h */
        expr_free(fe);
        return;
    }
    expr_free(fe);
    Expr* g = np_block_eval(np_fn2("Plus", expr_copy(f), expr_new_integer(-1)), spec);
    if (g && np_has_constant_addend(g, var)) {
        Expr* t = np_block_eval(np_apply1("Together", expr_copy(g)), spec);
        if (t && !np_has_constant_addend(t, var)) { expr_free(g); g = t; }
        else expr_free(t);
    }
    if (!g || !np_mentions(g, var) || !np_probe_head(g, spec, &out->peel)) {
        expr_free(g);
        out->peel = 0;
        out->summand = np_apply1("Log", expr_copy(f));
        return;
    }
    /* 2 ArcTanh[g (2 + g)^-1] */
    Expr* z = np_fn2("Times", expr_copy(g),
                     np_fn2("Power", np_fn2("Plus", expr_new_integer(2), g),
                            expr_new_integer(-1)));
    out->summand = np_fn2("Times", expr_new_integer(2), np_apply1("ArcTanh", z));
    out->twin = np_apply1("Log", expr_copy(f));
}

/* Times @@ Table[N[f (, digits)], ...] over the first `peel` factors. */
static Expr* np_head_product(Expr* f, Expr* spec, long peel, double digits) {
    Expr* nf = digits > 0 ? np_fn2("N", expr_copy(f), expr_new_integer((long)(digits + 0.5)))
                          : np_apply1("N", expr_copy(f));
    return eval_and_free(np_fn2("Apply", expr_new_symbol(SYM_Times),
                                np_table(nf, spec, 0, peel - 1)));
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

    /* Single-element iterators {c}, by the Product convention: a real c is
     * c terms, {dummy, 1, c}; anything else is NProduct::itraw (see
     * iter_numeric_count_specs). */
    {
        bool declined = false;
        Expr* rw = iter_numeric_count_specs(res, pos_end, "NProduct", "NProduct`Private`count$", &declined);
        if (declined) return NULL;
        if (rw) {
            Expr* r = builtin_nprod(rw);
            expr_free(rw);
            return r;
        }
    }

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

    /* Build the inner NSum[Log[eff_body], spec0, <mapped options>].  An
     * infinite single product sums the cancellation-free form of Log (see
     * np_log_summand); a finite or nested one keeps Log[f]. */
    bool infinite = false;
    if (nspecs == 1) {
        size_t hi_at = spec0->data.function.arg_count == 2 ? 1 : 2;
        Expr* hi = eval_and_free(expr_copy(spec0->data.function.args[hi_at]));
        infinite = np_is_infinity(hi);
        expr_free(hi);
    }
    NpLogSum ls = { NULL, NULL, 0 };
    Expr* tail_spec;
    if (infinite) {
        np_log_summand(eff_body, spec0, &ls);
        expr_free(eff_body);
    } else {
        ls.summand = np_apply1("Log", eff_body);
    }
    if (ls.peel > 0) {
        /* {var, imin + peel di, imax (, di)} */
        Expr *pv, *pmin, *pdi;
        np_spec_parts(spec0, &pv, &pmin, &pdi);
        size_t n = spec0->data.function.arg_count;
        Expr* it[4];
        size_t w = 0;
        it[w++] = expr_copy(pv);
        it[w++] = np_index(pmin, pdi, ls.peel);
        it[w++] = expr_copy(spec0->data.function.args[n == 2 ? 1 : 2]);
        if (n >= 4) it[w++] = expr_copy(pdi);
        tail_spec = expr_new_function(expr_new_symbol(SYM_List), it, w);
    } else {
        tail_spec = expr_copy(spec0);
    }
    Expr* logbody = ls.summand;

    Expr** nv = malloc(sizeof(Expr*) * 12);
    if (!nv) { expr_free(logbody); expr_free(ls.twin); expr_free(tail_spec); return NULL; }
    size_t nw = 0;
    nv[nw++] = logbody;
    nv[nw++] = tail_spec;
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

    /* the numeric Log-sum, or special; its messages say NProduct */
    Expr* s = nsum_eval_for_nproduct(nsum, ls.twin);
    expr_free(ls.twin);

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
    if (ls.peel > 0)
        result = eval_and_free(np_fn2("Times",
                     np_head_product(body, spec0, ls.peel,
                                     o.mpfr ? o.wdigits + NP_GUARD_DIGITS : 0.0),
                     result));

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
