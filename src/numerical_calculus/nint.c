/*
 * nint.c — NIntegrate[f, {x, xmin, xmax}, opts]   (see nint.h)
 *
 * Phase 1: one-dimensional integrals over a finite real interval at machine
 * precision, via globally-adaptive Gauss-Kronrod (gkadapt).  HoldAll: the
 * integrand and bounds are held, the bounds are evaluated to numbers, then the
 * integration variable is Block-localised and the integrand is evaluated /
 * numericalised at each sample point.  Subsequent phases layer endpoint
 * singularities, infinite ranges, complex contours, arbitrary precision,
 * multidimensional iteration, oscillatory and Monte-Carlo methods, Exclusions
 * and principal values on top of this same sampling machinery.
 *
 * Memory contract: never frees `res`; returns a fresh Expr* or NULL; restores
 * the variable binding on every return path.
 */

#include "nint.h"
#include "gkadapt.h"
#include "denint.h"
#include "oscint.h"
#include "mcint.h"
#include "cubature.h"

#include <complex.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ *
 *  Diagnostics                                                        *
 * ------------------------------------------------------------------ */

static void ni_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NIntegrate::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers                                   *
 * ------------------------------------------------------------------ */

static bool ni_to_double_real(Expr* e, double* out) {
    if (!e) return false;
    int64_t rn, rd;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational(e, &rn, &rd)) { *out = (double)rn / (double)rd; return true; }
            return false;
        default: return false;
    }
}

static bool ni_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (ni_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (ni_to_double_real(re, &r) && ni_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* ni_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    double ar = fabs(r), ai = fabs(i);
    /* Chop a component that is pure numerical noise relative to the other — in
     * particular the ~machine-epsilon imaginary part that Power[-x, n] leaves
     * behind when a real integrand is sampled at negative abscissae.  Genuine
     * complex (contour) results have components of comparable magnitude and are
     * untouched. */
    if (ai <= 1e-12 * ar) i = 0.0;
    else if (ar <= 1e-12 * ai) r = 0.0;
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors nsum.c)                      *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} NiBind;

static void ni_bind_snapshot(NiBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void ni_bind_clear_temp(SymbolDef* def) {
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
}

static void ni_bind_set(NiBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    ni_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);   /* copies value internally */
    expr_free(sym);
}

static void ni_bind_restore(NiBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    ni_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Integrand evaluation context                                       *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       body;    /* borrowed integrand (free in the variable)      */
    NiBind*     bind;
    NumericSpec spec;
    /* Abscissa map: the variable is bound to (x_scale·x + x_shift), used to
     * reflect a (-∞, b] range onto a half line and (later) to parametrise a
     * complex contour. */
    double _Complex x_scale;
    double _Complex x_shift;
} NiCtx;

static void ni_ctx_init(NiCtx* c, Expr* body, NiBind* bind, NumericSpec spec) {
    c->body = body; c->bind = bind; c->spec = spec;
    c->x_scale = 1.0; c->x_shift = 0.0;
}

/* Build the numeric leaf the variable is bound to from a mapped abscissa. */
static Expr* ni_abscissa_expr(const NiCtx* c, double _Complex z) {
    double _Complex w = c->x_scale * z + c->x_shift;
    if (cimag(w) == 0.0) return expr_new_real(creal(w));
    return make_complex(expr_new_real(creal(w)), expr_new_real(cimag(w)));
}

/* Evaluate the integrand with the variable bound to `value` (consumed).
 * Returns the numericalised result (caller frees) or NULL. */
static Expr* ni_eval_at(NiCtx* c, Expr* value) {
    ni_bind_set(c->bind, value);
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->body));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

/* Machine sample callback: f at a (possibly mapped) real abscissa x. */
static bool ni_sample_machine(void* vctx, double x, double _Complex* out) {
    NiCtx* c = (NiCtx*)vctx;
    Expr* num = ni_eval_at(c, ni_abscissa_expr(c, x));
    if (!num) return false;
    bool ok = ni_to_complex(num, out);
    if (ok && (!isfinite(creal(*out)) || !isfinite(cimag(*out)))) ok = false;
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
/* MPFR sample callback: f at the real abscissa x (mapped by the real part of
 * the affine map — MPFR is used for real integration ranges).  (out_re,out_im)
 * are pre-initialised at the working precision by the DE kernel. */
static bool ni_sample_mpfr(void* vctx, const mpfr_t x, mpfr_t out_re, mpfr_t out_im) {
    NiCtx* c = (NiCtx*)vctx;
    mpfr_t xm;
    mpfr_init2(xm, mpfr_get_prec(x));
    mpfr_mul_d(xm, x, creal(c->x_scale), MPFR_RNDN);
    mpfr_add_d(xm, xm, creal(c->x_shift), MPFR_RNDN);
    Expr* xe = expr_new_mpfr_copy(xm);
    mpfr_clear(xm);
    Expr* num = ni_eval_at(c, xe);
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}

/* Build a result from a guarded-precision (re,im) pair: round to target_bits,
 * chop a noise-level imaginary part, return real or Complex. */
static Expr* ni_mpfr_build(const mpfr_t re, const mpfr_t im, long target_bits) {
    mpfr_t r, i;
    mpfr_init2(r, (mpfr_prec_t)target_bits);
    mpfr_init2(i, (mpfr_prec_t)target_bits);
    mpfr_set(r, re, MPFR_RNDN);
    mpfr_set(i, im, MPFR_RNDN);
    double ar = fabs(mpfr_get_d(r, MPFR_RNDN));
    double ai = fabs(mpfr_get_d(i, MPFR_RNDN));
    Expr* out;
    if (mpfr_zero_p(i) || ai <= 1e-12 * ar) {
        out = expr_new_mpfr_copy(r);
    } else if (ar <= 1e-12 * ai) {
        mpfr_set_zero(r, 1);
        out = make_complex(expr_new_mpfr_copy(r), expr_new_mpfr_copy(i));
    } else {
        out = make_complex(expr_new_mpfr_copy(r), expr_new_mpfr_copy(i));
    }
    mpfr_clear(r);
    mpfr_clear(i);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Iterator-bound classification                                      *
 * ------------------------------------------------------------------ */

/* If `e` is an infinite quantity, set *dir to its (complex) unit direction and
 * return true: Infinity -> 1, -Infinity (Times[-1,Infinity]) -> -1,
 * I Infinity -> I, DirectedInfinity[d] -> d, ComplexInfinity -> 0. */
static bool ni_infinity_dir(Expr* e, double _Complex* dir) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == SYM_Infinity)        { *dir = 1.0; return true; }
        if (e->data.symbol == SYM_ComplexInfinity) { *dir = 0.0; return true; }
        return false;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_DirectedInfinity) {
            if (e->data.function.arg_count == 0) { *dir = 0.0; return true; }
            return ni_to_complex(e->data.function.args[0], dir);
        }
        if (h == SYM_Times) {
            double _Complex prod = 1.0; bool found = false;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                double _Complex d;
                Expr* a = e->data.function.args[i];
                if (ni_infinity_dir(a, &d)) { found = true; prod *= d; }
                else if (ni_to_complex(a, &d)) prod *= d;
                else return false;
            }
            if (found) { *dir = prod; return true; }
        }
    }
    return false;
}

/* Classification of an integration bound. */
typedef enum {
    NI_B_FIN,    /* finite real, value in z (imag 0)        */
    NI_B_CPLX,   /* finite complex, value in z              */
    NI_B_PINF,   /* +Infinity                               */
    NI_B_NINF,   /* -Infinity                               */
    NI_B_CRAY,   /* complex infinite ray, unit dir in dir   */
    NI_B_BAD     /* ComplexInfinity / symbolic / non-numeric */
} NiBoundKind;

typedef struct { NiBoundKind kind; double _Complex z; double _Complex dir; } NiBound;

static void ni_classify(Expr* e, NiBound* b) {
    b->z = 0.0; b->dir = 0.0;
    double _Complex dir;
    if (ni_infinity_dir(e, &dir)) {
        if (cimag(dir) != 0.0) { b->kind = NI_B_CRAY; b->dir = dir; return; }
        if (creal(dir) > 0.0) { b->kind = NI_B_PINF; return; }
        if (creal(dir) < 0.0) { b->kind = NI_B_NINF; return; }
        b->kind = NI_B_BAD;                          /* ComplexInfinity */
        return;
    }
    double _Complex z;
    if (ni_to_complex(e, &z)) {
        b->z = z;
        b->kind = (cimag(z) == 0.0) ? NI_B_FIN : NI_B_CPLX;
        return;
    }
    b->kind = NI_B_BAD;
}

/* Evaluate and numericalise a bound/contour-node spec (consuming a copy of `e`)
 * so symbolic constants like Pi, E, Sqrt[2] become finite reals, while leaving
 * Infinity / DirectedInfinity intact for ni_classify to recognise. */
static Expr* ni_num_endpoint(Expr* e) {
    Expr* ev = eval_and_free(expr_copy(e));
    if (!ev) return NULL;
    Expr* num = numericalize(ev, numeric_machine_spec());
    if (!num) return ev;            /* keep evaluated form if numericalize bails */
    expr_free(ev);
    return num;
}

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef enum {
    NI_AUTO = 0, NI_GK, NI_DE, NI_TRAP, NI_LEVIN,
    NI_MC, NI_QMC, NI_AMC, NI_PV,
    NI_UNIMPL          /* recognised name with no implementation yet            */
} NiMethod;

/* Methods implemented in the current build.  An explicitly requested method
 * outside this set causes NIntegrate to warn and stay unevaluated, so missing
 * implementations are visible rather than silently approximated.  Each phase
 * adds its method here as it lands. */
static bool ni_method_implemented(NiMethod m) {
    return m == NI_AUTO || m == NI_GK || m == NI_DE || m == NI_LEVIN
        || m == NI_MC || m == NI_QMC || m == NI_AMC || m == NI_PV;
}

typedef struct {
    NiMethod method;
    const char* method_name; /* requested Method name (for the warning); NULL=Automatic */
    bool   prec_mpfr;       /* WorkingPrecision selects MPFR                 */
    long   bits;            /* MPFR internal working precision (guarded)      */
    long   target_bits;     /* user-requested precision; result rounded back  */
    double acc_goal;        /* AccuracyGoal digits (-1 => Infinity)          */
    double prec_goal;       /* PrecisionGoal digits (-1 => Automatic)        */
    int    max_recursion;   /* -1 => Automatic                                */
    int    min_recursion;   /* default 0                                      */
    long   max_points;      /* -1 => Automatic                                */
    Expr*  exclusions;      /* borrowed; NULL if none                         */
} NiOpts;

static bool ni_is_known_option(const char* s) {
    return s == SYM_Method || s == SYM_WorkingPrecision
        || s == SYM_AccuracyGoal || s == SYM_PrecisionGoal
        || s == SYM_MaxRecursion || s == SYM_MinRecursion
        || s == SYM_MaxPoints || s == SYM_Exclusions
        || s == SYM_EvaluationMonitor;
}

static bool ni_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && ni_is_known_option(lhs->data.symbol);
}

/* Resolve a Method name to its engine.  Strategies/rules that this build
 * actually realises map to a live engine; recognised-but-unrealised rules map
 * to their semantic enum (gated by ni_method_implemented) so flipping them on
 * later is a one-line change; anything unknown is NI_UNIMPL. */
static NiMethod ni_method_from_string(const char* s) {
    /* GlobalAdaptive is exactly the adaptive Gauss-Kronrod strategy we run. */
    if (!strcmp(s, "GlobalAdaptive") || !strcmp(s, "GaussKronrodRule")) return NI_GK;
    if (!strcmp(s, "DoubleExponential")) return NI_DE;
    if (!strcmp(s, "Trapezoidal") || !strcmp(s, "TrapezoidalRule")) return NI_TRAP;
    if (!strcmp(s, "LevinRule")) return NI_LEVIN;
    if (!strcmp(s, "MonteCarlo")) return NI_MC;
    if (!strcmp(s, "QuasiMonteCarlo")) return NI_QMC;
    if (!strcmp(s, "AdaptiveMonteCarlo") || !strcmp(s, "AdaptiveQuasiMonteCarlo")) return NI_AMC;
    if (!strcmp(s, "PrincipalValue")) return NI_PV;
    /* LocalAdaptive, ClenshawCurtisRule, NewtonCotesRule, LobattoKronrodRule,
     * MultidimensionalRule, CartesianRule, MultipanelRule, RiemannRule, and any
     * unknown name: not implemented as a distinct rule. */
    return NI_UNIMPL;
}

static bool ni_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!ni_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool ni_parse_goal(Expr* v, double* out) {
    if (v->type == EXPR_SYMBOL
        && (v->data.symbol == SYM_Infinity || v->data.symbol == SYM_Automatic)) {
        *out = -1.0; return true;
    }
    double d;
    if (ni_to_double_real(v, &d) && d > 0.0) { *out = d; return true; }
    return false;
}

static bool ni_parse_int(Expr* v, int* out, bool allow_auto) {
    if (allow_auto && v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic) {
        *out = -1; return true;
    }
    if (v->type == EXPR_INTEGER && v->data.integer >= 0) { *out = (int)v->data.integer; return true; }
    return false;
}

static bool ni_apply_option(Expr* rule, NiOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Method) {
        Expr* m = rhs;
        /* Method -> {"name", subopts...} or Method -> "name" or symbol. */
        if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL
            && m->data.function.head->data.symbol == SYM_List
            && m->data.function.arg_count >= 1)
            m = m->data.function.args[0];
        if (m->type == EXPR_STRING) {
            o->method_name = m->data.string;
            o->method = ni_method_from_string(m->data.string);
            return true;
        }
        if (m->type == EXPR_SYMBOL) {
            if (m->data.symbol == SYM_Automatic) { o->method = NI_AUTO; o->method_name = NULL; return true; }
            o->method_name = m->data.symbol;
            o->method = ni_method_from_string(m->data.symbol);
            return true;
        }
        o->method = NI_AUTO;
        return true;
    }
    if (name == SYM_WorkingPrecision) {
        if (!ni_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            ni_warn("badopt", "invalid WorkingPrecision value"); return false;
        }
        return true;
    }
    if (name == SYM_AccuracyGoal)  return ni_parse_goal(rhs, &o->acc_goal);
    if (name == SYM_PrecisionGoal) return ni_parse_goal(rhs, &o->prec_goal);
    if (name == SYM_MaxRecursion)  return ni_parse_int(rhs, &o->max_recursion, true);
    if (name == SYM_MinRecursion)  return ni_parse_int(rhs, &o->min_recursion, true);
    if (name == SYM_MaxPoints) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) { o->max_points = -1; return true; }
        double d;
        if (ni_to_double_real(rhs, &d) && d >= 1.0) { o->max_points = (long)d; return true; }
        ni_warn("badopt", "MaxPoints must be a positive integer or Automatic"); return false;
    }
    if (name == SYM_Exclusions)       { o->exclusions = rhs; return true; }
    if (name == SYM_EvaluationMonitor) return true;   /* accepted, ignored */
    return false;
}

/* Iterator spec: a List whose first element is a symbol and with >= 3 parts
 * ({x, xmin, xmax} or a {x, x0, x1, ...} contour / split list). */
static bool ni_is_spec(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    return n >= 3 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* ------------------------------------------------------------------ *
 *  Convergence-tolerance helper                                       *
 * ------------------------------------------------------------------ */

/* Effective relative / absolute tolerances (machine).  PrecisionGoal default
 * for machine working precision is half the machine precision; AccuracyGoal
 * default is Infinity (no absolute target). */
static void ni_machine_tols(const NiOpts* o, double* reltol, double* abstol) {
    double pg = o->prec_goal;
    if (pg < 0.0) pg = 0.5 * NUMERIC_MACHINE_PRECISION_DIGITS;   /* ~8 digits */
    *reltol = pow(10.0, -pg);
    if (*reltol < 1e-14) *reltol = 1e-14;
    *abstol = (o->acc_goal > 0.0) ? pow(10.0, -o->acc_goal) : 0.0;
}

/* Map MaxRecursion to a tanh-sinh / sinh-sinh level count. */
static int ni_de_levels(const NiOpts* o) {
    int r = o->max_recursion;
    if (r < 0) return 13;
    if (r < 6) return 6;
    if (r > 24) return 24;
    return r;
}

#ifdef USE_MPFR
/* Relative tolerance for the MPFR double-exponential path.  PrecisionGoal
 * defaults to ~10 digits below the working precision (Mathematica's rule). */
static double ni_mpfr_reltol(const NiOpts* o) {
    long bits = o->target_bits > 0 ? o->target_bits : o->bits;
    double target_digits = numeric_bits_to_digits(bits);
    double pg = (o->prec_goal > 0.0) ? o->prec_goal : (target_digits - 10.0);
    if (pg < 1.0) pg = 1.0;
    double rt = pow(10.0, -pg);
    if (!(rt > 0.0)) rt = 1e-300;   /* double underflow guard at very high precision */
    return rt;
}

static int ni_mpfr_levels(const NiOpts* o) {
    int L = ni_de_levels(o);
    return L < 16 ? 16 : L;
}
#endif

/* Map MaxRecursion / MaxPoints to gkadapt resource caps. */
static void ni_gk_budget(const NiOpts* o, int* max_subdiv, long* max_eval) {
    int r = o->max_recursion;
    if (r < 0) {
        *max_subdiv = 800;                 /* Automatic */
    } else if (r == 0) {
        *max_subdiv = 0;                   /* no subdivision */
    } else {
        int cap = r < 22 ? (1 << r) : (1 << 22);
        if (cap > 4000000) cap = 4000000;
        *max_subdiv = cap;
    }
    if (o->min_recursion > 0 && *max_subdiv < o->min_recursion)
        *max_subdiv = o->min_recursion;
    *max_eval = (o->max_points > 0) ? o->max_points : 0;
}

/* ------------------------------------------------------------------ *
 *  One-dimensional driver (machine, finite real interval)             *
 * ------------------------------------------------------------------ */

/* Best-of-attempts accumulator: prefer a converged estimate, else the one with
 * the smaller error estimate. */
typedef struct { bool have, conv; double _Complex val; double err; } NiAtt;

static void ni_consider(NiAtt* best, NiAtt a) {
    if (!a.have) return;
    if (!best->have
        || (a.conv && !best->conv)
        || (a.conv == best->conv && a.err < best->err))
        *best = a;
}

static NiAtt ni_try_gk(NiCtx* ctx, double a, double b, const NiOpts* o,
                       double reltol, double abstol) {
    int max_subdiv; long max_eval;
    ni_gk_budget(o, &max_subdiv, &max_eval);
    GkResult R = gk_integrate_machine(ni_sample_machine, ctx, a, b,
                                      abstol, reltol, max_subdiv, max_eval, true);
    NiAtt out;
    out.have = (R.status != GK_NONNUMERIC);
    out.conv = (R.status == GK_OK);
    out.val = R.value;
    out.err = R.abs_err;
    return out;
}

static NiAtt ni_try_de(NiCtx* ctx, double a, double b, const NiOpts* o, double reltol) {
    double _Complex v; double err;
    bool conv = denint_tanhsinh_machine(ni_sample_machine, ctx, a, b, reltol,
                                        ni_de_levels(o), &v, &err);
    NiAtt out;
    out.have = isfinite(creal(v)) && isfinite(cimag(v));
    out.conv = conv;
    out.val = v;
    out.err = err;
    return out;
}

#define NI_OSC_MAX_PANELS 200000

/* Oscillatory panel quadrature over a finite interval (half-period panels). */
static NiAtt ni_try_osc(NiCtx* ctx, double a, double b, double reltol) {
    double _Complex v; double err;
    bool conv = osc_integrate_machine(ni_sample_machine, ctx, a, b, false,
                                      reltol, NI_OSC_MAX_PANELS, &v, &err);
    NiAtt out;
    out.have = conv && isfinite(creal(v)) && isfinite(cimag(v));
    out.conv = conv;
    out.val = v;
    out.err = err;
    return out;
}

/* Integrate over the finite real parameter interval [a,b] with the variable
 * already configured in `ctx` (directly, or through an affine abscissa map for
 * a complex segment).  Gauss-Kronrod samples strictly interior abscissae, so it
 * never evaluates a singular endpoint (no spurious 1/0 messages); for Automatic
 * it falls back to the endpoint-robust tanh-sinh rule when it cannot reach the
 * tolerance (the signature of an endpoint singularity or slow convergence). */
static NiAtt ni_core_finite(NiCtx* ctx, double a, double b, const NiOpts* o) {
    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    NiAtt best = { false, false, 0.0, INFINITY };
    if (o->method == NI_LEVIN) {
        ni_consider(&best, ni_try_osc(ctx, a, b, reltol));
    } else if (o->method == NI_DE) {
        ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
    } else {
        ni_consider(&best, ni_try_gk(ctx, a, b, o, reltol, abstol));
        if (!best.conv && (o->method == NI_AUTO || !best.have))
            ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
        /* Oscillatory fallback: neither GK (a subdivision per half-period) nor
         * tanh-sinh (assumes endpoint decay) converges on a many-period
         * integrand; integrate between the zeros instead. */
        if (!best.conv && o->method == NI_AUTO)
            ni_consider(&best, ni_try_osc(ctx, a, b, reltol));
    }
    return best;
}

static Expr* ni_run_1d_finite_real(Expr* body, const char* var,
                                   double a, double b, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());

    double sign = 1.0;
    if (a > b) { double t = a; a = b; b = t; sign = -1.0; }

#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        long bits = o->bits;
        mpfr_t am, bm, re, im;
        mpfr_init2(am, bits); mpfr_init2(bm, bits);
        mpfr_init2(re, bits); mpfr_init2(im, bits);
        mpfr_set_d(am, a, MPFR_RNDN); mpfr_set_d(bm, b, MPFR_RNDN);
        double abserr;
        bool conv = denint_tanhsinh_mpfr(ni_sample_mpfr, &ctx, am, bm, bits,
                                         ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                         re, im, &abserr);
        if (sign < 0) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
        ni_bind_restore(&bind);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(am, bm, re, im, (mpfr_ptr)0);
        if (out && !conv)
            ni_warn("ncvb", "high-precision integral did not converge "
                    "(error estimate %.3g)", abserr);
        return out;
    }
#endif

    NiAtt best = ni_core_finite(&ctx, a, b, o);
    ni_bind_restore(&bind);

    if (!best.have) return NULL;        /* no numeric estimate at all */
    if (!best.conv)
        ni_warn("ncvb", "failed to converge to prescribed accuracy "
                "(error estimate %.3g)", best.err);
    return ni_from_complex_d(best.val * sign);
}

/* Integral along a piecewise-linear contour through `nodes[0..nnodes-1]` in the
 * complex plane.  Each segment z(t) = z0 + t(z1-z0), t in [0,1], contributes
 * (z1-z0) * ∫_0^1 f(z(t)) dt via the affine abscissa map in NiCtx.  Real nodes
 * are the degenerate case: the contour collapses to the real axis and the
 * intermediate nodes simply split the integral (handling interior singularities
 * at those points through the tanh-sinh rule). */
static Expr* ni_run_contour(Expr* body, const char* var,
                            const double _Complex* nodes, int nnodes,
                            const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());

    double _Complex total = 0.0;
    bool all_conv = true, any_fail = false;
    for (int i = 0; i < nnodes - 1; i++) {
        double _Complex z0 = nodes[i], z1 = nodes[i + 1];
        if (z1 == z0) continue;                 /* zero-length segment */
        ctx.x_scale = z1 - z0;
        ctx.x_shift = z0;
        NiAtt a = ni_core_finite(&ctx, 0.0, 1.0, o);
        if (!a.have) { any_fail = true; break; }
        total += (z1 - z0) * a.val;
        if (!a.conv) all_conv = false;
    }
    ni_bind_restore(&bind);

    if (any_fail) return NULL;
    if (!all_conv)
        ni_warn("ncvb", "contour integral did not converge to prescribed accuracy");
    return ni_from_complex_d(total);
}

/* Integral along the complex ray  z(s) = z0 + dir·s,  s in [0, ∞), via exp-sinh
 * (dz = dir ds).  The overall result is multiplied by `sign` (for a reversed
 * orientation).  The real semi-infinite cases keep their dedicated path. */
static Expr* ni_run_ray(Expr* body, const char* var, double _Complex z0,
                        double _Complex dir, double sign, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    ctx.x_scale = dir;
    ctx.x_shift = z0;

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    double _Complex val; double abserr;
    bool conv = dequad_halfline_machine(ni_sample_machine, &ctx, 0.0, reltol,
                                        ni_de_levels(o), &val, &abserr);
    ni_bind_restore(&bind);

    val *= dir * sign;
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv)
        ni_warn("ncvb", "complex-ray integral did not converge "
                "(error estimate %.3g)", abserr);
    return ni_from_complex_d(val);
}

/* Semi-infinite range via exp-sinh.  `a` is the half-line lower limit in the
 * sampling variable; reflect=true binds the integrand variable to -u (mapping a
 * (-∞, b] range onto [-b, ∞)); the overall result is multiplied by `sign`. */
static Expr* ni_run_1d_halfline(Expr* body, const char* var, double a,
                                bool reflect, double sign, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    if (reflect) ctx.x_scale = -1.0;

#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        long bits = o->bits;
        mpfr_t am, re, im;
        mpfr_init2(am, bits); mpfr_init2(re, bits); mpfr_init2(im, bits);
        mpfr_set_d(am, a, MPFR_RNDN);
        double abserr;
        bool conv = dequad_halfline_mpfr(ni_sample_mpfr, &ctx, am, bits,
                                         ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                         re, im, &abserr);
        if (sign < 0) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
        ni_bind_restore(&bind);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(am, re, im, (mpfr_ptr)0);
        if (out && !conv)
            ni_warn("ncvb", "high-precision semi-infinite integral did not "
                    "converge (error estimate %.3g)", abserr);
        return out;
    }
#endif

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    double _Complex val; double abserr;
    bool conv;
    if (o->method == NI_LEVIN) {
        conv = osc_integrate_machine(ni_sample_machine, &ctx, a, 0.0, true,
                                     reltol, 600, &val, &abserr);
    } else {
        conv = dequad_halfline_machine(ni_sample_machine, &ctx, a, reltol,
                                       ni_de_levels(o), &val, &abserr);
        /* Oscillatory fallback: a slowly-decaying oscillatory tail (Sin[x]/x,
         * Bessel functions, …) defeats exp-sinh; integrate between the zeros
         * and accelerate the partial sums with Wynn's epsilon. */
        if (!conv && o->method == NI_AUTO) {
            double _Complex ov; double oerr;
            if (osc_integrate_machine(ni_sample_machine, &ctx, a, 0.0, true,
                                      reltol, 600, &ov, &oerr)) {
                val = ov; abserr = oerr; conv = true;
            }
        }
    }
    ni_bind_restore(&bind);

    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv)
        ni_warn("ncvb", "semi-infinite integral did not converge "
                "(error estimate %.3g)", abserr);
    return ni_from_complex_d(val * sign);
}

/* Doubly-infinite range via sinh-sinh. */
static Expr* ni_run_1d_wholeline(Expr* body, const char* var, double sign,
                                 const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());

#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        long bits = o->bits;
        mpfr_t re, im;
        mpfr_init2(re, bits); mpfr_init2(im, bits);
        double abserr;
        bool conv = denint_sinhsinh_mpfr(ni_sample_mpfr, &ctx, bits,
                                         ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                         re, im, &abserr);
        if (sign < 0) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
        ni_bind_restore(&bind);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(re, im, (mpfr_ptr)0);
        if (out && !conv)
            ni_warn("ncvb", "high-precision doubly-infinite integral did not "
                    "converge (error estimate %.3g)", abserr);
        return out;
    }
#endif

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    double _Complex val; double abserr;
    bool conv = denint_sinhsinh_machine(ni_sample_machine, &ctx, reltol,
                                        ni_de_levels(o), &val, &abserr);
    ni_bind_restore(&bind);

    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv)
        ni_warn("ncvb", "doubly-infinite integral did not converge "
                "(error estimate %.3g)", abserr);
    return ni_from_complex_d(val * sign);
}

/* ------------------------------------------------------------------ *
 *  Monte-Carlo cubature (high dimension / region integrands)          *
 * ------------------------------------------------------------------ */

/* Does the (held) integrand contain a region head (Boole / UnitStep)? */
static bool ni_body_has_region(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Boole || h == SYM_UnitStep) return true;
    }
    if (ni_body_has_region(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ni_body_has_region(e->data.function.args[i])) return true;
    return false;
}

typedef struct {
    Expr* body;
    NiBind* binds;
    size_t d;
    NumericSpec spec;
} NiMcCtx;

static bool ni_mc_sample(void* vctx, const double* x, size_t d, double _Complex* out) {
    NiMcCtx* c = (NiMcCtx*)vctx;
    for (size_t j = 0; j < d; j++) {
        Expr* xe = expr_new_real(x[j]);
        ni_bind_set(&c->binds[j], xe);
        expr_free(xe);
    }
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->body));
    if (!raw) return false;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    if (!num) return false;
    bool ok = ni_to_complex(num, out);
    expr_free(num);
    return ok;
}

/* Monte-Carlo over the box defined by `specs` (each {v, vmin, vmax}, constant
 * finite bounds).  Returns NULL (with a warning) if the box is not well-formed. */
static Expr* ni_run_mc(Expr* body, Expr** specs, size_t d, const NiOpts* o, bool quasi) {
    if (d == 0 || d > 24) return NULL;
    double* a = malloc(d * sizeof(double));
    double* b = malloc(d * sizeof(double));
    NiBind* binds = malloc(d * sizeof(NiBind));
    bool ok = true;
    for (size_t j = 0; j < d && ok; j++) {
        Expr* s = specs[j];
        if (s->data.function.arg_count != 3
            || s->data.function.args[0]->type != EXPR_SYMBOL) { ok = false; break; }
        Expr* lo = eval_and_free(expr_copy(s->data.function.args[1]));
        Expr* hi = eval_and_free(expr_copy(s->data.function.args[2]));
        if (!ni_to_double_real(lo, &a[j]) || !ni_to_double_real(hi, &b[j])) ok = false;
        expr_free(lo); expr_free(hi);
    }
    if (!ok) {
        ni_warn("mc", "Monte-Carlo requires constant finite real bounds");
        free(a); free(b); free(binds);
        return NULL;
    }

    for (size_t j = 0; j < d; j++)
        ni_bind_snapshot(&binds[j], specs[j]->data.function.args[0]->data.symbol);

    NiMcCtx ctx;
    ctx.body = body; ctx.binds = binds; ctx.d = d;
    ctx.spec = numeric_machine_spec();

    /* Monte-Carlo precision goal defaults low (≈2 digits), per Mathematica. */
    double reltol = (o->prec_goal > 0.0) ? pow(10.0, -o->prec_goal) : 1e-2;
    double abstol = (o->acc_goal > 0.0) ? pow(10.0, -o->acc_goal) : 0.0;
    long maxp = (o->max_points > 0) ? o->max_points : 2000000;

    double _Complex val; double abserr;
    bool conv = mc_integrate_machine(ni_mc_sample, &ctx, a, b, d, quasi,
                                     abstol, reltol, maxp, &val, &abserr);

    for (size_t j = 0; j < d; j++) ni_bind_restore(&binds[j]);
    free(a); free(b); free(binds);

    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv)
        ni_warn("maxp", "Monte-Carlo did not reach the accuracy goal "
                "(error estimate %.3g)", abserr);
    return ni_from_complex_d(val);
}

/* ------------------------------------------------------------------ *
 *  Adaptive Genz-Malik cubature (low-to-moderate dimension)           *
 * ------------------------------------------------------------------ */

/* Try adaptive cubature over the box defined by `specs` (each {v, vmin, vmax}
 * with constant finite real bounds).  Returns NULL — so the caller falls back
 * to iterated 1-D quadrature — when the box is not a constant finite real
 * rectangle (e.g. a variable-dependent or infinite bound) or the integrand is
 * non-numeric there.  A single adaptive process over the whole box avoids the
 * multiplicative cost of nested 1-D adaptive quadrature, and the rule never
 * samples the box corners, so a singularity on a corner is handled cleanly. */
static Expr* ni_try_cubature(Expr* body, Expr** specs, size_t d, const NiOpts* o) {
    if (d < 2 || d > 15) return NULL;
    double* a = malloc(d * sizeof(double));
    double* b = malloc(d * sizeof(double));
    NiBind* binds = malloc(d * sizeof(NiBind));
    bool ok = true;
    for (size_t j = 0; j < d && ok; j++) {
        Expr* s = specs[j];
        if (s->type != EXPR_FUNCTION || s->data.function.arg_count != 3
            || s->data.function.args[0]->type != EXPR_SYMBOL) { ok = false; break; }
        Expr* lo = ni_num_endpoint(s->data.function.args[1]);
        Expr* hi = ni_num_endpoint(s->data.function.args[2]);
        if (!ni_to_double_real(lo, &a[j]) || !ni_to_double_real(hi, &b[j])) ok = false;
        expr_free(lo); expr_free(hi);
    }
    if (!ok) { free(a); free(b); free(binds); return NULL; }

    for (size_t j = 0; j < d; j++)
        ni_bind_snapshot(&binds[j], specs[j]->data.function.args[0]->data.symbol);

    NiMcCtx ctx;
    ctx.body = body; ctx.binds = binds; ctx.d = d;
    ctx.spec = numeric_machine_spec();

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    long max_eval = (o->max_points > 0) ? o->max_points : 0;

    double _Complex val; double abserr;
    CubStatus st = cub_integrate_machine(ni_mc_sample, &ctx, a, b, d,
                                         abstol, reltol, max_eval, &val, &abserr);

    for (size_t j = 0; j < d; j++) ni_bind_restore(&binds[j]);
    free(a); free(b); free(binds);

    if (st == CUB_NONNUMERIC) return NULL;   /* fall back to iterated quadrature */
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (st == CUB_NOCONV)
        ni_warn("ncvb", "cubature did not reach the accuracy goal "
                "(error estimate %.3g)", abserr);
    return ni_from_complex_d(val);
}

/* Route a single {x, xmin, xmax} integral by the kinds of its bounds. */
static Expr* ni_dispatch_1d(Expr* body, const char* var,
                            Expr* amin, Expr* amax, const NiOpts* o) {
    NiBound lo, hi;
    ni_classify(amin, &lo);
    ni_classify(amax, &hi);
    NiBoundKind ka = lo.kind, kb = hi.kind;

    if (ka == NI_B_BAD || kb == NI_B_BAD) return NULL;

    bool lo_fin = (ka == NI_B_FIN || ka == NI_B_CPLX);
    bool hi_fin = (kb == NI_B_FIN || kb == NI_B_CPLX);

    /* Both endpoints finite. */
    if (lo_fin && hi_fin) {
        if (ka == NI_B_FIN && kb == NI_B_FIN)
            return ni_run_1d_finite_real(body, var, creal(lo.z), creal(hi.z), o);
        /* at least one complex endpoint: straight-line contour segment */
        double _Complex nodes[2] = { lo.z, hi.z };
        return ni_run_contour(body, var, nodes, 2, o);
    }

    /* Complex infinite ray from a finite endpoint. */
    if (lo_fin && kb == NI_B_CRAY) return ni_run_ray(body, var, lo.z, hi.dir, 1.0, o);
    if (ka == NI_B_CRAY && hi_fin) return ni_run_ray(body, var, hi.z, lo.dir, -1.0, o);

    /* Real semi-infinite / doubly-infinite ranges. */
    if (lo_fin && kb == NI_B_PINF) return ni_run_1d_halfline(body, var, creal(lo.z), false, 1.0, o);
    if (ka == NI_B_NINF && hi_fin) return ni_run_1d_halfline(body, var, -creal(hi.z), true, 1.0, o);
    if (ka == NI_B_NINF && kb == NI_B_PINF) return ni_run_1d_wholeline(body, var, 1.0, o);
    /* reversed orientations negate the result */
    if (ka == NI_B_PINF && hi_fin) return ni_run_1d_halfline(body, var, creal(hi.z), false, -1.0, o);
    if (lo_fin && kb == NI_B_NINF) return ni_run_1d_halfline(body, var, -creal(lo.z), true, -1.0, o);
    if (ka == NI_B_PINF && kb == NI_B_NINF) return ni_run_1d_wholeline(body, var, -1.0, o);
    return NULL;   /* (+∞,+∞), (-∞,-∞), ray/ray, ray/real-∞: degenerate */
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ *
 *  Exclusions and Cauchy principal value                              *
 * ------------------------------------------------------------------ */

/* Append a real exclusion point in (a,b) to pts[]. */
static void ni_add_point(double v, double a, double b, double* pts, int* n, int maxn) {
    if (!(v > a && v < b)) return;
    for (int i = 0; i < *n; i++) if (fabs(pts[i] - v) < 1e-12 * (1.0 + fabs(v))) return;
    if (*n < maxn) pts[(*n)++] = v;
}

/* Extract real singular points of an Exclusions specification lying in (a,b):
 * a number, a List of numbers, or an equation lhs==rhs solved for `var`. */
static void ni_extract_points(Expr* excl, const char* var, double a, double b,
                              double* pts, int* n, int maxn) {
    if (!excl) return;
    double v;
    if (ni_to_double_real(excl, &v)) { ni_add_point(v, a, b, pts, n, maxn); return; }
    if (excl->type == EXPR_FUNCTION && excl->data.function.head->type == EXPR_SYMBOL) {
        const char* h = excl->data.function.head->data.symbol;
        if (h == SYM_List) {
            for (size_t i = 0; i < excl->data.function.arg_count; i++)
                ni_extract_points(excl->data.function.args[i], var, a, b, pts, n, maxn);
            return;
        }
        if (h == SYM_Equal) {
            /* Solve[lhs == rhs, var] and collect numeric roots in (a,b). */
            Expr** sv = malloc(2 * sizeof(Expr*));
            sv[0] = expr_copy(excl);
            sv[1] = expr_new_symbol(var);
            Expr* sol = eval_and_free(expr_new_function(expr_new_symbol("Solve"), sv, 2));
            free(sv);
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol == SYM_List) {
                for (size_t i = 0; i < sol->data.function.arg_count; i++) {
                    Expr* rule_set = sol->data.function.args[i];  /* {var -> val} */
                    if (rule_set->type != EXPR_FUNCTION) continue;
                    for (size_t j = 0; j < rule_set->data.function.arg_count; j++) {
                        Expr* r = rule_set->data.function.args[j];
                        if (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
                            && r->data.function.head->data.symbol == SYM_Rule
                            && r->data.function.arg_count == 2) {
                            Expr* val = numericalize(r->data.function.args[1], numeric_machine_spec());
                            double rv;
                            if (val && ni_to_double_real(val, &rv)) ni_add_point(rv, a, b, pts, n, maxn);
                            expr_free(val);
                        }
                    }
                }
            }
            expr_free(sol);
            return;
        }
    }
}

static int ni_dbl_cmp(const void* x, const void* y) {
    double dx = *(const double*)x, dy = *(const double*)y;
    return (dx < dy) ? -1 : (dx > dy) ? 1 : 0;
}

/* Evaluate the (complex) integrand at a real abscissa with the variable bound. */
static bool ni_eval_complex_at(Expr* body, NiBind* bind, NumericSpec spec,
                               double x, double _Complex* out) {
    Expr* xe = expr_new_real(x);
    ni_bind_set(bind, xe);
    expr_free(xe);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(body));
    if (!raw) return false;
    Expr* num = numericalize(raw, spec);
    expr_free(raw);
    if (!num) return false;
    bool ok = ni_to_complex(num, out);
    expr_free(num);
    return ok && isfinite(creal(*out)) && isfinite(cimag(*out));
}

typedef struct { Expr* body; NiBind* bind; NumericSpec spec; double c; } NiPvCtx;

/* Symmetric sample F(t) = f(c+t) + f(c-t): the simple-pole singularity of a
 * Cauchy principal value cancels between the mirrored points, leaving a finite
 * integrand whose integral over (0,d] equals the PV over [c-d, c+d]. */
static bool ni_pv_sym(void* vctx, double t, double _Complex* out) {
    NiPvCtx* p = (NiPvCtx*)vctx;
    double _Complex v1, v2;
    if (!ni_eval_complex_at(p->body, p->bind, p->spec, p->c + t, &v1)) return false;
    if (!ni_eval_complex_at(p->body, p->bind, p->spec, p->c - t, &v2)) return false;
    *out = v1 + v2;
    return true;
}

static bool ni_pv_plain(void* vctx, double x, double _Complex* out) {
    NiPvCtx* p = (NiPvCtx*)vctx;
    return ni_eval_complex_at(p->body, p->bind, p->spec, x, out);
}

/* Cauchy principal value of ∫_a^b f dx with simple poles at the sorted points
 * pts[0..npts-1] in (a,b).  Each pole's subinterval uses the symmetric mirror
 * rule; the remainder is ordinary adaptive quadrature. */
static Expr* ni_run_pv(Expr* body, const char* var, double a, double b,
                       double* pts, int npts, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiPvCtx ctx;
    ctx.body = body; ctx.bind = &bind; ctx.spec = numeric_machine_spec();

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    int max_subdiv; long max_eval;
    ni_gk_budget(o, &max_subdiv, &max_eval);

    /* Subinterval boundaries: a, midpoints between consecutive poles, b. */
    double _Complex total = 0.0;
    bool ok = true;
    for (int i = 0; i < npts; i++) {
        double L = (i == 0) ? a : 0.5 * (pts[i - 1] + pts[i]);
        double R = (i == npts - 1) ? b : 0.5 * (pts[i] + pts[i + 1]);
        double c = pts[i];
        double d = fmin(c - L, R - c);
        ctx.c = c;
        GkResult Rs = gk_integrate_machine(ni_pv_sym, &ctx, 0.0, d,
                                           abstol, reltol, max_subdiv, max_eval, true);
        if (Rs.status == GK_NONNUMERIC) { ok = false; break; }
        total += Rs.value;
        if (c - d > L) {
            GkResult Rl = gk_integrate_machine(ni_pv_plain, &ctx, L, c - d,
                                               abstol, reltol, max_subdiv, max_eval, true);
            total += Rl.value;
        }
        if (c + d < R) {
            GkResult Rr = gk_integrate_machine(ni_pv_plain, &ctx, c + d, R,
                                               abstol, reltol, max_subdiv, max_eval, true);
            total += Rr.value;
        }
    }
    ni_bind_restore(&bind);
    if (!ok) return NULL;
    return ni_from_complex_d(total);
}

Expr* builtin_nintegrate(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && ni_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    if (pos_end < 2) return NULL;                      /* need f + >= 1 spec */
    Expr* body = res->data.function.args[0];
    for (size_t i = 1; i < pos_end; i++)
        if (!ni_is_spec(res->data.function.args[i])) return NULL;

    NiOpts o;
    o.method = NI_AUTO;
    o.method_name = NULL;
    o.prec_mpfr = false; o.bits = 0; o.target_bits = 0;
    o.acc_goal = -1.0; o.prec_goal = -1.0;
    o.max_recursion = -1; o.min_recursion = 0; o.max_points = -1;
    o.exclusions = NULL;
    for (size_t i = pos_end; i < argc; i++)
        if (!ni_apply_option(res->data.function.args[i], &o)) return NULL;

    /* Run the MPFR double-exponential path at extra guard precision and round
     * the result back to the requested precision (the integrand can lose digits
     * to cancellation at off-node sample points, as in NSum). */
    if (o.prec_mpfr && o.bits > 0) {
        o.target_bits = o.bits;
        o.bits = o.target_bits + 64;
    }

    /* An explicitly requested but not-yet-implemented Method must not be
     * silently approximated: warn and stay unevaluated so the gap is visible. */
    if (!ni_method_implemented(o.method)) {
        ni_warn("method", "the method \"%s\" is not implemented; "
                "returning unevaluated", o.method_name ? o.method_name : "?");
        return NULL;
    }

    size_t nspecs = pos_end - 1;

    /* Vector/tensor integrand: thread NIntegrate over the List components, each
     * keeping the full spec/option list (handles {f1, f2, ...} and matrices). */
    if (body->type == EXPR_FUNCTION
        && body->data.function.head->type == EXPR_SYMBOL
        && body->data.function.head->data.symbol == SYM_List) {
        size_t nc = body->data.function.arg_count;
        Expr** comps = malloc(nc * sizeof(Expr*));
        for (size_t c = 0; c < nc; c++) {
            Expr** v = malloc(argc * sizeof(Expr*));   /* component + specs + opts */
            size_t w = 0;
            v[w++] = expr_copy(body->data.function.args[c]);
            for (size_t i = 1; i < argc; i++) v[w++] = expr_copy(res->data.function.args[i]);
            comps[c] = expr_new_function(expr_new_symbol("NIntegrate"), v, w);
            free(v);
        }
        Expr* listr = expr_new_function(expr_new_symbol("List"), comps, nc);
        free(comps);
        return listr;
    }

    /* Monte-Carlo: explicitly requested, or chosen automatically for region
     * (Boole / UnitStep) integrands and high-dimensional boxes where adaptive
     * cubature is impractical. */
    bool mc_quasi = (o.method == NI_QMC);
    bool want_mc = (o.method == NI_MC || o.method == NI_QMC || o.method == NI_AMC);
    if (o.method == NI_AUTO && (ni_body_has_region(body) || nspecs >= 6)) {
        want_mc = true; mc_quasi = true;   /* quasi-MC is the better default */
    }
    if (want_mc)
        return ni_run_mc(body, &res->data.function.args[1], nspecs, &o, mc_quasi);

    /* Multidimensional over a constant rectangular box: a single adaptive
     * Genz-Malik cubature is far cheaper than nested 1-D adaptive quadrature
     * (whose cost is multiplicative in dimension).  Falls through to the
     * iterated path when the box is not a constant finite real rectangle
     * (dependent / infinite / complex bounds) or at arbitrary precision. */
    if (o.method == NI_AUTO && nspecs >= 2 && !o.prec_mpfr && !o.exclusions) {
        Expr* cub = ni_try_cubature(body, &res->data.function.args[1], nspecs, &o);
        if (cub) return cub;
    }

    /* Multidimensional: iterated 1D quadrature.  The integrand of the outer
     * variable is an inner NIntegrate over the remaining specs (HoldAll plus
     * variable localisation lets a dependent inner bound, e.g.
     * {y, -Sqrt[1-x^2], Sqrt[1-x^2]}, see the bound outer index). */
    Expr* eff_body = body;
    Expr* inner = NULL;
    if (nspecs >= 2) {
        size_t ninner = 1 + (nspecs - 1) + (argc - pos_end);
        Expr** v = malloc(ninner * sizeof(Expr*));
        size_t w = 0;
        v[w++] = expr_copy(body);
        for (size_t i = 2; i < pos_end; i++) v[w++] = expr_copy(res->data.function.args[i]);
        for (size_t i = pos_end; i < argc; i++) v[w++] = expr_copy(res->data.function.args[i]);
        inner = expr_new_function(expr_new_symbol("NIntegrate"), v, w);
        free(v);
        eff_body = inner;
    }

    Expr* spec0 = res->data.function.args[1];
    size_t nnodes = spec0->data.function.arg_count - 1;   /* {var, n0, n1, ...} */
    const char* var = spec0->data.function.args[0]->data.symbol;

    Expr* out = NULL;
    if (nnodes == 2) {
        /* {x, xmin, xmax}: real finite / infinite / complex segment / ray. */
        Expr* amin = ni_num_endpoint(spec0->data.function.args[1]);
        Expr* amax = ni_num_endpoint(spec0->data.function.args[2]);
        double a, b;
        bool real_finite = ni_to_double_real(amin, &a) && ni_to_double_real(amax, &b);
        double pts[64]; int npts = 0;
        if (real_finite && (o.method == NI_PV || o.exclusions)) {
            double lo = a < b ? a : b, hi = a < b ? b : a;
            ni_extract_points(o.exclusions, var, lo, hi, pts, &npts, 64);
            qsort(pts, npts, sizeof(double), ni_dbl_cmp);
        }
        if (real_finite && o.method == NI_PV) {
            /* Cauchy principal value about the Exclusions poles. */
            out = (npts > 0) ? ni_run_pv(eff_body, var, a, b, pts, npts, &o)
                             : ni_run_1d_finite_real(eff_body, var, a, b, &o);
        } else if (real_finite && npts > 0) {
            /* Split the real interval at the excluded points (singularities are
             * then handled at the panel endpoints by the tanh-sinh rule). */
            double _Complex nodes[66];
            int nn = 0;
            nodes[nn++] = a;
            for (int i = 0; i < npts; i++) nodes[nn++] = pts[i];
            nodes[nn++] = b;
            out = ni_run_contour(eff_body, var, nodes, nn, &o);
        } else {
            out = ni_dispatch_1d(eff_body, var, amin, amax, &o);
        }
        expr_free(amin);
        expr_free(amax);
    } else {
        /* {x, x0, x1, ..., xk}: piecewise-linear contour / singularity split.
         * Every node must be a finite (real or complex) number. */
        double _Complex* nodes = malloc(nnodes * sizeof(double _Complex));
        bool ok = true;
        for (size_t i = 0; i < nnodes && ok; i++) {
            Expr* nv = ni_num_endpoint(spec0->data.function.args[i + 1]);
            NiBound b;
            ni_classify(nv, &b);
            if (b.kind == NI_B_FIN || b.kind == NI_B_CPLX) nodes[i] = b.z;
            else ok = false;
            expr_free(nv);
        }
        if (ok) out = ni_run_contour(eff_body, var, nodes, (int)nnodes, &o);
        free(nodes);
    }
    expr_free(inner);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nintegrate_init(void) {
    symtab_add_builtin("NIntegrate", builtin_nintegrate);
    /* HoldAll: the integrand and iterator specs must not be pre-evaluated;
     * the variable is Block-localised internally. Not Listable. */
    symtab_get_def("NIntegrate")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
