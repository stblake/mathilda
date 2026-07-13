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
#include "dequad.h"
#include "oscint.h"
#include "oscde.h"
#include "mcint.h"
#include "cubature.h"
#include "ncrule.h"
#include "levincoll.h"

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
    /* Suppressed while an enclosing numeric probe holds the arithmetic-warning
     * mute (see arith_warnings_mute_push).  In particular, iterated
     * multidimensional NIntegrate evaluates each outer sample by running an
     * *inner* NIntegrate over the remaining variables; the inner integral's
     * convergence diagnostics are sampling noise — only the outer integral's
     * verdict should reach the user. */
    if (arith_warnings_muted()) return;
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

/* Abscissa map applied before binding the integration variable. */
typedef enum {
    NI_MAP_AFFINE = 0,  /* x_scale·z + x_shift (real axis, reflection, contour) */
    NI_MAP_EXP_LO,      /* end_a + span·e^{-z}: maps z∈[0,∞) onto a singular a   */
    NI_MAP_EXP_HI       /* end_a − span·e^{-z}: maps z∈[0,∞) onto a singular b   */
} NiMapMode;

typedef struct {
    Expr*       body;    /* borrowed integrand (free in the variable)      */
    NiBind*     bind;
    NumericSpec spec;
    /* Abscissa map: the variable is bound to (x_scale·x + x_shift), used to
     * reflect a (-∞, b] range onto a half line and (later) to parametrise a
     * complex contour. */
    double _Complex x_scale;
    double _Complex x_shift;
    /* Exponential endpoint map (NI_MAP_EXP_*): tames a singular finite endpoint
     * by spreading it onto a half line.  The sample is multiplied by the map
     * Jacobian span·e^{-z}, which vanishes at the singular end. */
    NiMapMode   map_mode;
    double      end_a;   /* singular endpoint coordinate                   */
    double      span;    /* b − a > 0                                      */
} NiCtx;

static void ni_ctx_init(NiCtx* c, Expr* body, NiBind* bind, NumericSpec spec) {
    c->body = body; c->bind = bind; c->spec = spec;
    c->x_scale = 1.0; c->x_shift = 0.0;
    c->map_mode = NI_MAP_AFFINE; c->end_a = 0.0; c->span = 0.0;
}

/* Build the numeric leaf the variable is bound to from an affine-mapped
 * abscissa (real axis, reflection, or a complex contour segment). */
static Expr* ni_abscissa_expr(const NiCtx* c, double _Complex z) {
    double _Complex w = c->x_scale * z + c->x_shift;
    if (cimag(w) == 0.0) return expr_new_real(creal(w));
    return make_complex(expr_new_real(creal(w)), expr_new_real(cimag(w)));
}

/* Real abscissa of the exponential-endpoint map at half-line parameter z >= 0. */
static double ni_exp_abscissa(const NiCtx* c, double z) {
    double off = c->span * exp(-z);
    return (c->map_mode == NI_MAP_EXP_LO) ? c->end_a + off : c->end_a - off;
}

/* Evaluate the integrand with the variable bound to `value` (consumed).
 * Returns the numericalised result (caller frees) or NULL. */
static Expr* ni_eval_at(NiCtx* c, Expr* value) {
    ni_bind_set(c->bind, value);
    expr_free(value);
    eval_clock_bump();
    /* Mute arithmetic diagnostics (e.g. Power::infy 1/0) raised by evaluating
     * the integrand at or near a singular abscissa: such samples are detected
     * as non-finite below and discarded, so the message is pure noise. */
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(c->body));
    arith_warnings_mute_pop();
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

/* Machine sample callback: f at a (possibly mapped) real abscissa x. */
static bool ni_sample_machine(void* vctx, double x, double _Complex* out) {
    NiCtx* c = (NiCtx*)vctx;
    /* Exponential-endpoint map: the Jacobian d x/d z = span·e^{-z} vanishes at
     * the singular end, taming the integrand there.  Past the underflow horizon
     * (e^{-z} = 0) the abscissa would collapse exactly onto the singular endpoint
     * — evaluating the integrand there raises a spurious 1/0; the contribution is
     * negligible, so truncate the tail instead. */
    if (c->map_mode != NI_MAP_AFFINE) {
        double jac = c->span * exp(-x);
        double w = ni_exp_abscissa(c, x);
        /* Once the offset span·e^{-z} falls below the endpoint's rounding unit,
         * the abscissa collapses onto the singular endpoint itself — by
         * underflow when the endpoint is 0, or by catastrophic cancellation
         * (1 − tiny == 1) otherwise.  Evaluating the integrand there raises a
         * spurious 1/0; the tail past this horizon is negligible, so truncate. */
        if (w == c->end_a || !(jac > 0.0)) return false;
        Expr* num = ni_eval_at(c, expr_new_real(w));
        if (!num) return false;
        bool ok = ni_to_complex(num, out);
        expr_free(num);
        if (!ok) return false;
        *out *= jac;
        if (!isfinite(creal(*out)) || !isfinite(cimag(*out))) return false;
        return true;
    }
    Expr* num = ni_eval_at(c, ni_abscissa_expr(c, x));
    if (!num) return false;
    bool ok = ni_to_complex(num, out);
    expr_free(num);
    if (!ok) return false;
    if (!isfinite(creal(*out)) || !isfinite(cimag(*out))) return false;
    return true;
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
        if (e->data.symbol.name == SYM_Infinity)        { *dir = 1.0; return true; }
        if (e->data.symbol.name == SYM_ComplexInfinity) { *dir = 0.0; return true; }
        return false;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
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
    NI_MC, NI_QMC, NI_AMC, NI_PV, NI_OSCSING,
    NI_RIEMANN, NI_NEWTONCOTES,
    NI_UNIMPL          /* recognised name with no implementation yet            */
} NiMethod;

/* Methods implemented in the current build.  An explicitly requested method
 * outside this set causes NIntegrate to warn and stay unevaluated, so missing
 * implementations are visible rather than silently approximated.  Each phase
 * adds its method here as it lands. */
static bool ni_method_implemented(NiMethod m) {
    return m == NI_AUTO || m == NI_GK || m == NI_DE || m == NI_LEVIN
        || m == NI_MC || m == NI_QMC || m == NI_AMC || m == NI_PV
        || m == NI_OSCSING || m == NI_TRAP || m == NI_RIEMANN
        || m == NI_NEWTONCOTES;
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
    /* Fixed-rule sub-options (RiemannRule / TrapezoidalRule / NewtonCotesRule). */
    int    rule_type;       /* Riemann sampling: 0 Left, 1 Right, 2 Midpoint   */
    bool   romberg;         /* Romberg (Richardson) extrapolation; default on  */
    int    nc_points;       /* Newton–Cotes points per panel (2..5); default 3 */
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
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && ni_is_known_option(lhs->data.symbol.name);
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
    if (!strcmp(s, "RiemannRule")) return NI_RIEMANN;
    if (!strcmp(s, "NewtonCotesRule") || !strcmp(s, "NewtonCotes")) return NI_NEWTONCOTES;
    if (!strcmp(s, "LevinRule")) return NI_LEVIN;
    /* Exponential endpoint transform + integration-between-the-zeros for an
     * oscillatory endpoint singularity (∫_0^1 Cos[Log[x]/x]/x dx). */
    if (!strcmp(s, "OscillatorySingularity")) return NI_OSCSING;
    if (!strcmp(s, "MonteCarlo")) return NI_MC;
    if (!strcmp(s, "QuasiMonteCarlo")) return NI_QMC;
    if (!strcmp(s, "AdaptiveMonteCarlo") || !strcmp(s, "AdaptiveQuasiMonteCarlo")) return NI_AMC;
    if (!strcmp(s, "PrincipalValue")) return NI_PV;
    /* LocalAdaptive, ClenshawCurtisRule, LobattoKronrodRule, MultidimensionalRule,
     * CartesianRule, MultipanelRule, and any unknown name: not implemented as a
     * distinct rule. */
    return NI_UNIMPL;
}

static bool ni_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
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
        && (v->data.symbol.name == SYM_Infinity || v->data.symbol.name == SYM_Automatic)) {
        *out = -1.0; return true;
    }
    double d;
    if (ni_to_double_real(v, &d) && d > 0.0) { *out = d; return true; }
    return false;
}

static bool ni_parse_int(Expr* v, int* out, bool allow_auto) {
    if (allow_auto && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) {
        *out = -1; return true;
    }
    if (v->type == EXPR_INTEGER && v->data.integer >= 0) { *out = (int)v->data.integer; return true; }
    return false;
}

/* A fixed-rule sub-option "name" -> value inside Method -> {"rule", subopts...}.
 * Recognises "Type" (Riemann sampling), "RombergQuadrature" (Trapezoidal /
 * Newton–Cotes extrapolation) and "Points" (Newton–Cotes panel order).  Unknown
 * sub-options are accepted and ignored, matching the engine's lenient option
 * handling. */
static void ni_apply_method_subopt(Expr* sub, NiOpts* o) {
    if (!sub || sub->type != EXPR_FUNCTION) return;
    Expr* h = sub->data.function.head;
    if (h->type != EXPR_SYMBOL) return;
    if ((h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed)
        || sub->data.function.arg_count != 2) return;
    Expr* k = sub->data.function.args[0];
    Expr* v = sub->data.function.args[1];
    if (k->type != EXPR_STRING) return;
    const char* key = k->data.string;

    if (!strcmp(key, "Type") && v->type == EXPR_STRING) {
        if      (!strcmp(v->data.string, "Left"))     o->rule_type = 0;
        else if (!strcmp(v->data.string, "Right"))    o->rule_type = 1;
        else if (!strcmp(v->data.string, "Midpoint")) o->rule_type = 2;
    } else if (!strcmp(key, "RombergQuadrature") && v->type == EXPR_SYMBOL) {
        if      (v->data.symbol.name == SYM_True)  o->romberg = true;
        else if (v->data.symbol.name == SYM_False) o->romberg = false;
    } else if (!strcmp(key, "Points") && v->type == EXPR_INTEGER) {
        long pts = v->data.integer;
        if (pts < 2) pts = 2;
        if (pts > 5) pts = 5;
        o->nc_points = (int)pts;
    }
}

static bool ni_apply_option(Expr* rule, NiOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_Method) {
        Expr* m = rhs;
        /* Method -> {"name", subopts...} or Method -> "name" or symbol. */
        if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL
            && m->data.function.head->data.symbol.name == SYM_List
            && m->data.function.arg_count >= 1) {
            for (size_t i = 1; i < m->data.function.arg_count; i++)
                ni_apply_method_subopt(m->data.function.args[i], o);
            m = m->data.function.args[0];
        }
        if (m->type == EXPR_STRING) {
            o->method_name = m->data.string;
            o->method = ni_method_from_string(m->data.string);
            return true;
        }
        if (m->type == EXPR_SYMBOL) {
            if (m->data.symbol.name == SYM_Automatic) { o->method = NI_AUTO; o->method_name = NULL; return true; }
            o->method_name = m->data.symbol.name;
            o->method = ni_method_from_string(m->data.symbol.name);
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
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) { o->max_points = -1; return true; }
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
        || e->data.function.head->data.symbol.name != SYM_List) return false;
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

/* ------------------------------------------------------------------ *
 *  Levin collocation: kernel detection + driver                       *
 * ------------------------------------------------------------------ */

/* True if the symbol named `name` occurs anywhere in e. */
static bool levin_occurs(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type == EXPR_FUNCTION) {
        if (levin_occurs(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (levin_occurs(e->data.function.args[i], name)) return true;
    }
    return false;
}

static bool levin_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool levin_is_head(const Expr* e, const char* head) {
    return e && e->type == EXPR_FUNCTION
        && levin_is_sym(e->data.function.head, head);
}

/* Extract the real phase h from an exponent of the form I·h (the argument of
 * Exp[...] or the exponent of Power[E, ...]).  Returns an owned phase Expr, or
 * NULL when the exponent is not a pure-imaginary multiple of an I-free factor
 * (a real part would make the kernel decay/grow, not oscillate). */
static Expr* levin_phase_from_exponent(const Expr* z) {
    if (!levin_is_head(z, SYM_Times)) return NULL;   /* need I·h */
    size_t n = z->data.function.arg_count;
    int icount = 0;
    for (size_t i = 0; i < n; i++)
        if (levin_is_sym(z->data.function.args[i], SYM_I)) icount++;
    if (icount != 1) return NULL;
    Expr** rest = malloc(sizeof(Expr*) * n);
    if (!rest) return NULL;
    size_t r = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = z->data.function.args[i];
        if (levin_is_sym(f, SYM_I)) continue;
        if (levin_occurs(f, SYM_I)) {   /* I nested in another factor */
            for (size_t k = 0; k < r; k++) expr_free(rest[k]);
            free(rest); return NULL;
        }
        rest[r++] = expr_copy(f);
    }
    Expr* h;
    if (r == 0)      { free(rest); return NULL; }    /* exponent was just I */
    else if (r == 1) { h = rest[0]; free(rest); }
    else             { h = expr_new_function(expr_new_symbol(SYM_Times), rest, r);
                       free(rest); }
    return h;
}

/* Recognise one oscillatory factor: Cos[h], Sin[h], Exp[I·h] or Power[E, I·h].
 * On success sets *k and returns the owned phase *g. */
static bool levin_match_osc_factor(const Expr* f, LevinKernel* k, Expr** g) {
    if (levin_is_head(f, SYM_Cos) && f->data.function.arg_count == 1) {
        *k = LEVIN_KERNEL_COS; *g = expr_copy(f->data.function.args[0]); return true;
    }
    if (levin_is_head(f, SYM_Sin) && f->data.function.arg_count == 1) {
        *k = LEVIN_KERNEL_SIN; *g = expr_copy(f->data.function.args[0]); return true;
    }
    if (levin_is_head(f, SYM_Exp) && f->data.function.arg_count == 1) {
        Expr* h = levin_phase_from_exponent(f->data.function.args[0]);
        if (!h) return false;
        *k = LEVIN_KERNEL_EXP; *g = h; return true;
    }
    if (levin_is_head(f, SYM_Power) && f->data.function.arg_count == 2
        && levin_is_sym(f->data.function.args[0], SYM_E)) {
        Expr* h = levin_phase_from_exponent(f->data.function.args[1]);
        if (!h) return false;
        *k = LEVIN_KERNEL_EXP; *g = h; return true;
    }
    return false;
}

/* Decompose `body` into amp·{cos g | sin g | e^{ig}} with exactly one
 * oscillatory factor whose phase contains `var`.  On success returns owned
 * *g_out (phase) and *amp_out (amplitude, literal 1 if none).  Returns false —
 * leaving outputs untouched — when no/multiple oscillatory factors are present
 * or the phase is constant in `var`. */
static bool levin_detect_kernel(const Expr* body, const char* var,
                                LevinKernel* kind, Expr** g_out, Expr** amp_out) {
    Expr* const single[1] = { (Expr*)body };
    Expr* const* fs;
    size_t nf;
    if (levin_is_head(body, SYM_Times)) {
        fs = body->data.function.args;
        nf = body->data.function.arg_count;
    } else { fs = single; nf = 1; }

    /* The oscillatory factor must oscillate in *this* axis: its phase contains
     * `var`.  A trig/exp factor whose phase is free of `var` (a constant, or an
     * oscillation in another integration variable, as in a separable product
     * Sin[1/x] Cos[1000 y]) is part of the amplitude for this axis. */
    int osc_idx = -1;
    LevinKernel k = LEVIN_KERNEL_EXP;
    Expr* g = NULL;
    for (size_t i = 0; i < nf; i++) {
        LevinKernel kk; Expr* gg = NULL;
        if (levin_match_osc_factor(fs[i], &kk, &gg)) {
            if (!levin_occurs(gg, var)) { expr_free(gg); continue; }  /* amplitude */
            if (osc_idx >= 0) { expr_free(gg); expr_free(g); return false; }
            osc_idx = (int)i; k = kk; g = gg;
        }
    }
    if (osc_idx < 0) return false;

    /* amp = product of the non-oscillatory factors (or the literal 1). */
    Expr* amp;
    if (nf == 1) {
        amp = expr_new_integer(1);
    } else {
        Expr** rest = malloc(sizeof(Expr*) * (nf - 1));
        size_t r = 0;
        for (size_t i = 0; i < nf; i++)
            if ((int)i != osc_idx) rest[r++] = expr_copy(fs[i]);
        amp = (r == 1) ? rest[0]
                       : expr_new_function(expr_new_symbol(SYM_Times), rest, r);
        free(rest);
    }
    *kind = k; *g_out = g; *amp_out = amp;
    return true;
}

/* g' = D[g, var], evaluated symbolically.  Returns owned g' or NULL when the
 * derivative cannot be taken in closed form (an unevaluated D/Derivative head
 * survives), so the caller can fall back. */
static Expr* levin_phase_derivative(const Expr* g, const char* var) {
    Expr** dargs = malloc(sizeof(Expr*) * 2);
    if (!dargs) return NULL;
    dargs[0] = expr_copy((Expr*)g);
    dargs[1] = expr_new_symbol(var);
    Expr* dcall = expr_new_function(expr_new_symbol(SYM_D), dargs, 2);
    free(dargs);
    /* The derivative must be taken with `var` FREE.  In the Automatic cascade a
     * prior quadrature rule may have left `var` bound to its last sample value,
     * which would collapse D[g, var] to 0 (or worse).  Detach any current
     * binding for the duration of the symbolic differentiation, then restore it
     * (the per-node samplers rebind it themselves). */
    SymbolDef* vd = symtab_get_def(var);
    Rule* saved_own = vd->own_values;
    vd->own_values = NULL;
    eval_clock_bump();
    Expr* gp = eval_and_free(dcall);
    vd->own_values = saved_own;
    eval_clock_bump();
    if (!gp) return NULL;
    if (levin_occurs(gp, SYM_D) || levin_occurs(gp, "Derivative")) {
        expr_free(gp); return NULL;
    }
    return gp;
}

#ifdef USE_MPFR
/* Recognise a well-aligned Fourier-type integrand F(x) = amp(x)·{sin|cos}(ω x)
 * whose oscillation has frequency ω (a positive real constant) and zero phase
 * offset, so the Ooura–Mori DE nodes land on the oscillation's zeros kπ/ω (sine)
 * or (k+½)π/ω (cosine).  On success sets *kind and *omega.  Anything else — a
 * complex exponential kernel, a non-constant or offset phase, a non-numeric
 * frequency — returns false so the caller uses the general (slower) method. */
static bool ni_fourier_detect(const Expr* body, const char* var,
                              OscDeKind* kind, double* omega) {
    LevinKernel k; Expr* g = NULL; Expr* amp = NULL;
    /* The held integrand is unevaluated, so its product is nested
     * (Times[Times[x,Sin[x]],Power[…]]); flatten a copy so the single-level
     * factor scan in levin_detect_kernel sees the trig factor. */
    Expr* fb = expr_copy((Expr*)body);
    eval_flatten_args(fb, SYM_Times);
    bool ok = levin_detect_kernel(fb, var, &k, &g, &amp);
    expr_free(fb);
    if (!ok) return false;
    expr_free(amp);
    if (k != LEVIN_KERNEL_SIN && k != LEVIN_KERNEL_COS) { expr_free(g); return false; }

    /* ω = d g / d var must be a nonzero real constant (linear phase). */
    Expr* gp = levin_phase_derivative(g, var);
    if (!gp || levin_occurs(gp, var)) { expr_free(g); expr_free(gp); return false; }
    double w;
    bool okw = ni_to_double_real(gp, &w);
    expr_free(gp);
    if (!okw || !(w != 0.0)) { expr_free(g); return false; }

    /* Phase offset g|_{var=0} must vanish (mod π handling is out of scope; a
     * misaligned phase falls back rather than risk a wrong "converged" answer). */
    Expr* rargs[2] = { expr_new_symbol(var), expr_new_integer(0) };
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    Expr* raargs[2] = { expr_copy(g), rule };
    Expr* g0 = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_ReplaceAll), raargs, 2));
    expr_free(g);
    double off;
    bool oko = g0 && ni_to_double_real(g0, &off);
    expr_free(g0);
    if (!oko || fabs(off) > 1e-9) return false;

    *omega = fabs(w);
    *kind = (k == LEVIN_KERNEL_SIN) ? OSCDE_SIN : OSCDE_COS;
    return true;
}
#endif /* USE_MPFR */

/* Levin collocation on a finite real interval.  Detects the oscillatory kernel
 * in ctx->body, differentiates the phase symbolically, and runs the collocation
 * engine through the existing numeric sampler.  Restricted to the plain real
 * axis; returns have=false (caller falls back) on any non-Levin input, a
 * singular phase endpoint, or weak oscillation. */
static NiAtt ni_try_levin_collocation(NiCtx* ctx, double a, double b,
                                      const NiOpts* o, double reltol) {
    (void)o;
    NiAtt none = { false, false, 0.0, INFINITY };
    if (ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return none;

    const char* var = ctx->bind->name;
    LevinKernel kind; Expr* g = NULL; Expr* amp = NULL;
    if (!levin_detect_kernel(ctx->body, var, &kind, &g, &amp)) return none;

    Expr* gp = levin_phase_derivative(g, var);
    if (!gp) { expr_free(g); expr_free(amp); return none; }

    /* Three sampler contexts sharing the binding/precision but evaluating the
     * amplitude, the phase, and its derivative respectively. */
    NiCtx camp = *ctx; camp.body = amp;
    NiCtx cg   = *ctx; cg.body   = g;
    NiCtx cgp  = *ctx; cgp.body  = gp;

    /* Spot-check the symbolic g' against a central difference of g at the
     * midpoint — a cheap guard against a mis-evaluated derivative. */
    double xm = 0.5 * (a + b), hh = 1e-6 * fmax(1.0, fabs(b - a));
    double _Complex gpm, gphi, gplo;
    if (ni_sample_machine(&cgp, xm, &gpm)
        && ni_sample_machine(&cg, xm + hh, &gphi)
        && ni_sample_machine(&cg, xm - hh, &gplo)) {
        double _Complex fd = (gphi - gplo) / (2.0 * hh);
        if (cabs(fd - gpm) > 1e-3 * (1.0 + cabs(gpm))) {
            expr_free(g); expr_free(amp); expr_free(gp); return none;
        }
    }

    LevinResult r = levin_collocation_machine(
        a, b, ni_sample_machine, &camp, ni_sample_machine, &cgp,
        ni_sample_machine, &cg, kind, reltol, 64);

    expr_free(g); expr_free(amp); expr_free(gp);
    NiAtt out;
    out.have = r.have; out.conv = r.conv; out.val = r.val; out.err = r.err;
    return out;
}

#ifdef USE_MPFR
/* Arbitrary-precision Levin collocation on a finite real interval.  `ctx` must
 * already carry the MPFR numeric spec.  Writes the result into (re,im) and the
 * convergence verdict into *conv; returns true if an estimate was produced. */
static bool ni_try_levin_collocation_mpfr(NiCtx* ctx, double a, double b,
                                          const NiOpts* o, mpfr_t re, mpfr_t im,
                                          bool* conv) {
    if (ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return false;
    const char* var = ctx->bind->name;
    LevinKernel kind; Expr* g = NULL; Expr* amp = NULL;
    if (!levin_detect_kernel(ctx->body, var, &kind, &g, &amp)) return false;
    Expr* gp = levin_phase_derivative(g, var);
    if (!gp) { expr_free(g); expr_free(amp); return false; }

    NiCtx camp = *ctx; camp.body = amp;
    NiCtx cg   = *ctx; cg.body   = g;
    NiCtx cgp  = *ctx; cgp.body  = gp;
    bool ok = levin_collocation_mpfr(
        a, b, o->bits, ni_sample_mpfr, &camp, ni_sample_mpfr, &cgp,
        ni_sample_mpfr, &cg, kind, ni_mpfr_reltol(o), 64, re, im, conv);
    expr_free(g); expr_free(amp); expr_free(gp);
    return ok;
}
#endif

/* Map MaxRecursion / MaxPoints to a fixed-rule level and evaluation budget.
 * Romberg-accelerated rules (Trapezoidal / Newton–Cotes) converge in a handful
 * of levels; the un-extrapolated Riemann sums need far more panels, so they get
 * more levels but a bounded evaluation count to stop an over-tight PrecisionGoal
 * from grinding indefinitely. */
static void ni_ncr_budget(const NiOpts* o, int* max_levels, long* max_eval) {
    int r = o->max_recursion;
    if (o->method == NI_RIEMANN) {
        *max_levels = (r < 0) ? 24 : (r > 30 ? 30 : (r < 6 ? 6 : r));
        *max_eval   = (o->max_points > 0) ? o->max_points : 4000000;
    } else {
        *max_levels = (r < 0) ? 20 : (r > 30 ? 30 : (r < 4 ? 4 : r));
        *max_eval   = (o->max_points > 0) ? o->max_points : 0;
    }
}

/* Equally-spaced composite rule (RiemannRule / TrapezoidalRule /
 * NewtonCotesRule) over a finite interval. */
static NiAtt ni_try_ncr(NiCtx* ctx, double a, double b, const NiOpts* o,
                        double reltol, double abstol) {
    NcrRuleKind kind;
    if      (o->method == NI_TRAP)        kind = NCR_TRAPEZOIDAL;
    else if (o->method == NI_NEWTONCOTES) kind = NCR_NEWTONCOTES;
    else kind = (o->rule_type == 1) ? NCR_RIEMANN_RIGHT
              : (o->rule_type == 2) ? NCR_RIEMANN_MIDPOINT
                                    : NCR_RIEMANN_LEFT;
    int max_levels; long max_eval;
    ni_ncr_budget(o, &max_levels, &max_eval);
    double _Complex v; double err;
    bool conv = ncr_integrate_machine(ni_sample_machine, ctx, a, b, kind,
                                      o->nc_points, o->romberg, reltol, abstol,
                                      max_levels, max_eval, &v, &err);
    NiAtt out;
    out.have = isfinite(creal(v)) && isfinite(cimag(v));
    out.conv = conv;
    out.val = v;
    out.err = err;
    return out;
}

/* Probe the integrand approaching one endpoint (`side` +1 → a, −1 → b) over a
 * geometric ladder of distances.  An endpoint is treated as singular when a
 * probe is non-numeric there, or the integrand's magnitude both becomes large
 * in absolute terms and grows far beyond its value a tenth of the way in — the
 * envelope of an algebraic / one-over-x or oscillatory singularity.  Requiring a
 * large absolute magnitude keeps a merely large relative swing of a bounded
 * oscillation (whose outer probe happened to land near a zero) from registering. */
static bool ni_endpoint_singular(NiCtx* c, double a, double b, int side) {
    double span = b - a;
    double mref = 0.0, mmax = 0.0;
    for (int k = 1; k <= 8; k++) {
        double d = span * pow(10.0, -(double)k);
        double x = (side > 0) ? a + d : b - d;
        double _Complex v;
        if (!ni_sample_machine(c, x, &v)) return true;   /* non-numeric at endpoint */
        double m = cabs(v);
        if (k == 1) mref = m;
        if (m > mmax) mmax = m;
    }
    return mmax > 1.0e3 && mmax > 1.0e2 * fmax(1.0, mref);
}

/* Handle a singular endpoint of the real interval [a,b] by the exponential
 * coordinate map  x = a + (b−a)e^{-t}  (or its mirror for a singular b), which
 * carries the singularity to t → ∞ where the map Jacobian (b−a)e^{-t} damps it,
 * then integrates the resulting half line.  An oscillatory singularity yields a
 * non-decaying oscillation on the half line, which the exp-sinh rule cannot
 * resolve but the integrate-between-the-zeros + Wynn extrapolation can.  This is
 * the engine behind Method -> "OscillatorySingularity".  When `force` is false
 * only endpoints detected singular are transformed; when true (explicit method)
 * both endpoints are tried and the better-converged result kept. */
static NiAtt ni_try_endpoint_sing(NiCtx* ctx, double a, double b,
                                  const NiOpts* o, bool force) {
    NiAtt best = { false, false, 0.0, INFINITY };
    /* Only meaningful on the plain real axis (not an affine contour segment). */
    if (!(b > a) || ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return best;

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    /* Transform only an endpoint detected singular: the exponential map damps a
     * singularity at the z → ∞ end only, so applying it to a regular endpoint
     * (the wrong mirror) would instead sample the *other*, singular end at z = 0
     * and converge to the wrong value.  With an explicit method (`force`) and no
     * endpoint detected singular, fall back to trying both mirrors. */
    int sides[2]; int ns = 0;
    if (ni_endpoint_singular(ctx, a, b, +1)) sides[ns++] = +1;
    if (ni_endpoint_singular(ctx, a, b, -1)) sides[ns++] = -1;
    if (ns == 0 && force) { sides[ns++] = +1; sides[ns++] = -1; }

    for (int s = 0; s < ns; s++) {
        ctx->span = b - a;
        if (sides[s] > 0) { ctx->map_mode = NI_MAP_EXP_LO; ctx->end_a = a; }
        else              { ctx->map_mode = NI_MAP_EXP_HI; ctx->end_a = b; }
        double _Complex v; double err;
        bool conv = dequad_halfline_machine(ni_sample_machine, ctx, 0.0, reltol,
                                            ni_de_levels(o), &v, &err);
        if (!conv) {
            double _Complex ov; double oerr;
            if (osc_integrate_machine(ni_sample_machine, ctx, 0.0, 0.0, true,
                                      reltol, 600, &ov, &oerr)) {
                v = ov; err = oerr; conv = true;
            }
        }
        ctx->map_mode = NI_MAP_AFFINE;
        NiAtt at;
        at.have = isfinite(creal(v)) && isfinite(cimag(v));
        at.conv = conv; at.val = v; at.err = err;
        ni_consider(&best, at);
    }
    return best;
}

/* Integrate over the finite real parameter interval [a,b] with the variable
 * already configured in `ctx` (directly, or through an affine abscissa map for
 * a complex segment).  Gauss-Kronrod samples strictly interior abscissae, so it
 * never evaluates a singular endpoint (no spurious 1/0 messages); for Automatic
 * it falls back to the endpoint-robust tanh-sinh rule, then — at a detected
 * oscillatory endpoint singularity — to the exponential-endpoint transform, and
 * finally to between-the-zeros panel quadrature for a many-period integrand. */
static NiAtt ni_core_finite(NiCtx* ctx, double a, double b, const NiOpts* o) {
    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    NiAtt best = { false, false, 0.0, INFINITY };
    if (o->method == NI_LEVIN) {
        /* Genuine Levin collocation for a detected f·{cos g|sin g|e^{ig}}
         * kernel.  When collocation does not apply or does not converge (a
         * phase singular at an endpoint, e.g. Sin[1/x] at 0; weak oscillation;
         * a non-Levin kernel) fall back through the same cascade Automatic
         * uses: adaptive Gauss-Kronrod, the endpoint-robust tanh-sinh rule, the
         * exponential endpoint-singularity transform (only for a *detected*
         * singular end — forcing the wrong mirror would converge to a wrong
         * value), and finally between-the-zeros panel quadrature. */
        ni_consider(&best, ni_try_levin_collocation(ctx, a, b, o, reltol));
        if (!best.conv) ni_consider(&best, ni_try_gk(ctx, a, b, o, reltol, abstol));
        if (!best.conv) ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
        if (!best.conv) ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, false));
        if (!best.conv) ni_consider(&best, ni_try_osc(ctx, a, b, reltol));
    } else if (o->method == NI_DE) {
        ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
    } else if (o->method == NI_OSCSING) {
        ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, true));
    } else if (o->method == NI_TRAP || o->method == NI_RIEMANN
               || o->method == NI_NEWTONCOTES) {
        ni_consider(&best, ni_try_ncr(ctx, a, b, o, reltol, abstol));
    } else {
        ni_consider(&best, ni_try_gk(ctx, a, b, o, reltol, abstol));
        if (!best.conv && (o->method == NI_AUTO || !best.have))
            ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
        /* Levin collocation: when the smooth/tanh-sinh rules fail, a detected
         * f·{cos|sin|e^{i}} kernel oscillating too fast for Gauss-Kronrod is
         * resolved exactly by collocation (self-gated; returns nothing if the
         * integrand is not of Levin form or the oscillation is too weak). */
        if (!best.conv && o->method == NI_AUTO)
            ni_consider(&best, ni_try_levin_collocation(ctx, a, b, o, reltol));
        /* Oscillatory endpoint singularity: spread the singular endpoint onto a
         * half line by the exponential map and integrate between the zeros. */
        if (!best.conv && o->method == NI_AUTO)
            ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, false));
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
        bool conv;
        /* An explicit fixed rule is honoured at high precision too (the default
         * MPFR path is the endpoint-robust tanh-sinh rule). */
        if (o->method == NI_TRAP || o->method == NI_RIEMANN
            || o->method == NI_NEWTONCOTES) {
            NcrRuleKind kind = (o->method == NI_TRAP) ? NCR_TRAPEZOIDAL
                : (o->method == NI_NEWTONCOTES) ? NCR_NEWTONCOTES
                : (o->rule_type == 1) ? NCR_RIEMANN_RIGHT
                : (o->rule_type == 2) ? NCR_RIEMANN_MIDPOINT : NCR_RIEMANN_LEFT;
            int max_levels; long max_eval;
            ni_ncr_budget(o, &max_levels, &max_eval);
            conv = ncr_integrate_mpfr(ni_sample_mpfr, &ctx, am, bm, kind,
                                      o->nc_points, o->romberg, bits,
                                      ni_mpfr_reltol(o), 0.0, max_levels, max_eval,
                                      re, im, &abserr);
        } else {
            /* Levin collocation at high precision for a detected oscillatory
             * kernel (explicit LevinRule, or Automatic when it converges);
             * otherwise the endpoint-robust tanh-sinh rule. */
            bool got = false; conv = false; abserr = 0.0;
            if (o->method == NI_LEVIN || o->method == NI_AUTO) {
                bool lconv = false;
                if (ni_try_levin_collocation_mpfr(&ctx, a, b, o, re, im, &lconv)
                    && (o->method == NI_LEVIN || lconv)) {
                    conv = lconv; got = true;
                }
            }
            if (!got)
                conv = denint_tanhsinh_mpfr(ni_sample_mpfr, &ctx, am, bm, bits,
                                            ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                            re, im, &abserr);
        }
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
        double mreltol = ni_mpfr_reltol(o);
        /* Geometry ctx: the zero-finding for the oscillatory fallback needs only
         * a few correct digits, so it samples at machine precision — a shallow
         * copy of ctx with the machine spec (borrowed pointers, never freed). */
        NiCtx gctx = ctx;
        gctx.spec = numeric_machine_spec();
        bool conv = false;

        /* Fast path: a well-aligned Fourier integrand ∫_0^∞ amp·{sin|cos}(ωx) dx
         * is done by the Ooura–Mori double-exponential rule, whose nodes fall on
         * the oscillation's zeros — a few thousand samples reach hundreds of
         * digits, where integrate-between-the-zeros + Wynn needs one costly
         * high-precision panel per half-period. */
        if ((o->method == NI_AUTO || o->method == NI_LEVIN)
            && !reflect && sign > 0.0 && a == 0.0) {
            OscDeKind fkind; double fomega;
            if (ni_fourier_detect(body, var, &fkind, &fomega))
                conv = oscde_fourier_mpfr(ni_sample_mpfr, &ctx, fomega, fkind,
                                          bits, mreltol, 16, re, im, &abserr);
        }

        if (!conv && o->method == NI_LEVIN) {
            conv = osc_integrate_mpfr(ni_sample_machine, &gctx, ni_sample_mpfr, &ctx,
                                      a, bits, mreltol, 600, ni_mpfr_levels(o),
                                      re, im, &abserr);
        } else if (!conv) {
            conv = dequad_halfline_mpfr(ni_sample_mpfr, &ctx, am, bits,
                                        mreltol, ni_mpfr_levels(o),
                                        re, im, &abserr);
            /* Oscillatory fallback: a slowly-decaying oscillatory tail
             * (x Sin[x]/(x^2+a), Bessel functions, …) defeats exp-sinh at any
             * precision.  Integrate between the zeros and accelerate the partial
             * sums with the MPFR Wynn epsilon — the high-precision counterpart of
             * the machine path's fallback. */
            if (!conv && o->method == NI_AUTO) {
                mpfr_t ore, oim;
                mpfr_init2(ore, bits); mpfr_init2(oim, bits);
                double oerr;
                if (osc_integrate_mpfr(ni_sample_machine, &gctx, ni_sample_mpfr, &ctx,
                                       a, bits, mreltol, 600, ni_mpfr_levels(o),
                                       ore, oim, &oerr)) {
                    mpfr_set(re, ore, MPFR_RNDN); mpfr_set(im, oim, MPFR_RNDN);
                    abserr = oerr; conv = true;
                }
                mpfr_clears(ore, oim, (mpfr_ptr)0);
            }
        }
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

    /* Oscillatory fallback: a non-decaying oscillatory integrand (Exp[I x^2],
     * Cos[x^2], Sin[x] …) defeats sinh-sinh, which presumes the tails decay and
     * otherwise diverges into garbage.  Split the line at 0 — integrating the
     * left half via the reflection x = -u — and integrate each half between the
     * zeros of the oscillation, accelerating the partial sums with Wynn's
     * epsilon (mirrors the half-line path).  If neither half registers as
     * oscillatory, osc_integrate reports failure and the sinh-sinh result
     * stands. */
    if (!conv && o->method == NI_AUTO) {
        double _Complex lv, rv; double lerr, rerr;
        ctx.x_scale = -1.0;   /* (-∞, 0] via x = -u, u in [0, ∞) */
        bool lc = osc_integrate_machine(ni_sample_machine, &ctx, 0.0, 0.0, true,
                                        reltol, 600, &lv, &lerr);
        ctx.x_scale = 1.0;    /* [0, ∞) */
        bool rc = osc_integrate_machine(ni_sample_machine, &ctx, 0.0, 0.0, true,
                                        reltol, 600, &rv, &rerr);
        if (lc && rc) { val = lv + rv; abserr = lerr + rerr; conv = true; }
    }

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
        const char* h = e->data.function.head->data.symbol.name;
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
    arith_warnings_mute_push();   /* see ni_eval_at: integrand-message noise */
    Expr* raw = eval_and_free(expr_copy(c->body));
    arith_warnings_mute_pop();
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
        ni_bind_snapshot(&binds[j], specs[j]->data.function.args[0]->data.symbol.name);

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
        ni_bind_snapshot(&binds[j], specs[j]->data.function.args[0]->data.symbol.name);

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

/* ------------------------------------------------------------------ *
 *  Multivariate Levin (dimension reduction over a rectangular box)    *
 * ------------------------------------------------------------------ */

/* Levin's multivariate method reduces a d-dimensional oscillatory integral one
 * axis at a time.  Reducing the oscillatory axis x_k by collocation turns the
 * inner integral into its boundary term  ∫ f e^{ig} dx_k = [p e^{ig}], so the
 * remaining (d-1)-dimensional integral has the *same* phase restricted to the
 * faces and a new amplitude p — crucially evaluated by a single small linear
 * solve, NOT a nested adaptive quadrature.  This avoids the combinatorial blow
 * up of iterated quadrature when an axis oscillates rapidly (every outer sample
 * would otherwise spawn a full inner integral).
 *
 * This implementation covers the 2-D case: one axis is reduced by 1-D Levin and
 * the remaining axis is integrated by the ordinary 1-D machine cascade, each of
 * whose samples costs one inner Levin solve.  Higher dimensions fall through to
 * the iterated path. */

/* Bind `b`'s variable to a real value.  ni_bind_set copies its argument but does
 * NOT take ownership, so the value Expr must be freed by the caller. */
static void nd_bind_real(NiBind* b, double v) {
    Expr* e = expr_new_real(v);
    ni_bind_set(b, e);
    expr_free(e);
}

/* Outer-axis sampler: bind the outer variable to xo and return the inner
 * integral over the reduction axis, evaluated by 1-D Levin collocation. */
typedef struct {
    NiBind* outer_bind;
    double  a_red, b_red;
    NiCtx*  camp; NiCtx* cg; NiCtx* cgp;   /* reduction-axis amp / g / g' */
    LevinKernel kind;
    double  reltol;
    LevinPrep* prep;   /* non-NULL: factored-once fast path (g' const in outer) */
} NdLevinCtx;

static bool nd_levin_outer_sample(void* vctx, double xo, double _Complex* out) {
    NdLevinCtx* c = (NdLevinCtx*)vctx;
    nd_bind_real(c->outer_bind, xo);
    if (c->prep) {
        /* Fast path: the collocation matrix is fixed; one back-substitution. */
        if (!levin_prepared_solve(c->prep, ni_sample_machine, c->camp,
                                  ni_sample_machine, c->cg, c->kind, out))
            return false;
        return isfinite(creal(*out)) && isfinite(cimag(*out));
    }
    LevinResult r = levin_collocation_machine(
        c->a_red, c->b_red,
        ni_sample_machine, c->camp,
        ni_sample_machine, c->cgp,
        ni_sample_machine, c->cg,
        c->kind, c->reltol, 64);
    if (!r.have || !r.conv) return false;
    if (!isfinite(creal(r.val)) || !isfinite(cimag(r.val))) return false;
    *out = r.val;
    return true;
}

/* True if the phase derivative g' (sampled at the reduction endpoints/midpoint)
 * is independent of the outer variable — the condition under which the
 * collocation matrix can be factored once and re-used across outer samples. */
static bool nd_gprime_const_in_outer(NiCtx* cgp, NiBind* ob,
                                     double a_red, double b_red,
                                     double a_out, double b_out) {
    double rs[3] = { a_red, 0.5 * (a_red + b_red), b_red };
    double os[2] = { a_out + 0.25 * (b_out - a_out),
                     a_out + 0.75 * (b_out - a_out) };
    double _Complex v0[3];
    nd_bind_real(ob, os[0]);
    for (int i = 0; i < 3; i++)
        if (!ni_sample_machine(cgp, rs[i], &v0[i])) return false;
    nd_bind_real(ob, os[1]);
    for (int i = 0; i < 3; i++) {
        double _Complex v1;
        if (!ni_sample_machine(cgp, rs[i], &v1)) return false;
        if (cabs(v1 - v0[i]) > 1e-9 * (1.0 + cabs(v0[i]))) return false;
    }
    return true;
}

/* Test whether axis `var_red` of `body` is a good reduction axis: a single
 * oscillatory factor in var_red whose phase and derivative stay finite over the
 * reduction interval for a few interior values of the outer variable, and whose
 * inner Levin solve converges at the box centre.  On success returns the owned
 * kernel pieces. */
static bool nd_axis_reducible(Expr* body, const char* var_red, const char* var_out,
                              double a_red, double b_red, double a_out, double b_out,
                              const NiOpts* o, LevinKernel* kind,
                              Expr** g_out, Expr** amp_out, Expr** gp_out) {
    LevinKernel k; Expr* g = NULL; Expr* amp = NULL;
    if (!levin_detect_kernel(body, var_red, &k, &g, &amp)) return false;
    Expr* gp = levin_phase_derivative(g, var_red);
    if (!gp) { expr_free(g); expr_free(amp); return false; }

    NiBind ob; ni_bind_snapshot(&ob, var_out);
    NiBind rb; ni_bind_snapshot(&rb, var_red);
    NiCtx cg, cgp; ni_ctx_init(&cg, g, &rb, numeric_machine_spec());
    ni_ctx_init(&cgp, gp, &rb, numeric_machine_spec());

    /* Phase regularity: g and g' finite at the reduction endpoints/midpoint for
     * a few interior outer abscissae. */
    bool ok = true;
    for (int io = 1; io <= 3 && ok; io++) {
        double xo = a_out + (b_out - a_out) * io / 4.0;
        nd_bind_real(&ob, xo);
        double rs[3] = { a_red, 0.5 * (a_red + b_red), b_red };
        for (int ir = 0; ir < 3 && ok; ir++) {
            double _Complex gv, gpv;
            if (!ni_sample_machine(&cg, rs[ir], &gv)
                || !ni_sample_machine(&cgp, rs[ir], &gpv)) ok = false;
        }
    }
    /* Inner Levin must actually converge at the box centre. */
    if (ok) {
        NiCtx camp; ni_ctx_init(&camp, amp, &rb, numeric_machine_spec());
        nd_bind_real(&ob, 0.5 * (a_out + b_out));
        LevinResult r = levin_collocation_machine(
            a_red, b_red, ni_sample_machine, &camp, ni_sample_machine, &cgp,
            ni_sample_machine, &cg, k, o->prec_goal > 0 ? 0.0 : 1e-9, 64);
        if (!r.have || !r.conv) ok = false;
    }
    ni_bind_restore(&rb);
    ni_bind_restore(&ob);

    if (!ok) { expr_free(g); expr_free(amp); expr_free(gp); return false; }
    *kind = k; *g_out = g; *amp_out = amp; *gp_out = gp;
    return true;
}

/* Two-dimensional Levin: reduce one axis by collocation, integrate the other by
 * the 1-D machine cascade.  Returns the result, or NULL to fall back. */
static Expr* ni_try_levin_nd(Expr* body, Expr** specs, size_t nspecs,
                             const NiOpts* o) {
    if (nspecs != 2 || o->prec_mpfr || o->exclusions) return NULL;

    /* Both specs must be a constant finite real interval {v, a, b}. */
    const char* var[2]; double lo[2], hi[2];
    for (int s = 0; s < 2; s++) {
        Expr* sp = specs[s];
        if (sp->data.function.arg_count != 3) return NULL;
        var[s] = sp->data.function.args[0]->data.symbol.name;
        Expr* amin = ni_num_endpoint(sp->data.function.args[1]);
        Expr* amax = ni_num_endpoint(sp->data.function.args[2]);
        bool okb = amin && amax && ni_to_double_real(amin, &lo[s])
                && ni_to_double_real(amax, &hi[s]);
        expr_free(amin); expr_free(amax);
        if (!okb || !(hi[s] > lo[s])) return NULL;
    }

    /* Choose the reduction axis: try each, prefer the first that qualifies. */
    for (int rsel = 0; rsel < 2; rsel++) {
        int osel = 1 - rsel;
        LevinKernel kind; Expr* g = NULL; Expr* amp = NULL; Expr* gp = NULL;
        if (!nd_axis_reducible(body, var[rsel], var[osel],
                               lo[rsel], hi[rsel], lo[osel], hi[osel],
                               o, &kind, &g, &amp, &gp))
            continue;

        /* Drive the outer (osel) axis with the reduction sampler. */
        NiBind ob; ni_bind_snapshot(&ob, var[osel]);
        NiBind rb; ni_bind_snapshot(&rb, var[rsel]);
        NiCtx camp, cg, cgp;
        ni_ctx_init(&camp, amp, &rb, numeric_machine_spec());
        ni_ctx_init(&cg,   g,   &rb, numeric_machine_spec());
        ni_ctx_init(&cgp,  gp,  &rb, numeric_machine_spec());
        double in_reltol = (o->prec_goal > 0) ? 0.0 : 1e-9;
        NdLevinCtx nd = { &ob, lo[rsel], hi[rsel], &camp, &cg, &cgp, kind,
                          in_reltol, NULL };

        /* If g' does not depend on the outer variable, factor the collocation
         * matrix once: every outer sample is then a cheap back-substitution
         * rather than a full O(n^3) solve, which is what makes a rapidly
         * oscillatory outer axis (many samples) tractable. */
        if (nd_gprime_const_in_outer(&cgp, &ob, lo[rsel], hi[rsel],
                                     lo[osel], hi[osel])) {
            double xc = 0.5 * (lo[osel] + hi[osel]);
            nd_bind_real(&ob, xc);
            LevinResult ref = levin_collocation_machine(
                lo[rsel], hi[rsel], ni_sample_machine, &camp, ni_sample_machine,
                &cgp, ni_sample_machine, &cg, kind, in_reltol, 64);
            for (int n = 8; n <= 64 && ref.have; n *= 2) {
                LevinPrep* P = levin_prepare_machine(lo[rsel], hi[rsel], n,
                                                     ni_sample_machine, &cgp);
                if (!P) continue;
                double _Complex tv;
                if (levin_prepared_solve(P, ni_sample_machine, &camp,
                                         ni_sample_machine, &cg, kind, &tv)
                    && cabs(tv - ref.val) <= 1e-8 * (1.0 + cabs(ref.val))) {
                    nd.prep = P; break;
                }
                levin_prepare_free(P);
            }
        }

        double reltol, abstol; ni_machine_tols(o, &reltol, &abstol);
        int max_subdiv; long max_eval; ni_gk_budget(o, &max_subdiv, &max_eval);
        /* Bound the outer evaluation count so a hard outer axis (an oscillatory
         * endpoint singularity such as Sin[1/x], whose zeros cluster without
         * bound) returns a best estimate in finite time rather than grinding.
         * The cached path is far cheaper per sample, so it gets a higher cap. */
        long cap = nd.prep ? 50000 : 6000;
        if (max_eval <= 0 || max_eval > cap) max_eval = cap;
        double _Complex val = 0.0; double err = INFINITY; bool conv = false, have = false;

        GkResult R = gk_integrate_machine(nd_levin_outer_sample, &nd,
                                          lo[osel], hi[osel], abstol, reltol,
                                          max_subdiv, max_eval, true);
        if (R.status != GK_NONNUMERIC) {
            val = R.value; err = R.abs_err; have = true; conv = (R.status == GK_OK);
        }
        if (!conv) {
            double _Complex v2; double e2;
            if (denint_tanhsinh_machine(nd_levin_outer_sample, &nd, lo[osel],
                                        hi[osel], reltol, ni_de_levels(o), &v2, &e2)
                && isfinite(creal(v2)) && isfinite(cimag(v2))) {
                val = v2; err = e2; have = true; conv = true;
            }
        }
        /* No between-the-zeros fallback here: each of its panels would trigger
         * an inner Levin solve, and for an outer axis whose phase is singular at
         * an endpoint (e.g. Sin[1/x] at 0) the zeros cluster without bound — a
         * combinatorial blow-up.  GK/tanh-sinh are budget-bounded; if they do
         * not converge we report the best estimate with a warning rather than
         * grinding indefinitely. */

        ni_bind_restore(&rb);
        ni_bind_restore(&ob);
        if (nd.prep) levin_prepare_free(nd.prep);
        expr_free(g); expr_free(amp); expr_free(gp);

        if (!have || !isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
        if (!conv)
            ni_warn("ncvb", "multivariate Levin did not reach the accuracy goal "
                    "(error estimate %.3g)", err);
        return ni_from_complex_d(val);
    }
    return NULL;   /* no reducible axis — fall back to iterated quadrature */
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
        const char* h = excl->data.function.head->data.symbol.name;
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
            Expr* sol = eval_and_free(expr_new_function(expr_new_symbol(SYM_Solve), sv, 2));
            free(sv);
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == SYM_List) {
                for (size_t i = 0; i < sol->data.function.arg_count; i++) {
                    Expr* rule_set = sol->data.function.args[i];  /* {var -> val} */
                    if (rule_set->type != EXPR_FUNCTION) continue;
                    for (size_t j = 0; j < rule_set->data.function.arg_count; j++) {
                        Expr* r = rule_set->data.function.args[j];
                        if (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
                            && r->data.function.head->data.symbol.name == SYM_Rule
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
    arith_warnings_mute_push();   /* see ni_eval_at: integrand-message noise */
    Expr* raw = eval_and_free(expr_copy(body));
    arith_warnings_mute_pop();
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
    o.rule_type = 0; o.romberg = true; o.nc_points = 3;
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
        && body->data.function.head->data.symbol.name == SYM_List) {
        size_t nc = body->data.function.arg_count;
        Expr** comps = malloc(nc * sizeof(Expr*));
        for (size_t c = 0; c < nc; c++) {
            Expr** v = malloc(argc * sizeof(Expr*));   /* component + specs + opts */
            size_t w = 0;
            v[w++] = expr_copy(body->data.function.args[c]);
            for (size_t i = 1; i < argc; i++) v[w++] = expr_copy(res->data.function.args[i]);
            comps[c] = expr_new_function(expr_new_symbol(SYM_NIntegrate), v, w);
            free(v);
        }
        Expr* listr = expr_new_function(expr_new_symbol(SYM_List), comps, nc);
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
    /* Multivariate Levin: reduce an oscillatory axis by collocation so the
     * remaining integral costs one inner solve per outer sample instead of a
     * nested adaptive quadrature.  Requested explicitly, or chosen by Automatic
     * (ni_try_levin_nd self-gates: it returns NULL unless an axis carries a
     * detectable, regular oscillatory kernel that converges), before the more
     * general cubature / iterated paths. */
    if ((o.method == NI_LEVIN || o.method == NI_AUTO)
        && nspecs == 2 && !o.prec_mpfr && !o.exclusions) {
        Expr* lv = ni_try_levin_nd(body, &res->data.function.args[1], nspecs, &o);
        if (lv) return lv;
    }

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
        inner = expr_new_function(expr_new_symbol(SYM_NIntegrate), v, w);
        free(v);
        eff_body = inner;
    }

    Expr* spec0 = res->data.function.args[1];
    size_t nnodes = spec0->data.function.arg_count - 1;   /* {var, n0, n1, ...} */
    const char* var = spec0->data.function.args[0]->data.symbol.name;

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
    symtab_set_docstring("NIntegrate",
        "NIntegrate[f, {x, a, b}] numerically integrates f over a..b. Bounds may "
        "be real, complex (a straight-line/contour) or infinite. Options: Method "
        "(\"GlobalAdaptive\", \"DoubleExponential\", \"TrapezoidalRule\", "
        "\"RiemannRule\", \"NewtonCotesRule\", \"LevinRule\", "
        "\"OscillatorySingularity\", \"MonteCarlo\", \"PrincipalValue\", ...), "
        "WorkingPrecision, PrecisionGoal, AccuracyGoal, MaxRecursion, MaxPoints, "
        "Exclusions. Fixed rules take Method sub-options \"Type\" "
        "(\"Left\"/\"Right\"/\"Midpoint\"), \"RombergQuadrature\" and \"Points\".");
}
