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
#include "fem.h"
#include "modular.h"   /* mth_block_enter/leave, mth_bind_push/pop: unwind-safe localisation */
#include "gkadapt.h"
#include "denint.h"
#include "dequad.h"
#include "nc_phase.h"      /* linear-phase groups, Bessel -> Hankel */
#include "oscint.h"
#include "oscde.h"
#include "mcint.h"
#include "cubature.h"
#include "ncrule.h"
#include "levincoll.h"
#include "compile/autocompile.h"

#include <complex.h>
#include <fenv.h>
#include <float.h>
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
#include "core.h"     /* leaf_count_internal */
#include "common.h"   /* common_scan_inexact */
#include "rationalize.h"  /* internal_contains_inexact: exact range endpoints */
#include "message.h"  /* mth_msg_suppress_push, mth_msg_fired_count */
#include "mtime.h"    /* mth_mono_seconds: the oscillatory-endpoint probe budget */
#include "eval.h"
#include "nc_accuracy.h"  /* shared AccuracyGoal/PrecisionGoal handling */
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
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(g_arith_warnings_muted, "NIntegrate", tag, fmt, ap);
    va_end(ap);
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

/* ------------------------------------------------------------------ *
 *  Escalation past the double range                                  *
 * ------------------------------------------------------------------ *
 *
 * A machine-precision integrand can leave the double range: x^x^x reaches
 * 6^46656 = 2.7e36305 on [1, 6], and a Laurent polynomial of degree 7000 on
 * the unit circle reaches 1e1428.  Machine arithmetic then hands back either an
 * arbitrary-precision Real beyond 1.8e308 or, for complex arithmetic,
 * Overflow[].  Neither converts to a double, so the sampler declined the
 * point; a quadrature rule that drops its largest samples returns whatever the
 * rest add up to (0. for the circle integral) or nothing at all (x^x^x stayed
 * unevaluated).
 *
 * The machine sampler instead records that a sample left the double range, the
 * machine pass is discarded with its diagnostics, and the integral is re-run on
 * the MPFR path at machine-equivalent precision (53 bits plus the usual guard),
 * whose numbers have MPFR's exponent range.  The result is a 53-bit Real, the
 * same representation machine arithmetic itself uses beyond the double range.
 *
 * g_ni_beyond is set by any machine sample that left the range; while
 * g_ni_escalated is false the goal diagnostics of such a pass are withheld.
 * builtin_nintegrate saves and restores both, so nested (iterated) integrals
 * keep their own verdicts. */
static bool g_ni_beyond = false;
/* Fewest bits carried by an MPFR-pass sample value (0: none seen).  An
 * integrand that returns machine numbers inside a WorkingPrecision run cannot
 * deliver the requested digits; builtin_nintegrate reports it (precw). */
static long g_ni_sample_bits = 0;
static bool g_ni_escalated = false;
/* Absolute error the endpoint cancellation guard could not remove (0: none).
 * A deep sample whose raised-precision re-evaluation is capped by an inexact
 * constant of the integrand (a 45-digit 1/2 in (E^x - 1 - x - x^2/2)/x^3 at
 * WorkingPrecision -> 40) keeps the noise its cancellation left; its size
 * times the node's distance to the endpoint (the order of its quadrature
 * weight) is added here, and builtin_nintegrate lowers the precision of the
 * result to what that error leaves (precw).  Saved and restored like
 * g_ni_sample_bits. */
static double g_ni_guard_noise_l2 = -INFINITY;   /* log2 of the sum */
/* Fewest bits an EXACT range endpoint numericalised to inside a
 * WorkingPrecision pass (0: every endpoint reached the working precision).
 * Set by ni_endpoint_mpfr for an exact spec node whose numerical value at the
 * working precision came back with fewer bits (a head with only a machine
 * implementation); builtin_nintegrate then rounds the result to those bits and
 * reports it (precw), rather than labelling an integral over a machine
 * endpoint with the working precision.  An endpoint given inexact (1.1) is
 * taken as its binary value, as N takes it, and is not counted.  Saved and
 * restored like g_ni_sample_bits. */
static long g_ni_endpoint_bits = 0;

/* True when the numericalised sample `e` is a number outside the double range
 * (an MPFR value beyond 2^1024) or carries Overflow[] anywhere. */
static bool ni_beyond_double(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
#ifdef USE_MPFR
        case EXPR_MPFR:
            return mpfr_number_p(e->data.mpfr) && !mpfr_zero_p(e->data.mpfr)
                && mpfr_get_exp(e->data.mpfr) > DBL_MAX_EXP;
#endif
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Overflow)
                return true;
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (ni_beyond_double(e->data.function.args[i])) return true;
            return false;
        }
        default: return false;
    }
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
    int         token;     /* binding frame (modular.h mth_bind_push) */
    bool        valid;
} NiBind;

static void ni_bind_snapshot(NiBind* b, const char* name) {
    b->name = name;
    /* Unwind-safe: a TimeConstrained timeout drains the frame. */
    b->token = mth_bind_push_name(name, MTH_BIND_ATTRS);
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
    mth_bind_pop(b->token);   /* frees the temporary value */
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
    /* Auto-compile fast path (machine-precision, real-abscissa integrals only).
     * Opt-in per entry point via ac_enabled; compiled lazily on first sample and
     * freed by the caller after ni_bind_restore. */
    bool          ac_enabled;
    bool          ac_tried;
    AutoCompiled* ac;
    /* Compiled MPFR fast path for a high-WorkingPrecision integral: the compiled
     * body evaluates f[x] at the working precision, avoiding the per-sample
     * tree-walk.  Same lazy-compile / faithful-degrade contract as `ac`. */
    AutoCompiled* ac_prec;
    /* The sampling-parameter interval of the current rule (±Inf when open),
     * used to tell an out-of-range sample in the bulk of the interval from one
     * at an endpoint singularity (ni_note_beyond). */
    double ep_lo, ep_hi;
#ifdef USE_MPFR
    /* Endpoint cancellation guard of the MPFR finite-interval rules (see
     * ni_sample_mpfr_guarded): the exact endpoints of the range being
     * integrated (borrowed, valid while gd_on), and per endpoint whether its
     * deep samples were found free of cancellation (1), cancelling at
     * gd_rate lost bits per bit of distance (2), losing a depth-independent
     * gd_rate bits (4), probed once at depth gd_L1 losing gd_lost1 bits (3),
     * or not yet probed (0). */
    bool        gd_on;
    mpfr_srcptr gd_a, gd_b;
    /* The exact integrand the guard's raised-precision evaluations use in
     * place of `body` (borrowed; NULL to use `body`), see NiOpts.pre_body. */
    Expr*       gd_exact;
    int         gd_state[2];
    double      gd_rate[2];
    double      gd_L1[2], gd_lost1[2];
#endif
} NiCtx;

static void ni_ctx_init(NiCtx* c, Expr* body, NiBind* bind, NumericSpec spec) {
    c->body = body; c->bind = bind; c->spec = spec;
    c->x_scale = 1.0; c->x_shift = 0.0;
    c->map_mode = NI_MAP_AFFINE; c->end_a = 0.0; c->span = 0.0;
    c->ac_enabled = false; c->ac_tried = false; c->ac = NULL; c->ac_prec = NULL;
    c->ep_lo = -INFINITY; c->ep_hi = INFINITY;
#ifdef USE_MPFR
    c->gd_on = false; c->gd_a = NULL; c->gd_b = NULL; c->gd_exact = NULL;
    c->gd_state[0] = c->gd_state[1] = 0;
    c->gd_rate[0] = c->gd_rate[1] = 0.0;
    c->gd_L1[0] = c->gd_L1[1] = c->gd_lost1[0] = c->gd_lost1[1] = 0.0;
#endif
}

/* Point a (value-copied) context at a different body.  The lazily-compiled `ac`
 * belongs to the OLD body, so it must be dropped — otherwise a copy made after
 * the parent already compiled would evaluate the new body with the wrong
 * program.  The derived sub-body keeps the interpreter path (ac stays off), so
 * there is nothing to free on the copy. */
static void ni_ctx_rebody(NiCtx* c, Expr* body) {
    c->body = body;
    c->ac_enabled = false; c->ac_tried = false; c->ac = NULL; c->ac_prec = NULL;
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
static void ni_ac_prepare(NiCtx* c) {
    if (c->ac_enabled && !c->ac_tried) {
        c->ac_tried = true;
        Expr* vsym = expr_new_symbol(c->bind->name);
        const Expr* vs[1] = { vsym };
        if (c->spec.mode == NUMERIC_MODE_MACHINE) {
            c->ac = autocompile_new(c->body, vs, 1);
        }
#ifdef USE_MPFR
        else if (c->spec.mode == NUMERIC_MODE_MPFR && c->spec.bits > 0) {
            c->ac_prec = autocompile_new_prec(c->body, vs, 1, (long)c->spec.bits);
            /* a complex-valued integrand of a real variable (a contour
             * parametrised by its angle) compiles in the complex domain */
            if (!c->ac_prec)
                c->ac_prec = autocompile_new_prec_z(c->body, vs, 1, (long)c->spec.bits);
        }
#endif
        expr_free(vsym);
    }
}

static Expr* ni_eval_at(NiCtx* c, Expr* value) {
    /* Compiled fast path (opt-in, machine precision only — the MPFR sampler also
     * routes through here and must keep full precision).  Compiled lazily so the
     * spec.mode is already settled.  Only a REAL abscissa uses it; a complex
     * (contour) abscissa or a non-finite compiled result — e.g. the interpreter
     * would produce a complex value from Sqrt of a negative — falls through to
     * the interpreter, which supplies that contribution. */
    ni_ac_prepare(c);
#ifdef USE_MPFR
    /* Compiled MPFR fast path: the abscissa `value` is already an EXPR_MPFR, and
     * the compiled body returns f[value] at the working precision — or NULL to
     * interpret this point (non-finite, or the interpreter would go complex where
     * a real was expected), exactly the machine path's contract below. */
    /* The compiled program loads its argument at the working precision; an
     * abscissa carrying more bits (the power map next to a non-zero endpoint)
     * would lose exactly the bits that hold its distance to the endpoint, so
     * it is interpreted. */
    if (c->ac_prec && !(value->type == EXPR_MPFR
                        && (long)mpfr_get_prec(value->data.mpfr) > (long)c->spec.bits)) {
        Expr* r = autocompiled_eval_mpfr(c->ac_prec, (const Expr* const*)&value);
        if (r) { expr_free(value); return r; }
    }
#endif
    if (c->ac) {
        double xin;
        if (ni_to_double_real(value, &xin)) {
            double _Complex z;
            if (autocompiled_eval_complex(c->ac, &xin, &z)
                && isfinite(creal(z)) && isfinite(cimag(z))) {
                expr_free(value);
                return (cimag(z) == 0.0)
                    ? expr_new_real(creal(z))
                    : make_complex(expr_new_real(creal(z)), expr_new_real(cimag(z)));
            }
        }
    }
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
#ifdef USE_MPFR
    /* The precision the integrand itself delivered, before numericalize
     * dresses a machine number up at the working precision. */
    if (c->spec.mode == NUMERIC_MODE_MPFR) {
        long nb = 0;
        Expr *nre, *nim;
        const Expr* parts[2] = { raw, NULL };
        if (is_complex(raw, &nre, &nim)) { parts[0] = nre; parts[1] = nim; }
        for (int i = 0; i < 2 && parts[i]; i++) {
            long pb = parts[i]->type == EXPR_REAL ? 53
                    : parts[i]->type == EXPR_MPFR ? (long)mpfr_get_prec(parts[i]->data.mpfr) : 0;
            if (pb && (!nb || pb < nb)) nb = pb;
        }
        if (nb && (!g_ni_sample_bits || nb < g_ni_sample_bits)) g_ni_sample_bits = nb;
    }
#endif
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

/* The compiled machine path WITHOUT the Expr round-trip.
 *
 * The general path builds an Expr for the abscissa, and ni_eval_at unboxes it
 * straight back to a double for the compiled body, boxes the result, and the
 * caller unboxes that again -- two allocations and two frees per sample point,
 * spent only to pass doubles through an Expr-shaped interface. On an
 * oscillatory integrand, which samples thousands of points per call, that
 * churn is most of the cost: expr_new_real, expr_free and _free dominate the
 * profile of NIntegrate[Sin[200 x], {x, 0, Pi}].
 *
 * Declining (false) leaves the caller on the general path, so every fallback
 * the interpreter provides is preserved unchanged -- in particular a non-finite
 * compiled result still routes to the interpreter, which supplies that
 * contribution (e.g. where the compiled body would go complex).
 *
 * c->ac is NULL until ni_eval_at has lazily built it, so the first sample of a
 * call takes the general path (building it) and the rest are fast. */
static bool ni_sample_compiled(NiCtx* c, double w, double _Complex* out) {
    if (!c->ac) return false;
    double _Complex z;
    if (!autocompiled_eval_complex(c->ac, &w, &z)) return false;
    if (!isfinite(creal(z)) || !isfinite(cimag(z))) return false;
    *out = z;
    return true;
}

/* Record a machine sample at parameter x that left the double range (see
 * g_ni_beyond).  A sample within 1e-10 of the span from an endpoint does not
 * count: that is an endpoint singularity (Cos[Log[x]/x]/x reaches 1e308 only
 * for x < 1e-308), whose samples the rules drop by design and whose weight is
 * negligible; the exponential endpoint map exists for exactly that case. */
static void ni_note_beyond(const NiCtx* c, const Expr* num, double x) {
    if (c->spec.mode != NUMERIC_MODE_MACHINE || c->map_mode != NI_MAP_AFFINE)
        return;
    if (isfinite(c->ep_lo) && isfinite(c->ep_hi)) {
        double tol = 1e-10 * (c->ep_hi - c->ep_lo);
        if (x - c->ep_lo <= tol || c->ep_hi - x <= tol) return;
    }
    if (ni_beyond_double(num)) g_ni_beyond = true;
}

/* Machine sample callback: f at a (possibly mapped) real abscissa x. */
static bool ni_sample_machine(void* vctx, double x, double _Complex* out) {
    NiCtx* c = (NiCtx*)vctx;
    /* Once a sample has left the double range the machine pass will be
     * discarded (see g_ni_beyond): decline every further sample so the rule
     * cascade unwinds at once instead of running each rule to completion
     * (x^x^x spent 11000 machine samples before the MPFR pass began). */
    if (g_ni_beyond && !g_ni_escalated && c->spec.mode == NUMERIC_MODE_MACHINE)
        return false;
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
        if (ni_sample_compiled(c, w, out)) {
            *out *= jac;
            return isfinite(creal(*out)) && isfinite(cimag(*out));
        }
        Expr* num = ni_eval_at(c, expr_new_real(w));
        if (!num) return false;
        bool ok = ni_to_complex(num, out);
        if (!ok || !isfinite(creal(*out)) || !isfinite(cimag(*out)))
            ni_note_beyond(c, num, x);
        expr_free(num);
        if (!ok) return false;
        *out *= jac;
        if (!isfinite(creal(*out)) || !isfinite(cimag(*out))) return false;
        return true;
    }
    /* Same affine map ni_abscissa_expr applies; a real abscissa is the only one
     * the compiled body accepts, and a complex contour point falls through. */
    double _Complex wz = c->x_scale * x + c->x_shift;
    if (cimag(wz) == 0.0 && ni_sample_compiled(c, creal(wz), out)) return true;

    Expr* num = ni_eval_at(c, ni_abscissa_expr(c, x));
    if (!num) return false;
    bool ok = ni_to_complex(num, out);
    if (!ok || !isfinite(creal(*out)) || !isfinite(cimag(*out)))
        ni_note_beyond(c, num, x);
    expr_free(num);
    if (!ok) return false;
    if (!isfinite(creal(*out)) || !isfinite(cimag(*out))) return false;
    return true;
}

#ifdef USE_MPFR
/* `e` with every arbitrary-precision number widened exactly to `bits` (a new
 * reference; NULL when `e` holds none, so the caller keeps `e`).  The value
 * of each number is unchanged: only the precision it caps the arithmetic at.
 * Nothing under a holding head is touched. */
static Expr* ni_widen_numbers(const Expr* e, long bits) {
    if (!e) return NULL;
    if (e->type == EXPR_MPFR) {
        if ((long)mpfr_get_prec(e->data.mpfr) >= bits) return NULL;
        mpfr_t w;
        mpfr_init2(w, (mpfr_prec_t)bits);
        mpfr_set(w, e->data.mpfr, MPFR_RNDN);          /* exact: more bits */
        return expr_new_mpfr_move(w);
    }
    if (e->type != EXPR_FUNCTION) return NULL;
    Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        uint32_t at = symtab_get_def(h->data.symbol.name)->attributes;
        if (at & (ATTR_HOLDALL | ATTR_HOLDFIRST | ATTR_HOLDREST | ATTR_HOLDALLCOMPLETE))
            return NULL;
    }
    size_t n = e->data.function.arg_count;
    Expr** a = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* y = ni_widen_numbers(e->data.function.args[i], bits);
        if (!y) continue;
        if (!a) {
            a = malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < n; j++) a[j] = NULL;
        }
        a[i] = y;
    }
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++)
        if (!a[i]) a[i] = expr_copy(e->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy(h), a, n);
    free(a);
    return out;
}

/* MPFR sample callback: f at the real abscissa x (mapped by the real part of
 * the affine map — MPFR is used for real integration ranges).  (out_re,out_im)
 * are pre-initialised at the working precision by the DE kernel. */
/* f at the abscissa xm with every number in the evaluation carried at `prec`
 * bits (the abscissa is widened exactly; the interpreter computes at the
 * precision of its operands).  The compiled program, built for the working
 * precision, is bypassed.  (re, im) are initialised by the caller. */
static bool ni_eval_mpfr_at_prec(NiCtx* c, const mpfr_t xm, long prec, mpfr_t re, mpfr_t im,
                                 long* delivered) {
    if ((long)mpfr_get_prec(xm) > prec) prec = (long)mpfr_get_prec(xm);
    mpfr_t xw;
    mpfr_init2(xw, (mpfr_prec_t)prec);
    mpfr_set(xw, xm, MPFR_RNDN);                      /* exact: prec >= prec(xm) */
    Expr* xe = expr_new_mpfr_copy(xw);
    mpfr_clear(xw);
    NumericSpec saved = c->spec;
    AutoCompiled* saved_ac = c->ac_prec;
    bool saved_en = c->ac_enabled;
    Expr* saved_body = c->body;
    c->spec.bits = prec;
    c->ac_prec = NULL;
    c->ac_enabled = false;
    if (c->gd_exact) c->body = c->gd_exact;
    long sbits = g_ni_sample_bits;
    g_ni_sample_bits = 0;
    /* A raised-precision re-evaluation is an internal probe: the user asked for
     * f at the working precision, not at p + 16384 bits, so what it reports is
     * sampling noise -- neither printed nor seen by Check[]. A node a tanh-sinh
     * rule placed 2^-200 inside an endpoint known only to the working
     * precision (a FindRoot root) can lie OUTSIDE the integrand's domain once
     * the abscissa is carried exactly: Surd[-y + Sqrt[-Log[g[y]]], 3] then
     * meets a complex argument and said Surd::preal, with a 160-digit number,
     * on every probe -- thousands of lines, half the run time. */
    unsigned long fired0 = mth_msg_fired_count();
    mth_msg_suppress_push();
    Expr* num = ni_eval_at(c, xe);
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired0);
    if (delivered) *delivered = g_ni_sample_bits;   /* 0: not an MPFR value */
    if (g_ni_sample_bits == 0 || (sbits && sbits < g_ni_sample_bits)) g_ni_sample_bits = sbits;
    c->spec = saved;
    c->ac_prec = saved_ac;
    c->ac_enabled = saved_en;
    c->body = saved_body;
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, re, im, &inexact);
    if (ok && (!mpfr_number_p(re) || !mpfr_number_p(im))) ok = false;
    expr_free(num);
    return ok;
}

/* log2 of the relative difference |u - v|/|v| of two complex MPFR values
 * (-INFINITY when equal; 0 when v is zero and u is not). */
static double ni_mpfr_rel_log2(const mpfr_t ur, const mpfr_t ui, const mpfr_t vr, const mpfr_t vi) {
    mpfr_prec_t p = mpfr_get_prec(vr);
    if (mpfr_get_prec(ur) > p) p = mpfr_get_prec(ur);
    mpfr_t dr, di;
    mpfr_inits2(p + 8, dr, di, (mpfr_ptr)0);
    mpfr_sub(dr, ur, vr, MPFR_RNDN);
    mpfr_sub(di, ui, vi, MPFR_RNDN);
    double ld = de_mpfr_log2_l1(dr, di);
    double lv = de_mpfr_log2_l1(vr, vi);
    mpfr_clears(dr, di, (mpfr_ptr)0);
    if (!isfinite(ld) && ld < 0) return -INFINITY;
    if (!isfinite(lv)) return 0.0;
    return ld - lv;
}

/* Endpoint cancellation guard.
 *
 * The interpreter evaluates at the precision of its operands and does not track
 * the digits a subtraction cancels.  tanh-sinh puts its last nodes at relative
 * distance 2^-L from an endpoint with L in the hundreds, and an integrand with a
 * removable singularity there loses about m L bits at such a node, m the order
 * of the cancellation: (E^x - 1 - x)/x^2 at x = 10^-40 and 30 digits evaluates
 * to -1e40, not 1/2, and NIntegrate returned -188.49 for 0.59962.  The node's
 * weight (~2^-L) does not save it, since the noise grows like 2^(mL).
 *
 * So a sample within 2^-20 of the span of an endpoint is checked by evaluating
 * it again with the abscissa widened by L + 32 bits (doubling the extra bits
 * until two evaluations agree to the target precision).  The bits the
 * working-precision value lost tell the endpoint's kind.  An endpoint whose
 * first deep sample lost nothing is marked clean and never checked again (a
 * genuine x^-1/2 singularity, or a regular integrand, pays one extra
 * evaluation).  Otherwise a second sample at least 16 bits deeper or shallower
 * decides: a loss that grows with the depth is a cancellation (lost bits per
 * bit of distance, the rate, then predicts the precision of every later deep
 * sample, evaluated directly at it plus a check); a loss that does not is the
 * integrand's own conditioning (E^E^E^x near x = 2 loses ~14 bits at any
 * depth), harmless within the guard bits and otherwise covered by a constant
 * extra precision. */
#define NI_GD_MIN_L     20.0
#define NI_GD_MAX_EXTRA 16384L
static bool ni_sample_mpfr_guarded(NiCtx* c, const mpfr_t xm, mpfr_t out_re, mpfr_t out_im,
                                   bool* handled) {
    *handled = false;
    if (!c->gd_on || !c->gd_a || !c->gd_b || c->spec.mode != NUMERIC_MODE_MPFR) return false;
    /* cheap screen in doubles: the bulk of the nodes is far from both ends
     * (and both ends clean is the common case) */
    if (c->gd_state[0] == 1 && c->gd_state[1] == 1) return false;
    {
        double xd = mpfr_get_d(xm, MPFR_RNDN);
        double ad = mpfr_get_d(c->gd_a, MPFR_RNDN), bd = mpfr_get_d(c->gd_b, MPFR_RNDN);
        double near = ldexp(fabs(bd - ad), -(int)NI_GD_MIN_L + 1);
        if (fabs(xd - ad) > near && fabs(bd - xd) > near) return false;
    }
    long p = c->spec.bits;
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(xm) + 8);
    mpfr_sub(t, xm, c->gd_a, MPFR_RNDN);
    double la = de_mpfr_log2_abs(t);
    mpfr_sub(t, c->gd_b, xm, MPFR_RNDN);
    double lb = de_mpfr_log2_abs(t);
    mpfr_sub(t, c->gd_b, c->gd_a, MPFR_RNDN);
    double ls = de_mpfr_log2_abs(t);
    mpfr_clear(t);
    int side = (la <= lb) ? 0 : 1;
    double L = ls - (side == 0 ? la : lb);
    if (!(L >= NI_GD_MIN_L) || !isfinite(L) || c->gd_state[side] == 1) return false;
    /* Once a loss has been measured, a sample whose loss is predicted to stay
     * well inside the guard bits (64) under the worse of the two readings --
     * growing with the depth at the measured rate, or constant -- needs no
     * second evaluation. */
    {
        double guard_bits = 64.0;
        int st = c->gd_state[side];
        double pred = st == 2 ? c->gd_rate[side] * L
                    : st == 3 ? (c->gd_lost1[side] * L / c->gd_L1[side] > c->gd_lost1[side]
                                 ? c->gd_lost1[side] * L / c->gd_L1[side] : c->gd_lost1[side])
                    : st == 4 ? c->gd_rate[side] - 32.0 : INFINITY;
        if (pred < guard_bits - 16.0) return false;
    }
    *handled = true;
    /* The probe evaluations' own precision is not the integrand's: a value that
     * cancelled to a 1-bit zero at the working precision is replaced by the
     * agreed one, so none of them may lower the delivered-precision record. */
    long sbits_saved = g_ni_sample_bits;

    long tb = p - 64 > 53 ? p - 64 : 53;              /* the target precision */
    double agree_l2 = -(double)tb - 8.0;
    mpfr_t pr, pi, cr, ci;
    mpfr_inits2((mpfr_prec_t)p, pr, pi, (mpfr_ptr)0);
    long extra;
    bool have_prev;
    if (c->gd_state[side] == 2 || c->gd_state[side] == 4) {
        extra = c->gd_state[side] == 2 ? (long)ceil(c->gd_rate[side] * L) + 32
                                       : (long)c->gd_rate[side];
        if (extra > NI_GD_MAX_EXTRA) extra = NI_GD_MAX_EXTRA;
        mpfr_set_prec(pr, (mpfr_prec_t)(p + extra));
        mpfr_set_prec(pi, (mpfr_prec_t)(p + extra));
        have_prev = ni_eval_mpfr_at_prec(c, xm, p + extra, pr, pi, NULL);
        extra *= 2;
    } else {
        /* the working-precision value through the ordinary (compiled) path */
        Expr* num = ni_eval_at(c, expr_new_mpfr_copy(xm));
        bool inexact;
        have_prev = num && get_approx_mpfr(num, pr, pi, &inexact)
                    && mpfr_number_p(pr) && mpfr_number_p(pi);
        expr_free(num);
        extra = (long)ceil(L) + 32;
    }
    mpfr_t v0r, v0i;                                  /* the first value, for the rate */
    mpfr_inits2(mpfr_get_prec(pr), v0r, v0i, (mpfr_ptr)0);
    if (have_prev) { mpfr_set(v0r, pr, MPFR_RNDN); mpfr_set(v0i, pi, MPFR_RNDN); }
    mpfr_inits2((mpfr_prec_t)p, cr, ci, (mpfr_ptr)0);
    bool ok = false;
    long db = 0, asked = 0;                           /* delivered / requested bits */
    bool capped = false;
    for (int it = 0; it < 8; it++) {
        mpfr_set_prec(cr, (mpfr_prec_t)(p + extra));
        mpfr_set_prec(ci, (mpfr_prec_t)(p + extra));
        long dbi = 0;
        bool got = ni_eval_mpfr_at_prec(c, xm, p + extra, cr, ci, &dbi);
        /* No number at a raised precision: the carried-exactly abscissa is
         * outside the integrand's domain (or past what it can evaluate), and
         * doubling the precision again -- up to NI_GD_MAX_EXTRA bits, eight
         * times -- resolves nothing; it only re-confirms the same non-number.
         * Decline the node now, as the loop would after its last doubling
         * (the rules truncate that side's tail at an unevaluable endpoint
         * node). */
        if (!got) break;
        db = dbi; asked = p + extra;
        /* held below the asked precision: raising it further cannot help */
        if (got && db > 0 && db < asked - 8) { capped = true; break; }
        if (got && have_prev && ni_mpfr_rel_log2(pr, pi, cr, ci) <= agree_l2) { ok = true; break; }
        if (got) {
            mpfr_set_prec(pr, mpfr_get_prec(cr)); mpfr_set_prec(pi, mpfr_get_prec(ci));
            mpfr_set(pr, cr, MPFR_RNDN); mpfr_set(pi, ci, MPFR_RNDN);
            have_prev = true;
        }
        if (extra >= NI_GD_MAX_EXTRA) { ok = got; break; }
        extra *= 2;
        if (extra > NI_GD_MAX_EXTRA) extra = NI_GD_MAX_EXTRA;
    }
    /* Capped: part of the integrand delivers fewer bits than asked whatever
     * the abscissa carries -- a machine number (its inexact arbitrary-precision
     * numbers are widened, see gd_exact) -- so raising the precision resolves
     * nothing: two such evaluations agree because they carry the same
     * rounding, and the agreed value went into the rule as if it were exact.
     * Its noise is measured instead by evaluating once more with the abscissa
     * rounded to xb = db - 12 bits: every operation then runs at xb bits, and
     * a cancellation's noise grows by the 2^12 the precision shrank, so
     * |f(db) - f(xb)| 2^-12 is the noise of the capped value -- while that
     * difference stays below the value; one of 1/16 of the value or more says
     * the sample has next to no digit left, and it is declined like an
     * unevaluable endpoint node (the rules truncate that side's tail there).
     * A kept sample's noise, times 2 L its distance (a tanh-sinh node 2^-L
     * from the endpoint weighs about 0.7 L h its distance), is the error it
     * leaves, recorded in g_ni_guard_noise_l2.  The abscissa must keep its
     * distance to the endpoint (L + 40 bits when that endpoint is not 0);
     * a node too deep for that is kept unmeasured, as before the guard. */
    if (capped) {
        ok = true;
        mpfr_srcptr ep = side == 0 ? c->gd_a : c->gd_b;
        double need = mpfr_zero_p(ep) ? 40.0 : L + 40.0;
        double v_l2 = de_mpfr_log2_l1(cr, ci);
        long xb = db - 12;
        if ((double)xb >= need && isfinite(v_l2)) {
            mpfr_t xl, qr, qi, d;
            mpfr_init2(xl, (mpfr_prec_t)xb);
            mpfr_set(xl, xm, MPFR_RNDN);
            mpfr_inits2((mpfr_prec_t)xb, qr, qi, (mpfr_ptr)0);
            if (ni_eval_mpfr_at_prec(c, xl, xb, qr, qi, NULL)) {
                mpfr_init2(d, (mpfr_prec_t)asked + 8);
                mpfr_sub(d, qr, cr, MPFR_RNDN);
                double n_l2 = de_mpfr_log2_abs(d);
                mpfr_sub(d, qi, ci, MPFR_RNDN);
                double ni_l2 = de_mpfr_log2_abs(d);
                mpfr_clear(d);
                if (ni_l2 > n_l2) n_l2 = ni_l2;
                if (isfinite(n_l2) && n_l2 > v_l2 - 4.0) {
                    ok = false;
                } else if (isfinite(n_l2) && n_l2 - 12.0 > v_l2 - (double)db + 16.0) {
                    double e_l2 = n_l2 - 12.0 + log2(L) + 1.0 + ls - L;
                    double g = g_ni_guard_noise_l2;
                    g_ni_guard_noise_l2 = !isfinite(g) ? e_l2
                        : (g > e_l2 ? g + log2(1.0 + exp2(e_l2 - g))
                                    : e_l2 + log2(1.0 + exp2(g - e_l2)));
                }
            }
            mpfr_clears(xl, qr, qi, (mpfr_ptr)0);
        }
        if (ok) {
            mpfr_set(out_re, cr, MPFR_RNDN);
            mpfr_set(out_im, ci, MPFR_RNDN);
        }
    } else if (ok) {
        /* bits the first evaluation lost, against the agreed value */
        double lost = 0.0;
        if (mpfr_number_p(v0r) && mpfr_number_p(v0i)) {
            double rl = ni_mpfr_rel_log2(v0r, v0i, cr, ci);
            double v0p = (double)mpfr_get_prec(v0r);
            lost = isfinite(rl) ? v0p + rl : 0.0;
            if (lost < 0.0) lost = 0.0;
        }
        double guard = (double)(p - tb);
        if (c->gd_state[side] == 0) {
            if (lost < 8.0) c->gd_state[side] = 1;
            else { c->gd_state[side] = 3; c->gd_L1[side] = L; c->gd_lost1[side] = lost; }
        } else if (c->gd_state[side] == 3 && fabs(L - c->gd_L1[side]) >= 16.0) {
            double growth = (lost - c->gd_lost1[side]) / (L - c->gd_L1[side]);
            double worst = lost > c->gd_lost1[side] ? lost : c->gd_lost1[side];
            if (growth < 0.25) {
                if (worst < guard - 16.0) c->gd_state[side] = 1;
                else { c->gd_state[side] = 4; c->gd_rate[side] = worst + 32.0; }
            } else {
                double r1 = c->gd_lost1[side] / c->gd_L1[side], r2 = lost / L;
                c->gd_state[side] = 2;
                c->gd_rate[side] = r1 > r2 ? r1 : r2;
            }
        }
        mpfr_set(out_re, cr, MPFR_RNDN);
        mpfr_set(out_im, ci, MPFR_RNDN);
    }
    mpfr_clears(pr, pi, cr, ci, v0r, v0i, (mpfr_ptr)0);
    g_ni_sample_bits = sbits_saved;
    return ok;
}

static bool ni_sample_mpfr(void* vctx, const mpfr_t x, mpfr_t out_re, mpfr_t out_im) {
    NiCtx* c = (NiCtx*)vctx;
    mpfr_t xm;
    mpfr_init2(xm, mpfr_get_prec(x));
    mpfr_mul_d(xm, x, creal(c->x_scale), MPFR_RNDN);
    mpfr_add_d(xm, xm, creal(c->x_shift), MPFR_RNDN);
    if (c->gd_on) {
        bool handled;
        bool ok = ni_sample_mpfr_guarded(c, xm, out_re, out_im, &handled);
        if (handled) { mpfr_clear(xm); return ok; }
    }
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

/* `out` (consumed) with its precision lowered to what an absolute error of
 * 2^noise_l2 leaves, when that is below target_bits: a real or complex MPFR
 * number is rounded to those bits and NIntegrate::precw reports it.  Any
 * other result is returned as it is. */
static Expr* ni_cap_precision(Expr* out, double noise_l2, long target_bits,
                               bool endpoint) {
    Expr *re = NULL, *im = NULL;
    bool cpx = is_complex(out, &re, &im);
    if (!cpx) re = out;
    if (re->type != EXPR_MPFR || (cpx && im->type != EXPR_MPFR)) return out;
    double mag = de_mpfr_log2_l1(re->data.mpfr, cpx ? im->data.mpfr : re->data.mpfr);
    if (!cpx) mag = de_mpfr_log2_abs(re->data.mpfr);
    if (!isfinite(mag)) return out;
    double acc = mag - noise_l2;
    if (!(acc < (double)target_bits)) return out;
    long nb = acc > 2.0 ? (long)floor(acc) : 2;
    if (endpoint)
        ni_warn("precw", "an endpoint of the integration range could not be "
                "evaluated to the WorkingPrecision; the result has about %.3g "
                "digits.", numeric_bits_to_digits(nb));
    else
        ni_warn("precw", "an endpoint cancellation in the integrand could not be "
                "resolved past the precision of its inexact constants; the result "
                "has about %.3g digits.", numeric_bits_to_digits(nb));
    mpfr_t r, i;
    mpfr_init2(r, (mpfr_prec_t)nb);
    mpfr_set(r, re->data.mpfr, MPFR_RNDN);
    Expr* nout;
    if (cpx) {
        mpfr_init2(i, (mpfr_prec_t)nb);
        mpfr_set(i, im->data.mpfr, MPFR_RNDN);
        nout = make_complex(expr_new_mpfr_move(r), expr_new_mpfr_move(i));
    } else {
        nout = expr_new_mpfr_move(r);
    }
    expr_free(out);
    return nout;
}

/* Build a result from a guarded-precision (re,im) pair: round to target_bits,
 * chop a noise-level imaginary part, return real or Complex. */
static Expr* ni_mpfr_build(const mpfr_t re, const mpfr_t im, long target_bits) {
    mpfr_t r, i;
    mpfr_init2(r, (mpfr_prec_t)target_bits);
    mpfr_init2(i, (mpfr_prec_t)target_bits);
    mpfr_set(r, re, MPFR_RNDN);
    mpfr_set(i, im, MPFR_RNDN);
    /* Noise-level test in log2 (1e-12 = 2^-39.86), so a result beyond the
     * double range keeps a genuine imaginary part (Inf <= 1e-12 Inf was true). */
    double ar = de_mpfr_log2_abs(r);
    double ai = de_mpfr_log2_abs(i);
    const double L_CHOP = -39.863137138648348;   /* log2(1e-12) */
    Expr* out;
    if (mpfr_zero_p(i) || ai <= L_CHOP + ar) {
        out = expr_new_mpfr_copy(r);
    } else if (ar <= L_CHOP + ai) {
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
    double acc_goal;        /* AccuracyGoal digits (nc_accuracy.h sentinels)  */
    bool   acc_goal_user;   /* AccuracyGoal given explicitly (else the default
                             * floor is rescaled to the integrand, see
                             * ni_scale_accuracy_floor)                       */
    double prec_goal;       /* PrecisionGoal digits (nc_accuracy.h sentinels) */
    int    max_recursion;   /* -1 => Automatic                                */
    int    min_recursion;   /* default 0                                      */
    long   max_points;      /* -1 => Automatic                                */
    Expr*  exclusions;      /* borrowed; NULL if none                         */
    /* Fixed-rule sub-options (RiemannRule / TrapezoidalRule / NewtonCotesRule). */
    int    rule_type;       /* Riemann sampling: 0 Left, 1 Right, 2 Midpoint   */
    bool   romberg;         /* Romberg (Richardson) extrapolation; default on  */
    int    nc_points;       /* Newton–Cotes points per panel (2..5); default 3 */
    /* The outer spec's nodes as written ({x, n0, n1, ...} args, borrowed), so
     * the MPFR paths take endpoints such as 2 Pi at the working precision
     * rather than as their double (see ni_endpoint_mpfr). */
    Expr* const* spec_nodes;
    size_t       spec_nnodes;
    /* At WorkingPrecision ni_dispatch_1d hands the rules pre_body, the
     * integrand with its exact constants numericalized once at `bits`
     * (ni_prenumericalize); exact_body is the integrand as given (both
     * borrowed, NULL when not in effect).  The endpoint cancellation guard
     * re-evaluates the deep nodes of pre_body through exact_body, since a
     * constant carried at the working precision caps the precision of every
     * sum it enters at that precision. */
    const Expr*  pre_body;
    Expr*        exact_body;
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

/* Delegates to the shared parser (nc_accuracy.h): Automatic -> NC_GOAL_AUTO
 * (-1), Infinity -> NC_GOAL_OFF (+Inf), MachinePrecision -> ~15.95 digits,
 * positive number -> itself, else false.  MachinePrecision is the newly-accepted
 * setting the old parser rejected.  The Infinity sentinel changed from -1 to
 * +Inf, but every consumer of acc_goal/prec_goal below still disables the
 * criterion for it: each computes its term as `(goal > 0) ? 10^-goal : 0` (and
 * pow(10, -Inf) == 0) or resolves `goal < 0` (Automatic) to a default, so +Inf
 * drops the term exactly as the old -1 did for Infinity. */
static bool ni_parse_goal(Expr* v, double* out) {
    return nc_parse_goal(v, out);
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

static bool ni_apply_option_value(const char* name, Expr* rhs, NiOpts* o);

static bool ni_apply_option(Expr* rule, NiOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    /* NIntegrate is HoldAll, so option values arrive unevaluated: the numeric
     * settings are evaluated first (PrecisionGoal -> i with i = 1 stayed
     * unevaluated).  Method keeps its literal form (the parsed options borrow
     * its strings), Exclusions is read structurally, and the monitors are
     * expressions in the integration variable. */
    if (name == SYM_WorkingPrecision || name == SYM_AccuracyGoal
        || name == SYM_PrecisionGoal || name == SYM_MaxRecursion
        || name == SYM_MinRecursion || name == SYM_MaxPoints) {
        Expr* v = eval_and_free(expr_copy(rhs));
        bool ok = v && ni_apply_option_value(name, v, o);
        expr_free(v);
        return ok;
    }
    return ni_apply_option_value(name, rhs, o);
}

static bool ni_apply_option_value(const char* name, Expr* rhs, NiOpts* o) {

    if (name == SYM_Method) {
        Expr* m = rhs;
        /* Method -> {"name", subopts...} or Method -> "name" or symbol. */
        bool is_list = m->type == EXPR_FUNCTION
            && m->data.function.head->type == EXPR_SYMBOL
            && m->data.function.head->data.symbol.name == SYM_List
            && m->data.function.arg_count >= 1;
        /* "Trapezoidal" is the plain nested trapezoidal rule, doubled until
         * two levels agree: the method for a periodic integrand over a period,
         * where it converges exponentially and a Romberg table (the
         * "TrapezoidalRule" default) would only mix in the aliased coarse
         * levels.  An explicit "RombergQuadrature" sub-option still wins. */
        Expr* mn = is_list ? m->data.function.args[0] : m;
        if (mn->type == EXPR_STRING && !strcmp(mn->data.string, "Trapezoidal"))
            o->romberg = false;
        if (is_list) {
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
    if (name == SYM_AccuracyGoal) {
        o->acc_goal_user = true;
        return ni_parse_goal(rhs, &o->acc_goal);
    }
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

/* Working precision of the run in decimal digits (machine, or the requested
 * MPFR target). */
static double ni_wp_digits(const NiOpts* o) {
    if (o->prec_mpfr) {
        long bits = o->target_bits > 0 ? o->target_bits : o->bits;
        return numeric_bits_to_digits(bits);
    }
    return NUMERIC_MACHINE_PRECISION_DIGITS;
}

#ifdef USE_MPFR
static double ni_mpfr_reltol(const NiOpts* o);
#endif

/* The goals the rules were actually asked to meet, in digits, for the
 * shortfall test of ni_warn_goal.  An Automatic PrecisionGoal resolves exactly
 * as the rules resolve it -- half the machine digits (ni_machine_tols, capped at
 * 14) at machine precision, WorkingPrecision - 10 under MPFR (ni_mpfr_reltol)
 * -- and an Automatic AccuracyGoal adds no absolute term, as in ni_machine_tols.
 * The shared nc_resolve_goal reads Automatic as WorkingPrecision - 2 (~14
 * digits at machine precision): testing the 8-digit Automatic run against that
 * both warned on results that had met their goal and stated a tolerance
 * (~1e-14 |I|) the run had never been asked for. */
static void ni_effective_goals(const NiOpts* o, double* acc_digits,
                               double* prec_digits) {
    double rt;
#ifdef USE_MPFR
    if (o->prec_mpfr) rt = ni_mpfr_reltol(o);
    else
#endif
    {
        double at;
        ni_machine_tols(o, &rt, &at);
    }
    *prec_digits = (o->prec_goal > 0.0) ? o->prec_goal
                 : (rt > 0.0 ? -log10(rt) : INFINITY);
    *acc_digits = (o->acc_goal > 0.0) ? o->acc_goal : INFINITY;
}

/* Standard AccuracyGoal/PrecisionGoal shortfall warning (nc_accuracy.h).  Warns
 * only when the achieved error estimate actually exceeds the combined tolerance
 * 10^-a + |result| 10^-p, and stays silent under an enclosing numeric probe's
 * warning mute (as ni_warn does), so an inner iterated integral's convergence
 * noise never reaches the user.  The caller keeps returning its best estimate. */
static void ni_warn_goal(const NiOpts* o, double achieved_err, double result_mag) {
    if (arith_warnings_muted()) return;
    /* a machine pass that left the double range is re-run in MPFR, which
     * gives its own verdict; an MPFR pass (any sample of a machine-precision
     * geometry probe may have set the flag) is the final one and warns */
    if (!o->prec_mpfr && g_ni_beyond && !g_ni_escalated) return;
    double a_dig, p_dig;
    ni_effective_goals(o, &a_dig, &p_dig);
    double tol = nc_combined_tol(a_dig, p_dig, result_mag, ni_wp_digits(o));
    if (achieved_err > tol) nc_warn_goal("NIntegrate", achieved_err, tol);
}

#ifdef USE_MPFR
/* The goal warning for an MPFR run that did not converge.  A result (or level
 * difference) beyond the double range cannot be compared through doubles, so
 * the shortfall is stated without the numbers; the driver's own convergence
 * test, which is scale-free, has already said the goal was missed. */
static void ni_warn_goal_mpfr(const NiOpts* o, bool in_range,
                              double achieved_err, double result_mag) {
    if (in_range) { ni_warn_goal(o, achieved_err, result_mag); return; }
    if (arith_warnings_muted()) return;
    mth_message("NIntegrate", "accgl",
                "requested accuracy goal not met; returning the best "
                "approximation obtained.");
}
#endif

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
/* A probe sample (the endpoint probes: oscillation, power law below the
 * double range).  The probes sample ever nearer an endpoint, where the
 * interpreter's fallback can cost without bound: Cos[Exp[1/x^2]] at 2^-12 has
 * a compiled value of NaN (e^(1.7e7) overflows), and the interpreter, which
 * carries the overflow on as an arbitrary-precision number, then needs a
 * 2^24-bit argument reduction for the cosine -- 15 s for one sample, and
 * minutes for the next.  So when the integrand compiled and its compiled value
 * is not finite because an intermediate overflowed (FE_OVERFLOW; a value that
 * merely goes complex -- the real program's Sqrt of a negative -- raises only
 * FE_INVALID and is still asked of the interpreter), the probe is told the
 * sample left the double range (-1) without evaluating it.  1: *v written;
 * 0: no numeric value there. */
static int ni_probe_machine(NiCtx* c, double x, double _Complex* v) {
    if (c->spec.mode == NUMERIC_MODE_MACHINE && c->map_mode == NI_MAP_AFFINE
        && creal(c->x_scale) == 1.0 && cimag(c->x_scale) == 0.0 && c->x_shift == 0.0) {
        ni_ac_prepare(c);
        if (c->ac) {
            feclearexcept(FE_OVERFLOW);
            double _Complex z;
            bool ok = autocompiled_eval_complex(c->ac, &x, &z)
                      && isfinite(creal(z)) && isfinite(cimag(z));
            bool over = fetestexcept(FE_OVERFLOW) != 0;
            if (ok) { *v = z; return 1; }
            if (over) return -1;
        }
    }
    return ni_sample_machine(c, x, v) ? 1 : 0;
}

typedef struct { bool have, conv; double _Complex val; double err; } NiAtt;

static void ni_consider(NiAtt* best, NiAtt a) {
    if (!a.have) return;
    if (!best->have
        || (a.conv && !best->conv)
        || (a.conv == best->conv && a.err < best->err))
        *best = a;
}

#define NI_GK_SING_DEPTH 4   /* halvings at an endpoint before tanh-sinh */

static NiAtt ni_try_gk(NiCtx* ctx, double a, double b, const NiOpts* o,
                       double reltol, double abstol) {
    int max_subdiv; long max_eval;
    ni_gk_budget(o, &max_subdiv, &max_eval);
    /* Endpoint singularity handler after NI_GK_SING_DEPTH halvings at an end
     * (gkadapt.h): an integrable x^alpha / Log endpoint is then settled by one
     * tanh-sinh panel, where bisecting it to the goal left ~3 digits beyond
     * PrecisionGoal (11 of 15 on (1-x^4)^(1/4)). */
    GkResult R = gk_integrate_machine_ex(ni_sample_machine, ctx, a, b,
                                         abstol, reltol, max_subdiv, max_eval, true,
                                         NI_GK_SING_DEPTH);
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
    int vtok = mth_bind_push_name(var, 0);   /* unwind-safe detach */
    eval_clock_bump();
    Expr* gp = eval_and_free(dcall);
    mth_bind_pop(vtok);
    eval_clock_bump();
    if (!gp) return NULL;
    if (levin_occurs(gp, SYM_D) || levin_occurs(gp, "Derivative")) {
        expr_free(gp); return NULL;
    }
    return gp;
}

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
    /* expr_copy only shares the node, and eval_flatten_args rewrites the args
     * of whatever node it is given whatever its head: on Exp[-x^2] it turned
     * the shared integrand into Power[E, -1, x^2] = E^((-1)^(x^2)), so every
     * WorkingPrecision half-line integral sampled after this probe evaluated
     * a different function (Exp[-x] on [0, oo) gave 4.6e619).  Flatten only a
     * Times, and only a private top-level node. */
    Expr* fb = expr_copy((Expr*)body);
    if (levin_is_head(fb, SYM_Times)) {
        fb = expr_unshare(fb);
        eval_flatten_args(fb, SYM_Times);
    }
    bool ok = levin_detect_kernel(fb, var, &k, &g, &amp);
    expr_free(fb);
    if (!ok) return false;
    /* The rule's nodes sit on the zeros of the one kernel, which presumes an
     * amplitude that does not oscillate itself: BesselJ[0, x] Cos[x/2] (two
     * frequencies) ground through every level, 26 s at 30 digits, before
     * failing.  Such an integrand is a sum of phase groups for the contour
     * rotation instead. */
    bool amp_osc = amp && nc_has_phase(amp, var);
    expr_free(amp);
    if (amp_osc) { expr_free(g); return false; }
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
    NiCtx camp = *ctx; ni_ctx_rebody(&camp, amp);
    NiCtx cg   = *ctx; ni_ctx_rebody(&cg,   g);
    NiCtx cgp  = *ctx; ni_ctx_rebody(&cgp,  gp);

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
                                          const NiOpts* o, double min_phase,
                                          bool* weak, mpfr_t re, mpfr_t im,
                                          bool* conv) {
    if (ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return false;
    const char* var = ctx->bind->name;
    LevinKernel kind; Expr* g = NULL; Expr* amp = NULL;
    if (!levin_detect_kernel(ctx->body, var, &kind, &g, &amp)) return false;
    Expr* gp = levin_phase_derivative(g, var);
    if (!gp) { expr_free(g); expr_free(amp); return false; }

    NiCtx camp = *ctx; ni_ctx_rebody(&camp, amp);
    NiCtx cg   = *ctx; ni_ctx_rebody(&cg,   g);
    NiCtx cgp  = *ctx; ni_ctx_rebody(&cgp,  gp);
    if (min_phase > 0.0) {
        /* Levin pays a dense complex MPFR solve per order (n up to 64): it is
         * the rule for MANY oscillations, where a quadrature rule needs nodes
         * in proportion to them.  Over a few periods -- the first lobe of
         * Sin[x^3], phase 0 .. Pi -- the endpoint-robust rules are far cheaper
         * (and Levin, with g' = 0 at an end, needed its top order), so the
         * caller tries them first and keeps Levin as the fallback.  The phase
         * variation is read at machine precision on a 16-panel grid. */
        NiCtx gm = cg;
        gm.spec = numeric_machine_spec();
        double tv = 0.0, prev = 0.0;
        for (int i = 0; i <= 16 && isfinite(tv); i++) {
            double _Complex v;
            if (!ni_sample_machine(&gm, a + (b - a) * i / 16.0, &v)) tv = INFINITY;
            else {
                if (i > 0) tv += fabs(creal(v) - prev);
                prev = creal(v);
            }
        }
        if (tv < min_phase) {
            *weak = true;
            expr_free(g); expr_free(amp); expr_free(gp); return false;
        }
    }
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
/* Level budget of the equally-spaced rules.  ncrule starts from 4 panels and
 * compares it with 8; MaxRecursion -> r allows 2^(r+3) panels, i.e. r + 1
 * doublings past the first rule (8193 nodes at r = 10, enough to integrate a
 * trigonometric polynomial of degree up to 8192 exactly by the trapezoidal
 * rule).  Automatic keeps the former caps of 2^24 (Riemann) and 2^20 panels. */
static void ni_ncr_budget(const NiOpts* o, int* max_levels, long* max_eval) {
    int r = o->max_recursion;
    if (o->method == NI_RIEMANN) {
        *max_levels = (r < 0) ? 22 : (r > 28 ? 29 : (r < 5 ? 6 : r + 1));
        *max_eval   = (o->max_points > 0) ? o->max_points : 4000000;
    } else {
        *max_levels = (r < 0) ? 18 : (r > 28 ? 29 : (r < 3 ? 4 : r + 1));
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
/* An oscillatory endpoint singularity -- a phase that grows without bound at
 * the endpoint (Cos[Log[x]/x], Sin[1/(x (1 - x))]) -- detected numerically,
 * so no particular form of the integrand is needed: the real or imaginary
 * part changes sign many times on a log-uniform grid of distances in two
 * shells next to the endpoint, 10^-11..10^-9 and 10^-6..10^-4 of the span.
 * An oscillation of fixed frequency, however fast (Sin[10^6 x]), is smooth
 * on the deep shell; only an unbounded phase still changes sign there.  The
 * shallow shell rejects rounding noise: a cancelling expression
 * ((E^x - 1 - x)/x^2 at machine precision) has random signs at 10^-10 but
 * none at 10^-5.  `probe` reports the signs of Re and Im of f at distance
 * d span from the endpoint; a shell stops early once it can no longer reach
 * the threshold, so a smooth integrand costs ~23 samples an endpoint. */
#define NI_OSE_POINTS   32
#define NI_OSE_MIN_SC   10
/* An unbounded phase sampled on a log grid gives effectively RANDOM signs
 * (successive phases differ by thousands of radians), so the count of sign
 * changes over the 31 transitions of a shell is ~Binomial(31, 1/2): below
 * NI_OSE_MIN_SC with probability ~1.5% per shell, ~3% per endpoint.  That
 * miss is a silent wrong answer, not a slow one: the endpoint falls through to
 * Gauss-Kronrod, which settles on an under-sampled value (x^a Sin[a/(2 - x)]
 * on [0, 2] was wrong by up to 10% at 3.7% of the parameters a in [0, 5],
 * with no message, and NMaximize over a picked exactly those spurious
 * peaks).  A shell that ends with at least NI_OSE_AMBIG changes -- a smooth
 * integrand has 0-2 -- is therefore re-sampled on a grid NI_OSE_REFINE times
 * finer and accepted at a quarter of its transitions: ~124 transitions with
 * mean 62 against a threshold of 31, a miss far below 1e-8.  Only the
 * ambiguous ~3% of genuine oscillations pay for it. */
#define NI_OSE_AMBIG    3
#define NI_OSE_REFINE   4

/* 1: signs written; 0: no numeric sample there; -1: stop probing. */
typedef int (*NiOseSign)(void* ctx, double d, int* sr, int* si);

/* Sample cost that keeps growing toward the endpoint stops the probe: an
 * integrand whose evaluation cost grows without bound there (Cos[Exp[1/x]],
 * whose huge phase is still a valid double cosine, computed through an MPFR
 * argument reduction of 10^4 .. 10^11 bits) is not one the lobe rule can
 * follow anyway, and the deep shell would cost hours.  The test is the TREND
 * -- a sample over NI_OSE_COST_MIN seconds costing NI_OSE_COST_RATIO times
 * the one NI_OSE_COST_LAG samples further out -- not the level: an
 * MPFR-valued parameter (a 60-digit a in x^a Sin[a/(2 - x)]) makes every
 * sample cost milliseconds, flatly, and must still be probed. */
#define NI_OSE_COST_RATIO 6.0
#define NI_OSE_COST_MIN   2e-3
#define NI_OSE_COST_LAG   8

/* Set by ni_osc_scan when the probe gave up (a sample left the double range,
 * or the cost trend tripped) rather than reaching a verdict: the endpoint may
 * then still oscillate, beyond what can be sampled there (ni_try_phase_coord).
 * Read by the caller right after ni_osc_endpoint_machine. */
static bool s_ni_ose_gaveup = false;

/* One log-uniform scan of npts distances over [10^lo, 10^hi] (outer edge
 * first), stopping once `need` sign changes are seen (return 1) or can no
 * longer be reached (0); -1 when the probe gave up or the cost trend tripped.
 * *best receives the larger of the Re / Im sign-change counts. */
static int ni_osc_scan(NiOseSign probe, void* pc, double lo, double hi,
                       int npts, int need, int* best, double* growth) {
    int scr = 0, sci = 0, pr = 0, pi = 0;
    double cost[NI_OSE_POINTS * NI_OSE_REFINE];
    double c0 = 0.0;
    int lag = NI_OSE_COST_LAG * (npts / NI_OSE_POINTS);
    *best = 0;
    for (int i = 0; i < npts; i++) {
        int left = npts - i;
        if (*best + left < need) return 0;
        double e = hi - (hi - lo) * ((double)i + 0.5) / (double)npts;
        int sr, si;
        double t0 = mth_mono_seconds();
        int st = probe(pc, pow(10.0, e), &sr, &si);
        double dt = mth_mono_seconds() - t0;
        cost[i] = dt;
        if (i == 0) c0 = dt;
        if (i >= lag && dt > NI_OSE_COST_MIN
            && dt > NI_OSE_COST_RATIO * cost[i - lag]) { s_ni_ose_gaveup = true; return -1; }
        if (st < 0) { s_ni_ose_gaveup = true; return -1; }   /* the probe gave up */
        if (st == 0) continue;
        if (sr != 0) { if (pr != 0 && sr != pr) scr++; pr = sr; }
        if (si != 0) { if (pi != 0 && si != pi) sci++; pi = si; }
        *best = scr > sci ? scr : sci;
        if (*best >= need) {
            /* how much dearer the samples got across the shell */
            if (dt > NI_OSE_COST_MIN && c0 > 0.0) *growth = dt / c0;
            return 1;
        }
    }
    return 0;
}

static bool ni_osc_shell(NiOseSign probe, void* pc, double lo, double hi,
                         double* growth) {
    *growth = 1.0;
    int best = 0;
    /* the coarse scan, run to its end only while it can still be ambiguous */
    int st = ni_osc_scan(probe, pc, lo, hi, NI_OSE_POINTS, NI_OSE_MIN_SC,
                         &best, growth);
    if (st != 0) return st > 0;
    if (best < NI_OSE_AMBIG) return false;
    int fine = NI_OSE_POINTS * NI_OSE_REFINE;
    return ni_osc_scan(probe, pc, lo, hi, fine, (fine - 1) / 4, &best, growth) > 0;
}

static bool ni_osc_endpoint_probe(NiOseSign probe, void* pc) {
    /* the deep shell is 10^5 times nearer: not entered when the samples of
     * the shallow one already grew dearer toward the endpoint */
    double g1, g2;
    s_ni_ose_gaveup = false;
    if (!ni_osc_shell(probe, pc, -6.0, -4.0, &g1)) return false;
    if (g1 > NI_OSE_COST_RATIO) return false;
    return ni_osc_shell(probe, pc, -11.0, -9.0, &g2);
}

typedef struct { NiCtx* c; double a, b; int side; } NiOseMachProbe;

/* A sample that leaves the double range stops the probe: the samples are
 * taken from the outside in, and one nearer the endpoint can cost without
 * bound (Cos[Exp[1/x]] at x = 10^-11 is the cosine of a number with 10^11
 * integer bits).  The flag stays set, so a machine pass is still re-run in
 * MPFR as any other out-of-range sample asks. */
static int ni_ose_machine_sign(void* vp, double d, int* sr, int* si) {
    NiOseMachProbe* q = (NiOseMachProbe*)vp;
    double span = q->b - q->a;
    double x = (q->side > 0) ? q->a + span * d : q->b - span * d;
    if (x == q->a || x == q->b) return 0;
    bool saved = g_ni_beyond;
    if (saved && !g_ni_escalated) return -1;   /* the pass is already void */
    g_ni_beyond = false;
    double _Complex v;
    int st = ni_probe_machine(q->c, x, &v);
    bool beyond = g_ni_beyond;
    g_ni_beyond = saved || beyond;
    if (beyond || st < 0) return -1;
    if (st == 0) return 0;
    *sr = (creal(v) > 0) - (creal(v) < 0);
    *si = (cimag(v) > 0) - (cimag(v) < 0);
    return 1;
}

/* The machine-precision test of one endpoint (side +1: a, -1: b). */
static bool ni_osc_endpoint_machine(NiCtx* c, double a, double b, int side) {
    NiOseMachProbe q = { c, a, b, side };
    return ni_osc_endpoint_probe(ni_ose_machine_sign, &q);
}

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
/* One exponential-endpoint piece: x = end_a +- span e^-z on z in [0, oo),
 * exp-sinh first, then between the zeros + Wynn. */
static NiAtt ni_endpoint_piece(NiCtx* ctx, NiMapMode mode, double end, double span,
                               const NiOpts* o, double reltol, bool osc) {
    ctx->span = span; ctx->map_mode = mode; ctx->end_a = end;
    double _Complex v = NAN; double err = INFINITY;
    /* exp-sinh only when the endpoint is not known to oscillate: on the
     * non-decaying oscillation the map makes of an oscillatory singularity
     * it spends every level and declines (two thirds of the cost of an
     * x^a Sin[a/(2 - x)] inside an optimiser) */
    bool conv = !osc && dequad_halfline_machine(ni_sample_machine, ctx, 0.0, reltol,
                                                ni_de_levels(o), &v, &err);
    if (!conv) {
        double _Complex ov; double oerr;
        /* the lobe extrapolation is asked for three digits beyond the goal
         * (to at most 1e-13): it converges geometrically, so this costs a few
         * lobes, and it met an Automatic 1e-8 with 8 digits where the rules
         * of the other paths overshoot to 11-12 */
        double ort = fmax(reltol * 1e-3, 1e-13);
        if (ort > reltol) ort = reltol;
        if (osc_integrate_machine(ni_sample_machine, ctx, 0.0, 0.0, true,
                                  ort, 600, &ov, &oerr)) {
            v = ov; err = oerr; conv = true;
        } else if (isfinite(oerr) && isfinite(creal(ov)) && isfinite(cimag(ov))) {
            v = ov; err = oerr;     /* the extrapolation, unconverged */
        }
    }
    ctx->map_mode = NI_MAP_AFFINE;
    NiAtt at;
    at.have = isfinite(creal(v)) && isfinite(cimag(v));
    at.conv = conv; at.val = v; at.err = err;
    return at;
}

static NiAtt ni_try_endpoint_sing(NiCtx* ctx, double a, double b,
                                  const NiOpts* o, bool force, int osc_lo, int osc_hi) {
    NiAtt best = { false, false, 0.0, INFINITY };
    /* Only meaningful on the plain real axis (not an affine contour segment). */
    if (!(b > a) || ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return best;

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    /* Transform only an endpoint detected singular -- by magnitude, or as an
     * oscillatory singularity of bounded amplitude (Cos[Log[x]/x]): the
     * exponential map damps a singularity at the z -> oo end only, so applying
     * it to a regular endpoint (the wrong mirror) would instead sample the
     * *other*, singular end at z = 0 and converge to the wrong value.  Both
     * ends singular (Sin[1/(x (1 - x))]): the range is split at its midpoint
     * and each half mapped at its own end.  With an explicit method (`force`)
     * and no endpoint detected singular, fall back to trying both mirrors. */
    /* osc_lo / osc_hi: the caller's oscillation verdict (1/0), or -1 */
    bool olo = osc_lo >= 0 ? osc_lo > 0 : ni_osc_endpoint_machine(ctx, a, b, +1);
    bool ohi = osc_hi >= 0 ? osc_hi > 0 : ni_osc_endpoint_machine(ctx, a, b, -1);
    bool lo = olo || ni_endpoint_singular(ctx, a, b, +1);
    bool hi = ohi || ni_endpoint_singular(ctx, a, b, -1);
    if (lo && hi) {
        double m = a + 0.5 * (b - a);
        NiAtt p1 = ni_endpoint_piece(ctx, NI_MAP_EXP_LO, a, m - a, o, reltol, olo);
        NiAtt p2 = ni_endpoint_piece(ctx, NI_MAP_EXP_HI, b, b - m, o, reltol, ohi);
        best.have = p1.have && p2.have;
        best.conv = p1.conv && p2.conv;
        best.val = p1.val + p2.val;
        best.err = p1.err + p2.err;
        return best;
    }
    if (lo) return ni_endpoint_piece(ctx, NI_MAP_EXP_LO, a, b - a, o, reltol, olo);
    if (hi) return ni_endpoint_piece(ctx, NI_MAP_EXP_HI, b, b - a, o, reltol, ohi);
    if (force) {
        ni_consider(&best, ni_endpoint_piece(ctx, NI_MAP_EXP_LO, a, b - a, o, reltol, false));
        ni_consider(&best, ni_endpoint_piece(ctx, NI_MAP_EXP_HI, b, b - a, o, reltol, false));
    }
    return best;
}

/* ------------------------------------------------------------------ *
 *  Phase coordinate at an oscillatory endpoint                        *
 * ------------------------------------------------------------------ */

/* An integrand amp(x) K(g(x)) -- K one of cos, sin, e^{i .} -- whose phase g
 * is monotone near an endpoint and unbounded there is, in the phase itself as
 * the coordinate, a Fourier integral on a half line:
 *
 *     Int_end^{x0} amp K(g) dx  =  Int_{v0}^{oo} amp(x(v)) K(s v) / |g'(x(v))| dv,
 *
 * s = +-1 the sign that sends v = s g to +oo at the endpoint and x(v) the
 * inverse of the phase.  The zeros of K are then equally spaced, and the
 * between-the-zeros + Wynn rule (oscint.c) needs a few dozen lobes.  It is the
 * rule for a phase that grows faster than any power of 1/distance --
 * Cos[Exp[1/x]] on [0, 1], whose phase overflows a double inside x < 1/709:
 * the exponential endpoint map turns that phase into e^(e^z), whose lobes the
 * zero finder cannot follow, and the integrand cannot even be sampled near the
 * endpoint (the oscillation probe gives up there).  In the phase coordinate the
 * same integral is cos(v)/(v Log[v]^2) on [e^2, oo): every sample stays in the
 * double range.  The inverse x(v) is computed numerically, by Newton's method
 * on the distance to the endpoint safeguarded by a bracket, so no particular
 * form of the phase is needed beyond a closed-form derivative. */
typedef struct {
    NiCtx cg, cgp;          /* the phase and its derivative (machine)       */
    NiCtx camp;             /* the amplitude                                 */
    LevinKernel kind;
    double end;             /* the oscillatory endpoint                      */
    int side;               /* +1: the lower endpoint, -1: the upper one     */
    double s;               /* v = s g -> +oo at the endpoint                */
    double d0;              /* distance of the split point from `end`        */
    double dlast;           /* the previous inverse (warm start)             */
    Expr *g, *gp, *amp;     /* owned bodies of cg, cgp, camp                 */
} NiPhaseMap;

/* v(d) = s g(end + side d); false when not a finite real.  *over is set
 * when the phase left the double range (it is then larger than any finite v:
 * the endpoint side of the bracket). */
static bool ni_pm_v(NiPhaseMap* pm, double d, double* v, bool* over) {
    *over = false;
    double x = pm->end + pm->side * d;
    bool saved = g_ni_beyond;
    g_ni_beyond = false;
    double _Complex z;
    bool ok = ni_sample_machine(&pm->cg, x, &z);
    bool beyond = g_ni_beyond;
    g_ni_beyond = saved;
    if (!ok) { *over = beyond; return false; }
    if (fabs(cimag(z)) > 1e-12 * (1.0 + fabs(creal(z)))) return false;
    *v = pm->s * creal(z);
    return isfinite(*v);
}

/* dv/dd at distance d (finite real), else false (*over as in ni_pm_v). */
static bool ni_pm_dv2(NiPhaseMap* pm, double d, double* dv, bool* over) {
    double _Complex z;
    bool saved = g_ni_beyond;
    g_ni_beyond = false;
    bool ok = ni_sample_machine(&pm->cgp, pm->end + pm->side * d, &z);
    *over = g_ni_beyond;
    g_ni_beyond = saved;
    if (!ok) return false;
    if (fabs(cimag(z)) > 1e-12 * (1.0 + fabs(creal(z)))) return false;
    *dv = pm->s * pm->side * creal(z);
    return isfinite(*dv) && *dv != 0.0;
}
static bool ni_pm_dv(NiPhaseMap* pm, double d, double* dv) {
    bool over;
    return ni_pm_dv2(pm, d, dv, &over);
}

/* The distance d in (0, d0] with v(d) = v (v >= v(d0)): Newton on d inside a
 * bracket [dlo, dhi] with v(dlo) >= v >= v(dhi), bisection (geometric, the
 * distances span many decades) whenever a step leaves it. */
static bool ni_pm_invert(NiPhaseMap* pm, double v, double* dout) {
    double dhi = pm->d0, dlo = pm->dlast;
    double vlo; bool over;
    /* walk the near side of the bracket toward the endpoint */
    for (int k = 0; ; k++) {
        if (k > 2100 || !(dlo > 0.0)) return false;
        if (ni_pm_v(pm, dlo, &vlo, &over)) {
            if (vlo >= v) break;
            dhi = dlo;
        } else if (!over) return false;
        else break;
        dlo *= 0.5;
    }
    double d = sqrt(dlo * dhi);
    for (int it = 0; it < 200; it++) {
        double vd, dv;
        if (!ni_pm_v(pm, d, &vd, &over)) {
            if (!over) return false;
            dlo = d;
        } else {
            if (vd == v) break;
            if (vd > v) dlo = d; else dhi = d;
            if (ni_pm_dv(pm, d, &dv)) {
                double dn = d - (vd - v) / dv;
                if (dn > dlo && dn < dhi) {
                    bool done = fabs(dn - d) <= 4e-16 * d;
                    d = dn;
                    if (done) break;
                    continue;
                }
            }
        }
        if (dhi - dlo <= 4e-16 * dhi) break;
        d = sqrt(dlo * dhi);
        if (!(d > dlo && d < dhi)) d = 0.5 * (dlo + dhi);
    }
    pm->dlast = d;
    *dout = d;
    return true;
}

/* The integrand in the phase coordinate. */
static bool ni_pm_sample(void* vp, double v, double _Complex* out) {
    NiPhaseMap* pm = (NiPhaseMap*)vp;
    double d, dv;
    if (!ni_pm_invert(pm, v, &d) || !ni_pm_dv(pm, d, &dv)) return false;
    double _Complex amp;
    if (!ni_sample_machine(&pm->camp, pm->end + pm->side * d, &amp)) return false;
    double u = pm->s * v;
    /* e^{iu} is integrated as its cosine and sine parts (pm->kind is COS or
     * SIN by then), so each half-line integrand has a real lobe structure */
    double _Complex k = (pm->kind == LEVIN_KERNEL_COS) ? cos(u) : sin(u);
    *out = amp * k / fabs(dv);
    return isfinite(creal(*out)) && isfinite(cimag(*out));
}

static NiAtt ni_core_finite(NiCtx* ctx, double a, double b, const NiOpts* o);

/* amp K(g) with K = cos | sin | e^{i .}: the kernel, the phase g, its closed
 * form derivative g' and the amplitude (all owned on success).  An amplitude
 * that oscillates on its own would leave a second phase: no. */
static bool ni_pm_detect(const Expr* body, const char* var, LevinKernel* kind,
                         Expr** g, Expr** gp, Expr** amp) {
    Expr* fb = expr_copy((Expr*)body);
    if (levin_is_head(fb, SYM_Times)) {
        fb = expr_unshare(fb);
        eval_flatten_args(fb, SYM_Times);
    }
    bool ok = levin_detect_kernel(fb, var, kind, g, amp);
    expr_free(fb);
    if (!ok) return false;
    if (nc_has_phase(*amp, var)) { expr_free(*g); expr_free(*amp); return false; }
    *gp = levin_phase_derivative(*g, var);
    if (!*gp) { expr_free(*g); expr_free(*amp); return false; }
    return true;
}

static void ni_pm_teardown(NiPhaseMap* pm) {
    NiCtx* cs[3] = { &pm->cg, &pm->cgp, &pm->camp };
    for (int i = 0; i < 3; i++) {
        autocompiled_free(cs[i]->ac);
        autocompiled_free(cs[i]->ac_prec);
        cs[i]->ac = NULL; cs[i]->ac_prec = NULL;
    }
    expr_free(pm->g); expr_free(pm->gp); expr_free(pm->amp);
    pm->g = pm->gp = pm->amp = NULL;
}

/* The machine phase map of the integrand of `ctx` (machine precision, plain
 * real axis) at endpoint `side` of [a, b]: false (nothing to tear down) when
 * the integrand is not amp K(g) with a phase monotone and unbounded there.
 * The split point x0 = end + side d0 is where the phase is monotone on
 * (end, x0] -- dv/dd of one sign at 48 distances, halving from x0 -- and
 * unbounded at the end (it leaves the double range, or grows past 10^4
 * lobes); tried at half, an eighth and a sixty-fourth of the span. */
static bool ni_pm_setup(NiCtx* ctx, double a, double b, int side, NiPhaseMap* pm) {
    if (!(b > a) || ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return false;
    if (!ni_pm_detect(ctx->body, ctx->bind->name, &pm->kind, &pm->g, &pm->gp, &pm->amp))
        return false;
    /* each derived body compiles on its own (freed by ni_pm_teardown) */
    pm->cg = *ctx;   ni_ctx_rebody(&pm->cg, pm->g);
    pm->cgp = *ctx;  ni_ctx_rebody(&pm->cgp, pm->gp);
    pm->camp = *ctx; ni_ctx_rebody(&pm->camp, pm->amp);
    pm->cg.ac_enabled = pm->cgp.ac_enabled = pm->camp.ac_enabled = true;
    pm->end = side > 0 ? a : b;
    pm->side = side;
    pm->s = 1.0;
    double span = b - a;
    static const double fr[3] = { 0.5, 0.125, 1.0 / 64.0 };
    for (int t = 0; t < 3; t++) {
        double d0 = fr[t] * span, v0, dv;
        bool over;
        pm->s = 1.0;
        if (!ni_pm_v(pm, d0, &v0, &over) || !ni_pm_dv(pm, d0, &dv)) continue;
        /* v must grow toward the endpoint: dv/dd < 0 */
        if (dv > 0.0) { pm->s = -1.0; v0 = -v0; }
        bool mono = true, unbounded = false;
        double vprev = v0, dd = d0;
        for (int j = 1; j <= 48 && mono; j++) {
            dd *= 0.5;
            double vj, dvj;
            if (!ni_pm_v(pm, dd, &vj, &over)) {
                if (over) unbounded = true;
                else mono = false;
                break;
            }
            if (!ni_pm_dv2(pm, dd, &dvj, &over)) {
                if (over) unbounded = true;
                else mono = false;
                break;
            }
            if (!(dvj < 0.0) || !(vj > vprev)) mono = false;
            vprev = vj;
            if (vj - v0 > 2e4 * M_PI) { unbounded = true; break; }
        }
        if (mono && unbounded) {
            pm->d0 = d0;
            pm->dlast = d0;
            return true;
        }
    }
    ni_pm_teardown(pm);
    return false;
}

/* The phase-coordinate rule at endpoint `side` of [a, b]; have=false when the
 * integrand is not amp K(g) with a phase monotone and unbounded there. */
static NiAtt ni_try_phase_coord(NiCtx* ctx, double a, double b, const NiOpts* o,
                                int side) {
    NiAtt out = { false, false, 0.0, INFINITY };
    NiPhaseMap pm;
    if (!ni_pm_setup(ctx, a, b, side, &pm)) return out;
    double v0; bool over;
    if (ni_pm_v(&pm, pm.d0, &v0, &over)) {
        double reltol, abstol;
        ni_machine_tols(o, &reltol, &abstol);
        /* as the exponential-map lobes (ni_endpoint_piece): three digits
         * beyond the goal, to at most 1e-13 */
        double ort = fmax(reltol * 1e-3, 1e-13);
        if (ort > reltol) ort = reltol;
        double _Complex ov; double oerr;
        bool conv;
        if (pm.kind == LEVIN_KERNEL_EXP) {
            /* amp e^{iu} = amp cos u + i amp sin u */
            double _Complex oc, os; double ec, es;
            pm.kind = LEVIN_KERNEL_COS;
            bool cc = osc_integrate_machine(ni_pm_sample, &pm, v0, 0.0, true,
                                            ort, 2000, &oc, &ec);
            pm.kind = LEVIN_KERNEL_SIN;
            pm.dlast = pm.d0;
            bool cs = osc_integrate_machine(ni_pm_sample, &pm, v0, 0.0, true,
                                            ort, 2000, &os, &es);
            conv = cc && cs;
            ov = oc + I * os;
            oerr = ec + es;
        } else {
            conv = osc_integrate_machine(ni_pm_sample, &pm, v0, 0.0, true,
                                         ort, 2000, &ov, &oerr);
        }
        if (isfinite(creal(ov)) && isfinite(cimag(ov)) && isfinite(oerr)) {
            double x0 = pm.end + side * pm.d0;
            NiAtt rest = side > 0 ? ni_core_finite(ctx, x0, b, o)
                                  : ni_core_finite(ctx, a, x0, o);
            if (rest.have) {
                out.have = true;
                out.conv = conv && rest.conv;
                out.val = ov + rest.val;
                out.err = oerr + rest.err;
            }
        }
    }
    ni_pm_teardown(&pm);
    return out;
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
        if (!best.conv) ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, false, -1, -1));
        if (!best.conv) ni_consider(&best, ni_try_osc(ctx, a, b, reltol));
    } else if (o->method == NI_DE) {
        ni_consider(&best, ni_try_de(ctx, a, b, o, reltol));
    } else if (o->method == NI_OSCSING) {
        ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, true, -1, -1));
    } else if (o->method == NI_TRAP || o->method == NI_RIEMANN
               || o->method == NI_NEWTONCOTES) {
        ni_consider(&best, ni_try_ncr(ctx, a, b, o, reltol, abstol));
    } else {
        /* An oscillatory endpoint singularity first, and its rule's verdict
         * is final: Gauss-Kronrod and tanh-sinh under-sample the ever-faster
         * oscillation and can settle on a wrong value (Cos[Log[x]/x] on
         * [0, 1]: 0.32397 for 0.32473, with no warning). */
        if (o->method == NI_AUTO && ctx->map_mode == NI_MAP_AFFINE
            && creal(ctx->x_scale) == 1.0 && cimag(ctx->x_scale) == 0.0
            && ctx->x_shift == 0.0 && b > a) {
            /* A probe that gave up at an endpoint the phase-coordinate rule
             * can take (a phase unbounded there, of closed-form derivative)
             * settles that end; the probe's out-of-range samples are then no
             * reason to discard the machine pass. */
            bool beyond0 = g_ni_beyond;
            int olo = ni_osc_endpoint_machine(ctx, a, b, +1);
            if (s_ni_ose_gaveup && !beyond0) {
                bool b1 = g_ni_beyond;
                g_ni_beyond = false;
                NiAtt pc = ni_try_phase_coord(ctx, a, b, o, +1);
                if (pc.have && !g_ni_beyond) return pc;
                g_ni_beyond = b1 || g_ni_beyond;
            }
            beyond0 = g_ni_beyond;
            int ohi = ni_osc_endpoint_machine(ctx, a, b, -1);
            if (s_ni_ose_gaveup && !beyond0) {
                bool b1 = g_ni_beyond;
                g_ni_beyond = false;
                NiAtt pc = ni_try_phase_coord(ctx, a, b, o, -1);
                if (pc.have && !g_ni_beyond) return pc;
                g_ni_beyond = b1 || g_ni_beyond;
            }
            if (olo || ohi) {
                NiAtt os = ni_try_endpoint_sing(ctx, a, b, o, false, olo, ohi);
                if (os.have) return os;
            }
        }
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
            ni_consider(&best, ni_try_endpoint_sing(ctx, a, b, o, false, -1, -1));
        /* Oscillatory fallback: neither GK (a subdivision per half-period) nor
         * tanh-sinh (assumes endpoint decay) converges on a many-period
         * integrand; integrate between the zeros instead. */
        if (!best.conv && o->method == NI_AUTO)
            ni_consider(&best, ni_try_osc(ctx, a, b, reltol));
    }
    return best;
}

#ifdef USE_MPFR
/* The working-precision value of the real endpoint whose double is d: the
 * spec node that rounds to d, numericalised at o->bits, else d itself (an
 * endpoint the driver made, such as an Exclusions split point).  Taking 2 Pi
 * as its double cut ~2.4e-16 off the period, a 5e-15 relative error in a
 * periodic integral at any WorkingPrecision.
 *
 * A node matches d when EITHER of its doubles does: the correctly rounded one
 * (mpfr_get_d of the working-precision value; the convention of the exact
 * interior singular points, ni_interior_singular_points) or its machine
 * numericalisation, which is how the driver made d from the spec in the first
 * place (ni_num_endpoint).  The two differ whenever the machine evaluation of
 * a nested exact constant is not correctly rounded -- Log[2 + Sqrt[3]] rounds
 * 2 + Sqrt[3] before the Log and lands one ulp off -- and matching on the
 * correctly rounded double alone then fell through to d itself: the range was
 * integrated over a machine endpoint (16 digits) while the result carried the
 * working precision. */
static void ni_endpoint_mpfr(const NiOpts* o, double d, mpfr_t out) {
    mpfr_set_d(out, d, MPFR_RNDN);
    if (!o->spec_nodes) return;
    NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = o->bits;
    mpfr_t re, im;
    mpfr_inits2(o->bits, re, im, (mpfr_ptr)0);
    for (size_t i = 0; i < o->spec_nnodes; i++) {
        Expr* ev = eval_and_free(expr_copy(o->spec_nodes[i]));
        Expr* num = ev ? numericalize(ev, sp) : NULL;
        bool inexact, hit = false;
        if (num && get_approx_mpfr(num, re, im, &inexact) && mpfr_zero_p(im)
            && mpfr_number_p(re)) {
            double cr = mpfr_get_d(re, MPFR_RNDN), md = NAN;
            Expr* mn = numericalize(ev, numeric_machine_spec());
            if (mn && !ni_to_double_real(mn, &md)) md = NAN;
            expr_free(mn);
            if (cr == d || md == d) {
                mpfr_set(out, re, MPFR_RNDN); hit = true;
            } else if (cr == -d || md == -d) {              /* reflected range */
                mpfr_neg(out, re, MPFR_RNDN); hit = true;
            }
            /* An exact endpoint that came back short of the working precision
             * must not let the result claim it.  Short either visibly (fewer
             * bits) or as a machine double padded to the working precision
             * (a head whose N rule only has a machine value): an exact
             * irrational never has a working-precision value that is exactly
             * a double, so a non-rational exact node whose value fits in 53
             * bits was rounded to one on the way. */
            if (hit && !internal_contains_inexact(ev)) {
                long nb = numeric_min_inexact_bits(num);
                if (!is_rational_like(ev) && !mpfr_zero_p(re)
                    && mpfr_min_prec(re) <= 53)
                    nb = 53;
                if (nb > 0 && nb < o->target_bits
                    && (g_ni_endpoint_bits == 0 || nb < g_ni_endpoint_bits))
                    g_ni_endpoint_bits = nb;
            }
        }
        expr_free(num);
        expr_free(ev);
        if (hit) break;
    }
    mpfr_clears(re, im, (mpfr_ptr)0);
}

/* ------------------------------------------------------------------ *
 *  Endpoint boundary layers at WorkingPrecision                       *
 * ------------------------------------------------------------------ *
 *
 * An integrand that climbs by many orders of magnitude into one endpoint
 * (x^x^x on [1, 6] gains 36000 decades, almost all within 1e-5 of x = 6) puts
 * its whole mass in a layer far thinner than the interval.  tanh-sinh does
 * cluster nodes at the endpoint, but in its t variable the layer becomes a
 * double-exponentially sharp step, and resolving that step to 100 digits took
 * step h = 1/512: 10575 samples.
 *
 * The layer is graded instead.  From the midpoint the distance to the endpoint
 * is divided by 16 while log2|f| keeps rising by more than 4 per step (faster
 * than any integrable power singularity, whose rise is < 4 per step); each
 * such point is a breakpoint, the last one sitting where the integrand has
 * become flat on its own scale.  Every piece then holds an integrand that is
 * analytic and changes smoothly on the piece's own length: Clenshaw-Curtis
 * resolves it in 33-257 nested Chebyshev points (tanh-sinh, which needs ~350
 * at 120 digits even for Exp[-u] on [0, 1], is the fallback for a piece still
 * too steep).  The pieces are integrated from the endpoint inwards and each
 * later piece stops once it is below reltol of the running total (an absolute
 * tolerance), so the far pieces, negligible by hundreds of digits, settle at
 * the coarsest level.  100 digits of x^x^x on [1, 6]: 1150 samples, was
 * 10575.  Returns false (nothing written) when no layer is detected. */
#define NI_LAYER_MAX_STEPS 60
#define NI_LAYER_SHIFT 4   /* distance ratio 2^4 per step; rise threshold 4 bits */

static double ni_mpfr_log2f(NiCtx* ctx, const mpfr_t x, mpfr_t vr, mpfr_t vi) {
    if (!ni_sample_mpfr(ctx, x, vr, vi)) return NAN;
    return de_mpfr_log2_l1(vr, vi);
}

static bool ni_mpfr_graded(NiCtx* ctx, const mpfr_t am, const mpfr_t bm,
                           const NiOpts* o, mpfr_t re, mpfr_t im,
                           double* abserr, bool* conv) {
    long bits = o->bits;
    mpfr_t span, half, mid, vr, vi;
    mpfr_inits2(bits, span, half, mid, vr, vi, (mpfr_ptr)0);
    mpfr_sub(span, bm, am, MPFR_RNDN);
    mpfr_div_2ui(half, span, 1, MPFR_RNDN);
    mpfr_add(mid, am, half, MPFR_RNDN);

    /* Breakpoints toward the steeper endpoint, nearest the midpoint first. */
    mpfr_t pts[NI_LAYER_MAX_STEPS];
    int npts = 0, side = 0;
    double lpeak = -INFINITY;
    double lmid = ni_mpfr_log2f(ctx, mid, vr, vi);
    for (int sd = 1; sd >= -1 && !isnan(lmid); sd -= 2) {
        mpfr_t d, x, cand[NI_LAYER_MAX_STEPS];
        mpfr_inits2(bits, d, x, (mpfr_ptr)0);
        mpfr_set(d, half, MPFR_RNDN);
        double lprev = lmid;
        int k = 0;
        /* sum of log2(rise_i / rise_{i+1}) over the steps: a layer's rise
         * c d (15/16) shrinks 16-fold per step, a power law's is constant */
        double rise_prev = 0.0, lratio = 0.0;
        while (k < NI_LAYER_MAX_STEPS) {
            mpfr_div_2ui(d, d, NI_LAYER_SHIFT, MPFR_RNDN);
            if (sd > 0) mpfr_sub(x, bm, d, MPFR_RNDN);
            else        mpfr_add(x, am, d, MPFR_RNDN);
            if (mpfr_equal_p(x, sd > 0 ? bm : am)) break;
            double l = ni_mpfr_log2f(ctx, x, vr, vi);
            if (isnan(l) || !(l - lprev > (NI_LAYER_SHIFT))) break;
            if (k > 0) lratio += log2(rise_prev / (l - lprev));
            rise_prev = l - lprev;
            mpfr_init2(cand[k], bits);
            mpfr_set(cand[k], x, MPFR_RNDN);
            k++;
            lprev = l;
        }
        /* A rise that stays (about) constant step after step is a power law
         * |x - c|^p, p <= -1, cut off at a small distance eps: the near-pole
         * 1/(eps^2 + x^2).  That is a scale feature of width ~1 in the log of
         * the distance, which the recentred rule (ni_mpfr_rc) takes
         * in ~250 samples, where these pieces need up to 257 Chebyshev points
         * each across log16(1/eps) of them.  Only a rise that keeps shrinking
         * (mean ratio >= 4 a step) is a layer. */
        if (k >= 2 && lratio / (k - 1) < 2.0) {
            for (int i = 0; i < k; i++) mpfr_clear(cand[i]);
            k = 0;
        }
        if (k >= 2 && lprev > lpeak) {
            for (int i = 0; i < npts; i++) mpfr_clear(pts[i]);
            for (int i = 0; i < k; i++) {
                mpfr_init2(pts[i], bits);
                mpfr_set(pts[i], cand[i], MPFR_RNDN);
            }
            npts = k; side = sd; lpeak = lprev;
        }
        for (int i = 0; i < k; i++) mpfr_clear(cand[i]);
        mpfr_clears(d, x, (mpfr_ptr)0);
    }
    if (npts == 0) {
        mpfr_clears(span, half, mid, vr, vi, (mpfr_ptr)0);
        return false;
    }

    /* Pieces from the layered endpoint inwards: [pts[K-1], end], [pts[K-2],
     * pts[K-1]], ..., [mid, pts[0]], [far end, mid] (mirrored for side < 0). */
    double reltol = ni_mpfr_reltol(o);
    int levels = ni_mpfr_levels(o);
    mpfr_t lo, hi, pr, pi, tol;
    mpfr_inits2(bits, lo, hi, pr, pi, tol, (mpfr_ptr)0);
    mpfr_set_zero(re, 1); mpfr_set_zero(im, 1);
    bool all = true;
    double err = 0.0;
    CcCache* cc = ccquad_cache_new(bits);
    for (int i = npts; i >= -1; i--) {
        /* inner point of the piece, toward the far end (i = -1: the far end) */
        mpfr_srcptr far_end = (side > 0) ? am : bm, near_end = (side > 0) ? bm : am;
        mpfr_srcptr inner = (i >= 1) ? pts[i - 1] : (i == 0 ? mid : far_end);
        mpfr_srcptr outer = (i == npts) ? near_end : (i >= 0 ? pts[i] : mid);
        if (side > 0) { mpfr_set(lo, inner, MPFR_RNDN); mpfr_set(hi, outer, MPFR_RNDN); }
        else          { mpfr_set(lo, outer, MPFR_RNDN); mpfr_set(hi, inner, MPFR_RNDN); }
        mpfr_srcptr at = NULL;
        if (i < npts) {
            mpfr_hypot(tol, re, im, MPFR_RNDN);
            mpfr_mul_d(tol, tol, reltol / 16.0, MPFR_RNDN);
            if (!mpfr_zero_p(tol)) at = tol;
        }
        double perr;
        /* Clenshaw-Curtis for the smooth piece; tanh-sinh where the piece is
         * still too steep for 256 Chebyshev points */
        if (!ccquad_mpfr(ni_sample_mpfr, ctx, lo, hi, bits, reltol, at, cc,
                         pr, pi, &perr)
            && !denint_tanhsinh_mpfr_abs(ni_sample_mpfr, ctx, lo, hi, bits, reltol,
                                         at, levels, pr, pi, &perr))
            all = false;
        if (!mpfr_number_p(pr) || !mpfr_number_p(pi)) { all = false; break; }
        mpfr_add(re, re, pr, MPFR_RNDN);
        mpfr_add(im, im, pi, MPFR_RNDN);
        err += perr;
    }
    ccquad_cache_free(cc);
    for (int i = 0; i < npts; i++) mpfr_clear(pts[i]);
    mpfr_clears(span, half, mid, vr, vi, lo, hi, pr, pi, tol, (mpfr_ptr)0);
    *abserr = err;
    *conv = all;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Strong algebraic endpoint singularities at WorkingPrecision         *
 * ------------------------------------------------------------------ *
 *
 * f ~ C |x - a|^p with -1 < p << 0 keeps a share ~ d^(1+p) of its integral
 * within distance d of the endpoint: for x^(-999/1000) E^-x on [0, 1] half of
 * the integral lies below x = 10^-301 and a tenth below 10^-1000.  No rule
 * that samples x can reach that (tanh-sinh's nodes stop near 10^-150 at 40
 * digits, a double stops at 10^-308), so the singular power is mapped away:
 * x = a + span u^(1/q), q = 1 + p, turns f dx into a function of u that is
 * regular at u = 0 (f(x) span u^(1/q - 1)/q -> C span^q/q), and tanh-sinh
 * integrates that.  p is read off log2|f| at three distances 2^-k1 span >
 * 2^-k2 span > 2^-k3 span from the endpoint (consistent slopes required), so
 * only a clean power law qualifies; an inexact p leaves a mild u^(p - p')
 * remainder, which tanh-sinh absorbs. */
typedef struct {
    NiCtx*      c;
    mpfr_srcptr end;    /* the singular endpoint                         */
    mpfr_srcptr span;   /* b - a                                         */
    int         side;   /* +1: endpoint a (x = a + ...), -1: endpoint b  */
    mpfr_t      invq;   /* 1/q                                           */
    mpfr_t      jac;    /* span/q                                        */
    mpfr_t      x, w;   /* scratch                                       */
} NiPowMap;

static bool ni_sample_powmap(void* vctx, const mpfr_t u, mpfr_t out_re, mpfr_t out_im) {
    NiPowMap* m = (NiPowMap*)vctx;
    if (mpfr_sgn(u) <= 0) return false;
    mpfr_pow(m->w, u, m->invq, MPFR_RNDN);               /* u^(1/q) */
    /* Next to a non-zero endpoint the abscissa carries the extra bits that
     * keep its distance from the endpoint exact (1 - u^10 at 164 bits would
     * leave 1 - x with only 164 - 33 log2(1/u) good bits, and the singular
     * factor (1 - x)^p turns that into noise); the integrand then sees the
     * true distance whatever precision it computes at. */
    mpfr_prec_t xp = mpfr_get_prec(m->w);
    if (!mpfr_zero_p(m->end) && mpfr_regular_p(m->w)) {
        long gap = (long)mpfr_get_exp(m->end)
                 - ((long)mpfr_get_exp(m->w) + (long)mpfr_get_exp(m->span));
        if (gap > 0) xp += (mpfr_prec_t)gap + 8;
        /* Capped at 8x the working precision: a node closer than that to the
         * endpoint is dropped (tail truncation), as every sampled rule drops
         * the abscissae that round onto a singular endpoint. */
        if (xp > 8 * mpfr_get_prec(m->w)) xp = 8 * mpfr_get_prec(m->w);
    }
    mpfr_set_prec(m->x, xp);
    mpfr_mul(m->x, m->w, m->span, MPFR_RNDN);
    if (m->side > 0) mpfr_add(m->x, m->end, m->x, MPFR_RNDN);
    else             mpfr_sub(m->x, m->end, m->x, MPFR_RNDN);
    if (mpfr_equal_p(m->x, m->end)) return false;         /* below the precision */
    if (!ni_sample_mpfr(m->c, m->x, out_re, out_im)) return false;
    mpfr_div(m->w, m->w, u, MPFR_RNDN);                   /* u^(1/q - 1) */
    mpfr_mul(m->w, m->w, m->jac, MPFR_RNDN);
    mpfr_mul(out_re, out_re, m->w, MPFR_RNDN);
    mpfr_mul(out_im, out_im, m->w, MPFR_RNDN);
    return true;
}

/* The power p of a clean |x - end|^p law at an endpoint, read off log|f| at
 * three distances (in MPFR: the slope must be good to the working precision,
 * since an error e in q = 1 + p leaves a factor u^(e/q) in the mapped
 * integrand); false when the two slopes disagree (not a power law there) or a
 * probe fails. */
static bool ni_mpfr_endpoint_power(NiCtx* ctx, mpfr_srcptr end, mpfr_srcptr span,
                                   int side, long bits, mpfr_t p_out) {
    /* the probes must stay representable next to a non-zero endpoint */
    long room = bits - 8;
    if (!mpfr_zero_p(end))
        room = bits - 8 - ((long)mpfr_get_exp(end) - (long)mpfr_get_exp(span));
    long k3 = room < 240 ? room : 240;
    if (k3 < 48) return false;
    long ks[3] = { k3 / 4, k3 / 2, k3 };
    mpfr_t x, vr, vi, l[3], s1, s2;
    mpfr_inits2(bits, x, vr, vi, l[0], l[1], l[2], s1, s2, (mpfr_ptr)0);
    bool ok = true;
    for (int i = 0; i < 3 && ok; i++) {
        mpfr_div_2si(x, span, ks[i], MPFR_RNDN);
        if (side > 0) mpfr_add(x, end, x, MPFR_RNDN); else mpfr_sub(x, end, x, MPFR_RNDN);
        if (!ni_sample_mpfr(ctx, x, vr, vi)) { ok = false; break; }
        mpfr_hypot(vr, vr, vi, MPFR_RNDN);
        if (!mpfr_regular_p(vr)) { ok = false; break; }
        mpfr_log2(l[i], vr, MPFR_RNDN);
    }
    if (ok) {
        mpfr_sub(s1, l[1], l[0], MPFR_RNDN);
        mpfr_div_si(s1, s1, -(ks[1] - ks[0]), MPFR_RNDN);   /* p from probes 0, 1 */
        mpfr_sub(s2, l[2], l[1], MPFR_RNDN);
        mpfr_div_si(s2, s2, -(ks[2] - ks[1]), MPFR_RNDN);   /* p from probes 1, 2 */
        double p1 = mpfr_get_d(s1, MPFR_RNDN), p2 = mpfr_get_d(s2, MPFR_RNDN);
        ok = fabs(p1 - p2) <= 1e-3 * (1.0 + fabs(p2));
        if (ok) mpfr_set(p_out, s2, MPFR_RNDN);
    }
    mpfr_clears(x, vr, vi, l[0], l[1], l[2], s1, s2, (mpfr_ptr)0);
    return ok;
}

static bool ni_mpfr_powsing(NiCtx* ctx, const mpfr_t am, const mpfr_t bm,
                            const NiOpts* o, mpfr_t re, mpfr_t im,
                            double* abserr, bool* conv) {
    long bits = o->bits;
    mpfr_t span;
    mpfr_init2(span, bits);
    mpfr_sub(span, bm, am, MPFR_RNDN);
    bool done = false;
    for (int side = 1; side >= -1 && !done; side -= 2) {
        mpfr_t pm;
        mpfr_init2(pm, bits);
        if (!ni_mpfr_endpoint_power(ctx, side > 0 ? am : bm, span, side, bits, pm)) {
            mpfr_clear(pm);
            continue;
        }
        double p = mpfr_get_d(pm, MPFR_RNDN);
        if (!(p < -0.5 && p > -1.0)) { mpfr_clear(pm); continue; }
        NiPowMap m;
        m.c = ctx; m.end = side > 0 ? am : bm; m.span = span; m.side = side;
        mpfr_inits2(bits, m.invq, m.jac, m.x, m.w, (mpfr_ptr)0);
        mpfr_add_ui(m.invq, pm, 1, MPFR_RNDN);               /* q */
        mpfr_clear(pm);
        mpfr_div(m.jac, span, m.invq, MPFR_RNDN);            /* span/q */
        mpfr_ui_div(m.invq, 1, m.invq, MPFR_RNDN);           /* 1/q */
        mpfr_t u0, u1;
        mpfr_inits2(bits, u0, u1, (mpfr_ptr)0);
        mpfr_set_ui(u0, 0, MPFR_RNDN);
        mpfr_set_ui(u1, 1, MPFR_RNDN);
        *conv = denint_tanhsinh_mpfr(ni_sample_powmap, &m, u0, u1, bits,
                                     ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                     re, im, abserr);
        /* side < 0: x = b - span u^(1/q) runs from b to a as u goes 0 -> 1,
         * and dx = -(span/q) u^(1/q-1) du, so the orientation cancels */
        mpfr_clears(u0, u1, m.invq, m.jac, m.x, m.w, (mpfr_ptr)0);
        done = mpfr_number_p(re) && mpfr_number_p(im);
    }
    mpfr_clear(span);
    return done;
}

/* The integrand as the machine-precision geometry samplers of the MPFR
 * rules see it: every arbitrary-precision number replaced by its double (a
 * new tree; NULL when there is none, and `e` itself serves).  The zero
 * finder and the endpoint probe need a few digits, but a 60-digit parameter
 * (a = SetPrecision[0.64, 60] in x^a Sin[a/(2 - x)], an optimiser's
 * iterate) made each of their samples an MPFR evaluation through the
 * interpreter, and the 50-digit integral 70 times slower than with a = 16/25.
 * Nothing under a holding head is touched. */
static Expr* ni_machine_twin(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_MPFR) return expr_new_real(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
    if (e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        uint32_t at = symtab_get_def(h->data.symbol.name)->attributes;
        if (at & (ATTR_HOLDALL | ATTR_HOLDFIRST | ATTR_HOLDREST | ATTR_HOLDALLCOMPLETE))
            return NULL;
    }
    size_t n = e->data.function.arg_count;
    Expr** a = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* y = ni_machine_twin(e->data.function.args[i]);
        if (y && !a) {
            a = malloc(sizeof(Expr*) * n);
            for (size_t k = 0; k < i; k++) a[k] = expr_copy((Expr*)e->data.function.args[k]);
        }
        if (a) a[i] = y ? y : expr_copy((Expr*)e->data.function.args[i]);
    }
    if (!a) return NULL;
    Expr* out = expr_new_function(expr_copy((Expr*)h), a, n);
    free(a);
    return out;
}

/* A machine geometry context for an MPFR rule: the machine spec, no endpoint
 * guard, the machine twin of the body (held in *twin, freed by the caller
 * with ni_geom_done) and the compiled fast path. */
static void ni_geom_ctx(NiCtx* g, const NiCtx* ctx, Expr** twin) {
    *g = *ctx;
    *twin = ni_machine_twin(ctx->body);
    ni_ctx_rebody(g, *twin ? *twin : ctx->body);
    g->spec = numeric_machine_spec();
    g->gd_on = false;
    g->ac_enabled = true;
}

static void ni_geom_done(NiCtx* g, Expr* twin) {
    autocompiled_free(g->ac);
    autocompiled_free(g->ac_prec);
    g->ac = NULL; g->ac_prec = NULL;
    expr_free(twin);
}

/* ------------------------------------------------------------------ *
 *  Oscillatory endpoint singularities at WorkingPrecision              *
 * ------------------------------------------------------------------ *
 *
 * x^-1 Cos[Log[x]/x] on [0, 1] oscillates infinitely often as x -> 0 (with a
 * growing amplitude; the integral converges only conditionally).  No rule
 * that samples x resolves that: tanh-sinh, Clenshaw-Curtis and the graded
 * pieces each see an under-sampled oscillation, and the cascade returned
 * -0.1186 for the true 0.3234 at WorkingPrecision -> 60.  The machine path
 * has the right rule (ni_try_endpoint_sing): the exponential map
 * x = a + span e^-t spreads the singular endpoint over a half line, where the
 * oscillation becomes a sequence of lobes between successive zeros, each
 * integrated on its own and the partial sums extrapolated by Wynn's epsilon.
 * This is its MPFR counterpart (osc_integrate_mpfr): the zeros are located
 * with the machine sampler (a lobe boundary needs only a few digits, and
 * consecutive lobes share it, so the lobes tile the half line exactly), each
 * lobe is integrated at the working precision, and the extrapolated sum
 * decides convergence.  Once such an endpoint is detected this rule's result
 * is the answer: a rule that cannot see the oscillation cannot certify a
 * value either, so an extrapolation that does not settle is returned
 * unconverged (NIntegrate::accgl), never replaced by one of theirs.
 *
 * Detection (ni_osc_endpoint_probe) samples the guarded MPFR integrand, so
 * cancellation noise next to a removable singularity is not mistaken for an
 * oscillation. */
#define NI_OSE_PANELS   2000

typedef struct {
    NiCtx*      c;
    mpfr_srcptr end, span;
    int         side;
    mpfr_t      x, vr, vi;
} NiOseMpfrProbe;

static int ni_ose_mpfr_sign(void* vp, double d, int* sr, int* si) {
    NiOseMpfrProbe* q = (NiOseMpfrProbe*)vp;
    mpfr_set_d(q->x, d, MPFR_RNDN);
    mpfr_mul(q->x, q->x, q->span, MPFR_RNDN);
    if (q->side > 0) mpfr_add(q->x, q->end, q->x, MPFR_RNDN);
    else             mpfr_sub(q->x, q->end, q->x, MPFR_RNDN);
    if (!ni_sample_mpfr(q->c, q->x, q->vr, q->vi)) return 0;
    *sr = mpfr_sgn(q->vr); *si = mpfr_sgn(q->vi);
    return 1;
}

static bool ni_mpfr_osc_endpoint(NiCtx* ctx, mpfr_srcptr end, mpfr_srcptr span,
                                 int side, long bits) {
    /* the deep shell must stay representable next to a non-zero endpoint */
    if (!mpfr_zero_p(end)
        && (long)mpfr_get_exp(end) - (long)mpfr_get_exp(span) > bits - 48)
        return false;
    /* Machine samples first, MPFR to confirm: an MPFR sample next to the
     * endpoint can cost without bound (Cos[Exp[1/x]] at x = 10^-11 asks for
     * the cosine of a number with 10^11 integer bits), where the double
     * overflows at once and declines. */
    NiCtx gm; Expr* twin;
    ni_geom_ctx(&gm, ctx, &twin);
    double ed = mpfr_get_d(end, MPFR_RNDN), sd = mpfr_get_d(span, MPFR_RNDN);
    bool mhit = side > 0 ? ni_osc_endpoint_machine(&gm, ed, ed + sd, +1)
                         : ni_osc_endpoint_machine(&gm, ed - sd, ed, -1);
    ni_geom_done(&gm, twin);
    if (!mhit) return false;
    NiOseMpfrProbe q;
    q.c = ctx; q.end = end; q.span = span; q.side = side;
    mpfr_inits2(bits, q.x, q.vr, q.vi, (mpfr_ptr)0);
    bool hit = ni_osc_endpoint_probe(ni_ose_mpfr_sign, &q);
    mpfr_clears(q.x, q.vr, q.vi, (mpfr_ptr)0);
    return hit;
}

/* The exponential endpoint map in MPFR: t >= 0 -> x = end + side span e^-t,
 * dx = span e^-t dt in both orientations (see ni_try_endpoint_sing). */
typedef struct {
    NiCtx*      c;
    mpfr_srcptr end, span;
    int         side;
    mpfr_t      off, x;   /* scratch */
} NiExpMap;

static bool ni_sample_expmap(void* vctx, const mpfr_t t, mpfr_t out_re, mpfr_t out_im) {
    NiExpMap* m = (NiExpMap*)vctx;
    if (mpfr_sgn(t) < 0) return false;
    mpfr_neg(m->off, t, MPFR_RNDN);
    mpfr_exp(m->off, m->off, MPFR_RNDN);
    mpfr_mul(m->off, m->off, m->span, MPFR_RNDN);          /* span e^-t */
    if (!mpfr_regular_p(m->off)) return false;
    /* as ni_sample_powmap: extra bits keep the distance to a non-zero
     * endpoint exact, capped at 8x the working precision */
    mpfr_prec_t wp = mpfr_get_prec(m->off), xp = wp;
    if (!mpfr_zero_p(m->end)) {
        long gap = (long)mpfr_get_exp(m->end) - (long)mpfr_get_exp(m->off);
        if (gap > 0) xp += (mpfr_prec_t)gap + 8;
        if (xp > 8 * wp) xp = 8 * wp;
    }
    mpfr_set_prec(m->x, xp);
    if (m->side > 0) mpfr_add(m->x, m->end, m->off, MPFR_RNDN);
    else             mpfr_sub(m->x, m->end, m->off, MPFR_RNDN);
    if (mpfr_equal_p(m->x, m->end)) return false;
    if (!ni_sample_mpfr(m->c, m->x, out_re, out_im)) return false;
    mpfr_mul(out_re, out_re, m->off, MPFR_RNDN);
    mpfr_mul(out_im, out_im, m->off, MPFR_RNDN);
    return true;
}

/* One piece [end, end + span] (side > 0) or [end - span, end] (side < 0)
 * whose endpoint `end` is the oscillatory singularity, added to (re, im).
 * Returns whether the extrapolation settled; *any whether a lobe was
 * integrated at all (no oscillation geometry found: nothing is added). */
static bool ni_mpfr_oscsing_piece(NiCtx* ctx, mpfr_srcptr end, mpfr_srcptr span,
                                  int side, const NiOpts* o, mpfr_t re, mpfr_t im,
                                  double* abserr, bool* any) {
    long bits = o->bits;
    NiCtx g; Expr* twin;                /* geometry: machine samples, same map */
    ni_geom_ctx(&g, ctx, &twin);
    g.map_mode = side > 0 ? NI_MAP_EXP_LO : NI_MAP_EXP_HI;
    g.end_a = mpfr_get_d(end, MPFR_RNDN);
    g.span = mpfr_get_d(span, MPFR_RNDN);
    NiExpMap m;
    /* The lobe samples bypass the endpoint cancellation guard: the map keeps
     * the distance to the endpoint exact by construction (the abscissa
     * carries the extra bits), and the guard, built for the affine nodes,
     * declined samples of an integrand with an inexact parameter there
     * (x^a Sin[a/(2 - x)] with a 60-digit a), leaving their zeros at double
     * precision and the extrapolation floored at 1e-33 for all 2000 lobes. */
    NiCtx cl = *ctx;
    cl.gd_on = false;
    m.c = &cl; m.end = end; m.span = span; m.side = side;
    mpfr_inits2(bits, m.off, m.x, (mpfr_ptr)0);
    mpfr_t pr, pi;
    mpfr_inits2(bits, pr, pi, (mpfr_ptr)0);
    double perr = INFINITY;
    bool conv = osc_integrate_mpfr(ni_sample_machine, &g, ni_sample_expmap, &m,
                                   0.0, bits, ni_mpfr_reltol(o), NI_OSE_PANELS,
                                   ni_mpfr_levels(o), pr, pi, &perr);
    *any = (isfinite(perr) || !mpfr_zero_p(pr) || !mpfr_zero_p(pi))
        && mpfr_number_p(pr) && mpfr_number_p(pi);
    if (*any) {
        mpfr_add(re, re, pr, MPFR_RNDN);
        mpfr_add(im, im, pi, MPFR_RNDN);
        *abserr += perr;
    } else {
        conv = false;
    }
    mpfr_clears(pr, pi, m.off, m.x, (mpfr_ptr)0);
    ni_geom_done(&g, twin);
    return conv;
}

/* The phase-coordinate rule at WorkingPrecision (see ni_try_phase_coord).
 * The half line runs in w = v - v0 >= 0, v0 = s g(x0) taken at the working
 * precision, so w = 0 is the split point x0 (a double) exactly and the rest
 * [x0, other end] is an ordinary finite piece.  The lobe geometry is the
 * machine phase map's, on the integrand's machine twin; a lobe sample inverts
 * the phase at the working precision by Newton's method from the machine
 * inverse (16 digits, doubling per step). */
typedef struct {
    NiPhaseMap* pm;          /* machine phase map (geometry, start points)  */
    double v0d;              /* v0 as a double, for the geometry samples    */
    NiCtx cg, cgp, camp;     /* MPFR contexts: phase, derivative, amplitude */
    mpfr_srcptr v0;
    long bits;
    mpfr_t v, x, step, gr, gi, dr, di, ar, ai, k;
} NiPhaseMapMp;

static bool ni_pmm_geo(void* vp, double w, double _Complex* out) {
    NiPhaseMapMp* q = (NiPhaseMapMp*)vp;
    return ni_pm_sample(q->pm, q->v0d + w, out);
}

static bool ni_pmm_sample(void* vp, const mpfr_t w, mpfr_t out_re, mpfr_t out_im) {
    NiPhaseMapMp* q = (NiPhaseMapMp*)vp;
    NiPhaseMap* pm = q->pm;
    if (mpfr_sgn(w) < 0) return false;
    mpfr_add(q->v, q->v0, w, MPFR_RNDN);
    double d;
    if (!ni_pm_invert(pm, mpfr_get_d(q->v, MPFR_RNDN), &d)) return false;
    mpfr_set_d(q->x, pm->end + pm->side * d, MPFR_RNDN);
    /* Newton on s g(x) = v: x -= (s g - v) / (s g') */
    bool done = false;
    for (int it = 0; it < 12 && !done; it++) {
        if (!ni_sample_mpfr(&q->cg, q->x, q->gr, q->gi)
            || !ni_sample_mpfr(&q->cgp, q->x, q->dr, q->di)
            || mpfr_zero_p(q->dr)) return false;
        if (pm->s < 0) { mpfr_neg(q->gr, q->gr, MPFR_RNDN); mpfr_neg(q->dr, q->dr, MPFR_RNDN); }
        mpfr_sub(q->step, q->gr, q->v, MPFR_RNDN);
        mpfr_div(q->step, q->step, q->dr, MPFR_RNDN);
        mpfr_sub(q->x, q->x, q->step, MPFR_RNDN);
        done = mpfr_zero_p(q->step) || mpfr_zero_p(q->x)
            || (long)mpfr_get_exp(q->x) - (long)mpfr_get_exp(q->step) > q->bits - 4;
    }
    if (!done) return false;
    if (!ni_sample_mpfr(&q->cgp, q->x, q->dr, q->di) || mpfr_zero_p(q->dr)
        || !ni_sample_mpfr(&q->camp, q->x, q->ar, q->ai)) return false;
    mpfr_abs(q->dr, q->dr, MPFR_RNDN);
    /* K(u), u = s v */
    if (pm->kind == LEVIN_KERNEL_COS) mpfr_cos(q->k, q->v, MPFR_RNDN);
    else {
        mpfr_sin(q->k, q->v, MPFR_RNDN);
        if (pm->s < 0) mpfr_neg(q->k, q->k, MPFR_RNDN);
    }
    mpfr_mul(out_re, q->ar, q->k, MPFR_RNDN);
    mpfr_mul(out_im, q->ai, q->k, MPFR_RNDN);
    mpfr_div(out_re, out_re, q->dr, MPFR_RNDN);
    mpfr_div(out_im, out_im, q->dr, MPFR_RNDN);
    return mpfr_number_p(out_re) && mpfr_number_p(out_im);
}

static bool ni_mpfr_finite(NiCtx* ctx, double a, double b, const NiOpts* o,
                           mpfr_t re, mpfr_t im, double* abserr);

/* Returns 1 with the result in (re, im) (*conv the verdict), else 0 (the
 * rule does not apply; nothing written).  `side` as ni_try_phase_coord; a, b
 * the doubles of the range. */
static int ni_mpfr_phase_coord(NiCtx* ctx, double a, double b, int side,
                               const NiOpts* o, mpfr_t re, mpfr_t im,
                               double* abserr, bool* conv) {
    if (o->method != NI_AUTO || !(b > a) || ctx->map_mode != NI_MAP_AFFINE
        || creal(ctx->x_scale) != 1.0 || cimag(ctx->x_scale) != 0.0
        || ctx->x_shift != 0.0) return 0;
    long bits = o->bits;
    NiCtx gm; Expr* twin;
    ni_geom_ctx(&gm, ctx, &twin);
    NiPhaseMap pm;
    int out = 0;
    if (!ni_pm_setup(&gm, a, b, side, &pm)) { ni_geom_done(&gm, twin); return 0; }
    LevinKernel kmp; Expr *gmp = NULL, *gpmp = NULL, *ampmp = NULL;
    if ((pm.kind == LEVIN_KERNEL_COS || pm.kind == LEVIN_KERNEL_SIN)
        && ni_pm_detect(ctx->body, ctx->bind->name, &kmp, &gmp, &gpmp, &ampmp)) {
        if (kmp == pm.kind) {
            NiPhaseMapMp q;
            q.pm = &pm; q.bits = bits;
            q.cg = *ctx;   ni_ctx_rebody(&q.cg, gmp);
            q.cgp = *ctx;  ni_ctx_rebody(&q.cgp, gpmp);
            q.camp = *ctx; ni_ctx_rebody(&q.camp, ampmp);
            q.cg.gd_on = q.cgp.gd_on = q.camp.gd_on = false;
            q.cg.ac_enabled = q.cgp.ac_enabled = q.camp.ac_enabled = true;
            mpfr_inits2(bits, q.v, q.x, q.step, q.gr, q.gi, q.dr, q.di, q.ar, q.ai,
                        q.k, (mpfr_ptr)0);
            mpfr_t v0;
            mpfr_init2(v0, bits);
            double x0 = pm.end + side * pm.d0;
            mpfr_set_d(q.x, x0, MPFR_RNDN);
            if (ni_sample_mpfr(&q.cg, q.x, v0, q.gi) && mpfr_number_p(v0)) {
                if (pm.s < 0) mpfr_neg(v0, v0, MPFR_RNDN);
                q.v0 = v0;
                q.v0d = mpfr_get_d(v0, MPFR_RNDN);
                mpfr_t pr, pi;
                mpfr_inits2(bits, pr, pi, (mpfr_ptr)0);
                double perr = INFINITY;
                pm.dlast = pm.d0;
                bool c1 = osc_integrate_mpfr(ni_pmm_geo, &q, ni_pmm_sample, &q, 0.0,
                                             bits, ni_mpfr_reltol(o), NI_OSE_PANELS,
                                             ni_mpfr_levels(o), pr, pi, &perr);
                if (isfinite(perr) && mpfr_number_p(pr) && mpfr_number_p(pi)) {
                    mpfr_t rr, ri;
                    mpfr_inits2(bits, rr, ri, (mpfr_ptr)0);
                    double rerr = 0.0;
                    bool c2 = side > 0 ? ni_mpfr_finite(ctx, x0, b, o, rr, ri, &rerr)
                                       : ni_mpfr_finite(ctx, a, x0, o, rr, ri, &rerr);
                    if (mpfr_number_p(rr) && mpfr_number_p(ri)) {
                        mpfr_add(re, pr, rr, MPFR_RNDN);
                        mpfr_add(im, pi, ri, MPFR_RNDN);
                        *abserr = perr + rerr;
                        *conv = c1 && c2;
                        out = 1;
                    }
                    mpfr_clears(rr, ri, (mpfr_ptr)0);
                }
                mpfr_clears(pr, pi, (mpfr_ptr)0);
            }
            mpfr_clear(v0);
            mpfr_clears(q.v, q.x, q.step, q.gr, q.gi, q.dr, q.di, q.ar, q.ai, q.k,
                        (mpfr_ptr)0);
            NiCtx* cs[3] = { &q.cg, &q.cgp, &q.camp };
            for (int i = 0; i < 3; i++) {
                autocompiled_free(cs[i]->ac);
                autocompiled_free(cs[i]->ac_prec);
            }
        }
        expr_free(gmp); expr_free(gpmp); expr_free(ampmp);
    }
    ni_pm_teardown(&pm);
    ni_geom_done(&gm, twin);
    return out;
}

/* Returns -1 when neither endpoint is an oscillatory singularity (nothing
 * written), 0 when one is but no lobe could be integrated (nothing written;
 * the caller runs its cascade but must not trust its verdict), else 1 with
 * the result in (re, im) and the extrapolation's verdict in *conv.  With both
 * endpoints singular the range is split at its midpoint. */
static int ni_mpfr_oscsing(NiCtx* ctx, mpfr_srcptr am, mpfr_srcptr bm,
                           const NiOpts* o, mpfr_t re, mpfr_t im,
                           double* abserr, bool* conv) {
    long bits = o->bits;
    mpfr_t span, half;
    mpfr_inits2(bits, span, half, (mpfr_ptr)0);
    mpfr_sub(span, bm, am, MPFR_RNDN);
    s_ni_ose_gaveup = false;
    bool lo = ni_mpfr_osc_endpoint(ctx, am, span, +1, bits);
    bool glo = s_ni_ose_gaveup;
    s_ni_ose_gaveup = false;
    bool hi = ni_mpfr_osc_endpoint(ctx, bm, span, -1, bits);
    bool ghi = s_ni_ose_gaveup;
    int out = -1;
    /* the machine probe gave up at an end (its samples left the double
     * range): the phase-coordinate rule, when the phase is unbounded there */
    if (!lo && !hi && (glo || ghi)
        && ni_mpfr_phase_coord(ctx, mpfr_get_d(am, MPFR_RNDN), mpfr_get_d(bm, MPFR_RNDN),
                               glo ? +1 : -1, o, re, im, abserr, conv)) {
        mpfr_clears(span, half, (mpfr_ptr)0);
        return 1;
    }
    if (lo || hi) {
        mpfr_t sr, si;
        mpfr_inits2(bits, sr, si, (mpfr_ptr)0);
        mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1);
        double err = 0.0;
        bool ok, any;
        if (lo && hi) {
            mpfr_div_2ui(half, span, 1, MPFR_RNDN);
            bool a2 = false;
            ok = ni_mpfr_oscsing_piece(ctx, am, half, +1, o, sr, si, &err, &any);
            if (any) ok = ni_mpfr_oscsing_piece(ctx, bm, half, -1, o, sr, si, &err, &a2) && ok;
            any = any && a2;
        } else if (lo) {
            ok = ni_mpfr_oscsing_piece(ctx, am, span, +1, o, sr, si, &err, &any);
        } else {
            ok = ni_mpfr_oscsing_piece(ctx, bm, span, -1, o, sr, si, &err, &any);
        }
        if (any) {
            mpfr_set(re, sr, MPFR_RNDN); mpfr_set(im, si, MPFR_RNDN);
            *abserr = err; *conv = ok;
            out = 1;
        } else {
            out = 0;
        }
        mpfr_clears(sr, si, (mpfr_ptr)0);
    }
    mpfr_clears(span, half, (mpfr_ptr)0);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Narrow features at WorkingPrecision: the recentred DE rule          *
 * ------------------------------------------------------------------ *
 *
 * A narrow interior peak (exp(-x^2 - b/x^2), x^n e^-x), a near-pole
 * (a Lorentzian 1/((x - x0)^2 + eps^2)), or a scale feature next to an
 * endpoint (1/(eps^2 + x^2)) holds its mass where the standard DE rules are
 * coarse; the rule recentred on it (denint.h) resolves it in a few hundred
 * samples.  When one recentring does not settle -- a second feature the
 * centred map under-resolves (two Lorentzians) -- the range is split at the
 * located feature and each piece is done the same way: the feature is then a
 * scale feature at a piece's endpoint and the other one the piece's own
 * dominant feature.  A piece with no feature takes the standard rule.  Depth
 * NI_RC_SPLIT_DEPTH bounds the pieces at 2^depth.
 *
 * Returns -1 when no feature was located on the whole range (nothing written:
 * the caller's cascade runs), else whether the pieces converged. */
#define NI_RC_SPLIT_DEPTH 3

/* The abscissa of the map variable y: a + e^y, or a + (b - a)/(1 + e^-y). */
static void ni_rc_abscissa(DenRecentreKind kind, mpfr_srcptr a, mpfr_srcptr b,
                           double y, mpfr_t xs) {
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(xs));
    if (kind == DEN_RC_HALFLINE) {
        mpfr_set_d(t, y, MPFR_RNDN);
        mpfr_exp(t, t, MPFR_RNDN);
        mpfr_add(xs, a, t, MPFR_RNDN);
    } else {
        mpfr_set_d(t, -y, MPFR_RNDN);
        mpfr_exp(t, t, MPFR_RNDN);
        mpfr_add_ui(t, t, 1, MPFR_RNDN);
        mpfr_sub(xs, b, a, MPFR_RNDN);
        mpfr_div(xs, xs, t, MPFR_RNDN);
        mpfr_add(xs, a, xs, MPFR_RNDN);
    }
    mpfr_clear(t);
}

/* The standard rule on a piece: tanh-sinh, or exp-sinh on a half line. */
static bool ni_mpfr_rc_standard(NiCtx* ctx, DenRecentreKind kind, mpfr_srcptr a,
                                mpfr_srcptr b, const NiOpts* o,
                                mpfr_t re, mpfr_t im, double* abserr) {
    if (kind == DEN_RC_FINITE)
        return denint_tanhsinh_mpfr(ni_sample_mpfr, ctx, a, b, o->bits,
                                    ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                    re, im, abserr);
    return dequad_halfline_mpfr(ni_sample_mpfr, ctx, a, o->bits, ni_mpfr_reltol(o),
                                ni_mpfr_levels(o), re, im, abserr);
}

static int ni_mpfr_rc(NiCtx* ctx, DenRecentreKind kind, mpfr_srcptr a, mpfr_srcptr b,
                      const NiOpts* o, int depth, mpfr_t re, mpfr_t im, double* abserr) {
    long bits = o->bits;
    double y0, w;
    if (!denint_recentre_locate_mpfr(ni_sample_mpfr, ctx, kind, a, b, bits, &y0, &w)) {
        if (depth == 0) return -1;
        return ni_mpfr_rc_standard(ctx, kind, a, b, o, re, im, abserr);
    }
    /* A broad feature (width >= 1/2 in y) hugging an endpoint (within e^-4
     * of the length, or of the half line's origin) is a scale feature there:
     * the near-pole 1/(eps^2 + x^2), or the half of a peak a split left at
     * the cut.  Alone it settles under the full level budget; when it does
     * not, the piece is cut e^4 further in, which leaves the feature its own
     * piece and the rest only its e^-4 remnant (splitting at its centre would
     * just move the cut).  Any other feature is split at its centre when one
     * recentring has not settled by h = 1/32 -- a single well-centred
     * feature settles by h = 1/16, so it shares the piece with another. */
    bool hugs = w >= 0.5 && (kind == DEN_RC_FINITE ? fabs(y0) >= 4.0 : y0 <= -4.0);
    bool split = depth < NI_RC_SPLIT_DEPTH;
    if (denint_recentre_mpfr(ni_sample_mpfr, ctx, kind, a, b, bits, y0, w,
                             ni_mpfr_reltol(o), split && !hugs ? 5 : 7, re, im, abserr))
        return 1;
    if (!split)
        return depth == 0 ? 0 : ni_mpfr_rc_standard(ctx, kind, a, b, o, re, im, abserr);
    double ys = !hugs ? y0 : (y0 > 0.0 ? y0 - 4.0 : y0 + 4.0);
    mpfr_t xs, r2, i2;
    mpfr_inits2(bits, xs, r2, i2, (mpfr_ptr)0);
    ni_rc_abscissa(kind, a, b, ys, xs);
    int ok = 0;
    if (mpfr_number_p(xs) && mpfr_greater_p(xs, a) && (!b || mpfr_less_p(xs, b))) {
        double e1 = 0.0, e2 = 0.0;
        int c1 = ni_mpfr_rc(ctx, DEN_RC_FINITE, a, xs, o, depth + 1, re, im, &e1);
        int c2 = ni_mpfr_rc(ctx, kind, xs, b, o, depth + 1, r2, i2, &e2);
        mpfr_add(re, re, r2, MPFR_RNDN);
        mpfr_add(im, im, i2, MPFR_RNDN);
        *abserr = e1 + e2;
        ok = (c1 == 1 && c2 == 1);
    }
    mpfr_clears(xs, r2, i2, (mpfr_ptr)0);
    return ok;
}

/* ∫_a^b f at the MPFR working precision, ctx already in MPFR mode: an explicit
 * fixed rule, Levin collocation for a detected oscillatory kernel, else the
 * endpoint-robust tanh-sinh rule.  (re, im) are initialised by the caller at
 * o->bits; returns whether the rule converged, with its level difference in
 * *abserr. */
static bool ni_mpfr_finite(NiCtx* ctx, double a, double b, const NiOpts* o,
                           mpfr_t re, mpfr_t im, double* abserr) {
    long bits = o->bits;
    mpfr_t am, bm;
    mpfr_init2(am, bits); mpfr_init2(bm, bits);
    ni_endpoint_mpfr(o, a, am);
    ni_endpoint_mpfr(o, b, bm);
    bool conv;
    *abserr = 0.0;
    /* Endpoint cancellation guard for every rule below (ni_sample_mpfr_guarded);
     * only for the plain real abscissa map, whose x is the integration variable. */
    bool gd_saved = ctx->gd_on;
    mpfr_srcptr gda_saved = ctx->gd_a, gdb_saved = ctx->gd_b;
    Expr* gde_saved = ctx->gd_exact;
    Expr* gd_widened = NULL;
    int gs0 = ctx->gd_state[0], gs1 = ctx->gd_state[1];
    double gr0 = ctx->gd_rate[0], gr1 = ctx->gd_rate[1];
    double gl0 = ctx->gd_L1[0], gl1 = ctx->gd_L1[1];
    double gx0 = ctx->gd_lost1[0], gx1 = ctx->gd_lost1[1];
    if (ctx->map_mode == NI_MAP_AFFINE && creal(ctx->x_scale) == 1.0
        && creal(ctx->x_shift) == 0.0) {
        ctx->gd_on = true; ctx->gd_a = am; ctx->gd_b = bm;
        /* The guard's raised-precision evaluations take the integrand as
         * given, not its prenumericalized copy (whose constants carry only
         * the working precision), with every inexact number of it widened
         * exactly: a number holds the arithmetic it enters to its own
         * precision, and a 45-digit 1/2 in (E^x - 1 - x - x^2/2)/x^3 at
         * WorkingPrecision -> 40 left the deep nodes 10^127 for 1/6. */
        Expr* base = (o->pre_body && ctx->body == o->pre_body) ? o->exact_body : ctx->body;
        gd_widened = ni_widen_numbers(base, bits + NI_GD_MAX_EXTRA + 64);
        ctx->gd_exact = gd_widened ? gd_widened : (base != ctx->body ? base : NULL);
        ctx->gd_state[0] = ctx->gd_state[1] = 0;
        ctx->gd_rate[0] = ctx->gd_rate[1] = 0.0;
    }
    if (o->method == NI_TRAP || o->method == NI_RIEMANN
        || o->method == NI_NEWTONCOTES) {
        NcrRuleKind kind = (o->method == NI_TRAP) ? NCR_TRAPEZOIDAL
            : (o->method == NI_NEWTONCOTES) ? NCR_NEWTONCOTES
            : (o->rule_type == 1) ? NCR_RIEMANN_RIGHT
            : (o->rule_type == 2) ? NCR_RIEMANN_MIDPOINT : NCR_RIEMANN_LEFT;
        int max_levels; long max_eval;
        ni_ncr_budget(o, &max_levels, &max_eval);
        conv = ncr_integrate_mpfr(ni_sample_mpfr, ctx, am, bm, kind,
                                  o->nc_points, o->romberg, bits,
                                  ni_mpfr_reltol(o), 0.0, max_levels, max_eval,
                                  re, im, abserr);
    } else {
        bool got = false; conv = false;
        /* Automatic: Levin first only over many periods (>= 4); a weakly
         * oscillatory integrand goes to the quadrature rules first and to Levin
         * only if they do not converge. */
        const double strong = 8.0 * M_PI;
        bool weak = false;
        /* An oscillatory endpoint singularity first: once one is found, no
         * later rule's convergence verdict can be trusted (see
         * ni_mpfr_oscsing). */
        bool untrusted = false;
        if (o->method == NI_AUTO || o->method == NI_OSCSING) {
            int os = ni_mpfr_oscsing(ctx, am, bm, o, re, im, abserr, &conv);
            if (os == 1) got = true;
            else if (os == 0) untrusted = true;
        }
        if (!got && (o->method == NI_LEVIN || o->method == NI_AUTO)) {
            bool lconv = false;
            if (ni_try_levin_collocation_mpfr(ctx, a, b, o,
                                              o->method == NI_AUTO ? strong : 0.0,
                                              &weak, re, im, &lconv)
                && (o->method == NI_LEVIN || lconv)) {
                conv = lconv; got = true;
            }
        }
        if (!got && (o->method == NI_AUTO || o->method == NI_DE))
            got = ni_mpfr_powsing(ctx, am, bm, o, re, im, abserr, &conv);
        if (!got && (o->method == NI_AUTO || o->method == NI_DE))
            got = ni_mpfr_graded(ctx, am, bm, o, re, im, abserr, &conv);
        /* An integrand analytic on the interval's own scale (x PolyLog[3, -x]
         * on [0, 1], the first lobe of Sin[x^3]) converges geometrically under
         * Clenshaw-Curtis: 33-129 nested Chebyshev samples where tanh-sinh,
         * built for endpoint singularities, pays 150-350 at 30-40 digits.  A
         * bounded probe (N <= 128, kept only when it converged) runs once the
         * cheap detectors of a power-law endpoint and of a boundary layer
         * (x^x^x) have found neither; it stops as soon as a doubling gains
         * under 8 bits, so on anything steeper it declines within 33 samples
         * and tanh-sinh runs as before. */
        if (!got && o->method == NI_AUTO) {
            CcCache* cc = ccquad_cache_new(bits);
            double cerr;
            /* asked for every target digit, not just PrecisionGoal: the rule
             * it replaces overshoots the goal to about full precision, and a
             * converging Chebyshev sum gets there within its last doubling */
            long tb = o->target_bits > 0 ? o->target_bits : bits;
            double prt = pow(2.0, -(double)tb);
            if (prt > ni_mpfr_reltol(o) || !(prt > 0.0)) prt = ni_mpfr_reltol(o);
            if (ccquad_mpfr_upto(ni_sample_mpfr, ctx, am, bm, bits, prt,
                                 NULL, cc, 4, re, im, &cerr)) {
                got = true; conv = true; *abserr = cerr;
            }
            ccquad_cache_free(cc);
        }
        /* a narrow peak, near-pole or endpoint scale feature: the DE rule
         * recentred on it */
        if (!got && o->method == NI_AUTO
            && ni_mpfr_rc(ctx, DEN_RC_FINITE, am, bm, o, 0, re, im, abserr) == 1) {
            got = true; conv = true;
        }
        if (!got)
            conv = denint_tanhsinh_mpfr(ni_sample_mpfr, ctx, am, bm, bits,
                                        ni_mpfr_reltol(o), ni_mpfr_levels(o),
                                        re, im, abserr);
        if (!conv && weak) {
            mpfr_t lre, lim;
            mpfr_init2(lre, bits); mpfr_init2(lim, bits);
            bool lconv = false, w2 = false;
            if (ni_try_levin_collocation_mpfr(ctx, a, b, o, 0.0, &w2, lre, lim, &lconv)
                && lconv) {
                mpfr_set(re, lre, MPFR_RNDN); mpfr_set(im, lim, MPFR_RNDN);
                conv = true;
            }
            mpfr_clears(lre, lim, (mpfr_ptr)0);
        }
        if (untrusted) conv = false;
    }
    ctx->gd_on = gd_saved; ctx->gd_a = gda_saved; ctx->gd_b = gdb_saved;
    ctx->gd_exact = gde_saved;
    expr_free(gd_widened);
    ctx->gd_state[0] = gs0; ctx->gd_state[1] = gs1;
    ctx->gd_rate[0] = gr0; ctx->gd_rate[1] = gr1;
    ctx->gd_L1[0] = gl0; ctx->gd_L1[1] = gl1;
    ctx->gd_lost1[0] = gx0; ctx->gd_lost1[1] = gx1;
    mpfr_clears(am, bm, (mpfr_ptr)0);
    return conv;
}
#endif

#ifdef USE_MPFR
/* Machine precision cannot integrate a power singularity |x - a|^p whose
 * mass lies below the double range: a share ~2^(-1022 q), q = 1 + p, sits
 * closer to the endpoint than the smallest normal double, which is half of
 * the integral for x^(-999/1000) E^-x on [0, 1] (the rules returned 131 and
 * 538 for 999.2).  Detected like the MPFR path does (consistent log-log slopes
 * at three distances), such an integral is handed to the MPFR pass, whose
 * power substitution takes the singular factor exactly. */
static bool ni_machine_power_below_range(NiCtx* c, double a, double b) {
    double span = b - a;
    /* A probe sample that leaves the double range is no verdict on the
     * integral: the rules that sample there themselves raise g_ni_beyond, but
     * the probe must not -- 2^-12 from 0 is already past the double range of
     * Cos[Exp[1/x]], which the phase-coordinate rule integrates at machine
     * precision without sampling there (the flag sent it to an MPFR pass that
     * took the cosine of 10^1000-sized phases, 8 s for 1 digit). */
    bool beyond0 = g_ni_beyond;
    bool hit = false;
    for (int side = 1; side >= -1 && !hit; side -= 2) {
        double end = side > 0 ? a : b;
        double l[3];
        const int ks[3] = { 12, 24, 36 };
        bool ok = true;
        for (int i = 0; i < 3 && ok; i++) {
            double x = end + side * ldexp(span, -ks[i]);
            double _Complex v;
            /* an intermediate past the double range (ni_probe_machine) is no
             * power law, and the deeper samples would only cost more */
            if (x == end || ni_probe_machine(c, x, &v) <= 0 || !(cabs(v) > 0.0)) ok = false;
            else l[i] = log2(cabs(v));
        }
        if (!ok) continue;
        double p1 = -(l[1] - l[0]) / 12.0, p2 = -(l[2] - l[1]) / 12.0;
        if (!(fabs(p1 - p2) <= 1e-3 * (1.0 + fabs(p2)))) continue;
        double q = 1.0 + p2;
        if (q > 0.0 && 1022.0 * q < 60.0) hit = true;
    }
    g_ni_beyond = beyond0;
    return hit;
}
#endif

/* The default AccuracyGoal is a ~16-digit ABSOLUTE floor, 10^-15.95: the
 * combined tolerance 10^-a + |I| 10^-p is met as soon as the error is below it.
 * That is a statement about the integrand's units, not about the integral: for
 * an integrand of size 1e-54 (Gamma[z]^3 on a circle about z = -20, a residue
 * of 3.19e-54) the floor is met by any estimate at all, and the rules stopped
 * after their first panel, three digits off.  When the user has not set
 * AccuracyGoal, the floor is therefore taken relative to the integrand's own
 * scale -- max |f| over seven interior points (sampled through the integral's own,
 * compiled, sampler) times the length of the range --
 * whenever that scale is below 1 (a larger integrand keeps the absolute floor,
 * which its relative goal dominates anyway).  The scaled floor plays the same
 * role at every magnitude: it lets an integral that cancels to ~0 against an
 * O(1) integrand stop, without letting a tiny integrand stop early. */
static void ni_scale_accuracy_floor(NiCtx* ctx, double a, double b, NiOpts* o) {
    if (o->acc_goal_user || o->prec_mpfr || !(o->acc_goal > 0.0) || !isfinite(o->acc_goal))
        return;
    double span = fabs(b - a);
    if (!(span > 0.0) || !isfinite(span)) return;
    /* middle first: a sample already of size >= 1/span settles it */
    static const double fr[7] = { 0.5, 0.2113, 0.7887, 0.0731, 0.9269, 0.3473, 0.6527 };
    double mag = 0.0;
    bool any = false;
    arith_warnings_mute_push();
    for (int i = 0; i < 7 && !(mag * span >= 1.0); i++) {
        double _Complex v;
        if (ni_sample_machine(ctx, a + fr[i] * (b - a), &v)) {
            double m = cabs(v);
            if (m > mag) mag = m;
            any = true;
        }
    }
    arith_warnings_mute_pop();
    double scale = mag * span;
    if (!any || !(scale > 0.0) || scale >= 1.0) return;
    o->acc_goal += -log10(scale);
}

static Expr* ni_run_1d_finite_real(Expr* body, const char* var,
                                   double a, double b, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    ctx.ac_enabled = true;

    double sign = 1.0;
    if (a > b) { double t = a; a = b; b = t; sign = -1.0; }
    ctx.ep_lo = a; ctx.ep_hi = b;

#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        long bits = o->bits;
        mpfr_t re, im;
        mpfr_init2(re, bits); mpfr_init2(im, bits);
        double abserr;
        bool conv = ni_mpfr_finite(&ctx, a, b, o, re, im, &abserr);
        if (sign < 0) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
        ni_bind_restore(&bind);
        double rmag = hypot(mpfr_get_d(re, MPFR_RNDN), mpfr_get_d(im, MPFR_RNDN));
        bool in_range = isfinite(rmag) && isfinite(abserr);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(re, im, (mpfr_ptr)0);
        if (out && !conv) ni_warn_goal_mpfr(o, in_range, abserr, rmag);
        autocompiled_free(ctx.ac);        /* MPFR branch returns here, before the */
        autocompiled_free(ctx.ac_prec);   /* machine-path free at the function end */
        return out;
    }
#endif

#ifdef USE_MPFR
    if (o->method == NI_AUTO && !g_ni_escalated
        && ni_machine_power_below_range(&ctx, a, b)) {
        g_ni_beyond = true;               /* re-run on the MPFR path */
        ni_bind_restore(&bind);
        autocompiled_free(ctx.ac);
        autocompiled_free(ctx.ac_prec);
        return NULL;
    }
#endif
    NiOpts so = *o;
    ni_scale_accuracy_floor(&ctx, a, b, &so);
    NiAtt best = ni_core_finite(&ctx, a, b, &so);
    ni_bind_restore(&bind);
    autocompiled_free(ctx.ac);
    autocompiled_free(ctx.ac_prec);

    if (!best.have) return NULL;        /* no numeric estimate at all */
    if (!best.conv) ni_warn_goal(&so, best.err, cabs(best.val));
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

#ifdef USE_MPFR
    /* Real nodes at WorkingPrecision: the segments are real intervals, so each
     * runs on the MPFR finite-interval path (the machine segment rules below
     * would silently compute a high-precision request in doubles). */
    bool real_nodes = true;
    for (int i = 0; i < nnodes; i++)
        if (cimag(nodes[i]) != 0.0) real_nodes = false;
    if (o->prec_mpfr && real_nodes) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        ctx.ac_enabled = true;
        long bits = o->bits;
        mpfr_t tr, ti, re, im;
        mpfr_inits2(bits, tr, ti, re, im, (mpfr_ptr)0);
        mpfr_set_zero(tr, 1); mpfr_set_zero(ti, 1);
        bool conv = true, ok = true, in_range = true;
        double err_sum = 0.0;
        for (int i = 0; i < nnodes - 1 && ok; i++) {
            double a = creal(nodes[i]), b = creal(nodes[i + 1]);
            if (a == b) continue;
            bool rev = a > b;
            double abserr;
            if (!ni_mpfr_finite(&ctx, rev ? b : a, rev ? a : b, o, re, im, &abserr))
                conv = false;
            if (!mpfr_number_p(re) || !mpfr_number_p(im)) { ok = false; break; }
            if (rev) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
            mpfr_add(tr, tr, re, MPFR_RNDN);
            mpfr_add(ti, ti, im, MPFR_RNDN);
            err_sum += abserr;
        }
        ni_bind_restore(&bind);
        double rmag = hypot(mpfr_get_d(tr, MPFR_RNDN), mpfr_get_d(ti, MPFR_RNDN));
        in_range = isfinite(rmag) && isfinite(err_sum);
        Expr* out = ok ? ni_mpfr_build(tr, ti, o->target_bits) : NULL;
        mpfr_clears(tr, ti, re, im, (mpfr_ptr)0);
        if (out && !conv) ni_warn_goal_mpfr(o, in_range, err_sum, rmag);
        autocompiled_free(ctx.ac);
        autocompiled_free(ctx.ac_prec);
        return out;
    }
#endif

    double _Complex total = 0.0;
    double total_err = 0.0;
    bool all_conv = true, any_fail = false;
    for (int i = 0; i < nnodes - 1; i++) {
        double _Complex z0 = nodes[i], z1 = nodes[i + 1];
        if (z1 == z0) continue;                 /* zero-length segment */
        if (cimag(z0) == 0.0 && cimag(z1) == 0.0) {
            /* A real segment is integrated in x itself, not through the affine
             * t-map: the cascade's oscillatory-endpoint rule applies only on
             * the plain real axis, so with the map NIntegrate[x^a Sin[a/(2 - x)],
             * {x, 0, 1, 2}] lost it on [1, 2] and settled on a Gauss-Kronrod
             * value 6e-5 off (a = 0.138) that {x, 0, 2} got right. */
            double ra = creal(z0), rb = creal(z1), sg = 1.0;
            if (ra > rb) { double t = ra; ra = rb; rb = t; sg = -1.0; }
            ctx.x_scale = 1.0;
            ctx.x_shift = 0.0;
            ctx.ep_lo = ra; ctx.ep_hi = rb;
            NiAtt a = ni_core_finite(&ctx, ra, rb, o);
            if (!a.have) { any_fail = true; break; }
            total += sg * a.val;
            total_err += a.err;
            if (!a.conv) all_conv = false;
            continue;
        }
        ctx.x_scale = z1 - z0;
        ctx.x_shift = z0;
        ctx.ep_lo = 0.0; ctx.ep_hi = 1.0;
        NiAtt a = ni_core_finite(&ctx, 0.0, 1.0, o);
        if (!a.have) { any_fail = true; break; }
        total += (z1 - z0) * a.val;
        total_err += cabs(z1 - z0) * a.err;     /* segment error scaled like its value */
        if (!a.conv) all_conv = false;
    }
    ni_bind_restore(&bind);

    if (any_fail) return NULL;
    if (!all_conv) ni_warn_goal(o, total_err, cabs(total));
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
    if (!conv) ni_warn_goal(o, abserr, cabs(val));
    return ni_from_complex_d(val);
}

/* Semi-infinite range via exp-sinh.  `a` is the half-line lower limit in the
 * sampling variable; reflect=true binds the integrand variable to -u (mapping a
 * (-∞, b] range onto [-b, ∞)); the overall result is multiplied by `sign`. */
#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Oscillatory half lines at WorkingPrecision: contour rotation       *
 * ------------------------------------------------------------------ *
 *
 * An integrand on [a, oo) that is a finite sum of linear-phase groups,
 *     f(x) = sum_j A_j(x) e^(lam_j x),  A_j non-oscillatory,
 * (Sin, Cos and complex exponentials through TrigToExp; BesselJ and BesselY
 * through the scaled Hankel functions, nc_phase.c) is integrated without
 * resolving a single oscillation.  Beyond a base point X0 every group with
 * Im lam != 0 is moved onto the vertical ray X0 + i s t, s = sign(Im lam),
 * where its factor e^(lam x) decays like e^(-|Im lam| t):
 *     int_X0^oo A e^(lam x) dx = i s int_0^oo A(X0 + i s t) e^(lam (X0 + i s t)) dt,
 * exact by Cauchy's theorem when A is analytic and of sub-exponential growth in
 * the quarter plane Re x >= X0, s Im x >= 0.  A group with Im lam = 0 stays on
 * the real axis; [a, X0] is an ordinary finite integral of f.
 *
 * Partial-sum extrapolation between the zeros, the previous route, has no
 * model for several frequencies at once: J_0(x) Cos[x/2] (frequencies 1/2 and
 * 3/2) ran past 100 s at 30 digits; and on a single frequency it stopped near
 * 21 digits.  On the rotated rays exp-sinh converges double-exponentially.
 *
 * Analyticity is established, not assumed: each amplitude may contain only
 * sums, products, the placeholder Bessel K of a linear argument, and powers and
 * logarithms whose base is a numeric constant or a polynomial in x (or a
 * non-negative integer power of an allowed expression).  The polynomials' roots
 * (all of them, complex ones included) are the only possible singularities;
 * X0 is placed to the right of every one, and an integrand with such a root on
 * the real axis inside the range is left to the other rules. */
#define NI_ROT_MAX_GROUPS 16
#define NI_SING_MAX_DEGREE 24   /* polynomial singular arguments (also below) */

static Expr* ni_ev(const char* head, Expr** args, size_t n, const char* var);
static Expr* ni_ev2(const char* head, Expr* a, Expr* b, const char* var);
static bool ni_free_of_var(const Expr* e, const char* var) { return !levin_occurs(e, var); }

/* Degree of u as a polynomial in var with numeric coefficients, read off the
 * tree (sums, products, non-negative integer powers, numbers and the numeric
 * constants), or -1 when u is not one.  A structural test, so the singular
 * point scans cost no evaluation unless a root is actually sought. */
static int ni_numeric_poly_degree(const Expr* u, const char* var) {
    if (!u) return -1;
    switch (u->type) {
    case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT: case EXPR_MPFR: return 0;
    case EXPR_SYMBOL: {
        const char* n = u->data.symbol.name;
        if (n == var) return 1;
        return (n == SYM_Pi || n == SYM_E || n == SYM_I || n == SYM_EulerGamma
                || n == SYM_Catalan || n == SYM_GoldenRatio || n == SYM_Degree) ? 0 : -1;
    }
    case EXPR_FUNCTION: break;
    default: return -1;
    }
    const Expr* hd = u->data.function.head;
    if (hd->type != EXPR_SYMBOL) return -1;
    const char* h = hd->data.symbol.name;
    size_t n = u->data.function.arg_count;
    if (h == SYM_Rational || h == SYM_Complex) return 0;
    if (h == SYM_Plus || h == SYM_Times) {
        int d = 0;
        for (size_t i = 0; i < n; i++) {
            int di = ni_numeric_poly_degree(u->data.function.args[i], var);
            if (di < 0) return -1;
            d = (h == SYM_Plus) ? (di > d ? di : d) : d + di;
            if (d > 4 * NI_SING_MAX_DEGREE) return -1;
        }
        return d;
    }
    if (h == SYM_Power && n == 2) {
        const Expr* b = u->data.function.args[0];
        const Expr* e = u->data.function.args[1];
        if (ni_free_of_var(b, var) && ni_free_of_var(e, var))
            return (ni_numeric_poly_degree(b, var) == 0 && ni_numeric_poly_degree(e, var) == 0)
                   ? 0 : -1;                        /* a numeric constant, e.g. 2^(1/3) */
        if (e->type != EXPR_INTEGER || e->data.integer < 0 || e->data.integer > NI_SING_MAX_DEGREE)
            return -1;
        int db = ni_numeric_poly_degree(b, var);
        return db < 0 ? -1 : db * (int)e->data.integer;
    }
    if (ni_free_of_var(u, var)) {
        /* Sqrt[2] is Power; anything else free of var must be a numeric function
         * of numeric constants (Log[2], Sin[1]) */
        for (size_t i = 0; i < n; i++)
            if (ni_numeric_poly_degree(u->data.function.args[i], var) != 0) return -1;
        return 0;
    }
    return -1;
}

typedef struct { double maxre; bool ok; const char* var; int nroots; } NiRotScan;

static void ni_rot_poly_roots(Expr* u, NiRotScan* sc, double a) {
    const char* var = sc->var;
    int deg = ni_numeric_poly_degree(u, var);
    if (deg == 0) return;
    if (deg < 0 || deg > NI_SING_MAX_DEGREE) { sc->ok = false; return; }
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                                 (Expr*[]){ expr_copy(u), expr_new_integer(0) }, 2);
    Expr* sol = ni_ev(SYM_Solve, (Expr*[]){ eq, expr_new_symbol(var) }, 2, var);
    if (!sol || sol->type != EXPR_FUNCTION || sol->data.function.head->type != EXPR_SYMBOL
        || sol->data.function.head->data.symbol.name != SYM_List) {
        expr_free(sol); sc->ok = false; return;
    }
    for (size_t j = 0; j < sol->data.function.arg_count && sc->ok; j++) {
        Expr* rs = sol->data.function.args[j];
        Expr* r = (rs->type == EXPR_FUNCTION && rs->data.function.arg_count == 1)
                  ? rs->data.function.args[0] : NULL;
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != 2) { sc->ok = false; break; }
        Expr* num = numericalize(r->data.function.args[1], numeric_machine_spec());
        double _Complex z;
        if (!num || !ni_to_complex(num, &z)) sc->ok = false;
        else {
            /* a real singularity on the path itself: not this method's */
            if (fabs(cimag(z)) <= 1e-12 * (1.0 + fabs(creal(z))) && creal(z) > a) sc->ok = false;
            if (creal(z) > sc->maxre) sc->maxre = creal(z);
            sc->nroots++;
        }
        expr_free(num);
    }
    expr_free(sol);
}

/* Walk an amplitude: every head must be in the analytic whitelist. */
static void ni_rot_scan(Expr* e, NiRotScan* sc, double a) {
    if (!sc->ok || !e) return;
    if (e->type != EXPR_FUNCTION) return;            /* numbers and symbols */
    if (ni_free_of_var(e, sc->var)) return;         /* a constant subexpression */
    const Expr* hd = e->data.function.head;
    if (hd->type != EXPR_SYMBOL) { sc->ok = false; return; }
    const char* h = hd->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Plus || h == SYM_Times) {
        for (size_t i = 0; i < n; i++) ni_rot_scan(e->data.function.args[i], sc, a);
        return;
    }
    if (h == SYM_Power && n == 2) {
        Expr* b = e->data.function.args[0];
        Expr* p = e->data.function.args[1];
        if (!ni_free_of_var(p, sc->var)) {
            /* b^(var-dependent): only a constant base, e^(g x) with g analytic */
            if (!ni_free_of_var(b, sc->var)) { sc->ok = false; return; }
            ni_rot_scan(p, sc, a);
            return;
        }
        if (p->type == EXPR_INTEGER && p->data.integer >= 0) { ni_rot_scan(b, sc, a); return; }
        ni_rot_poly_roots(b, sc, a);
        return;
    }
    if (h == SYM_Log && n == 1) { ni_rot_poly_roots(e->data.function.args[0], sc, a); return; }
    if (h == SYM_BesselK && n == 2 && ni_free_of_var(e->data.function.args[0], sc->var)) {
        /* K_nu(w): branch point at w = 0 only; w is linear here */
        ni_rot_poly_roots(e->data.function.args[1], sc, a);
        return;
    }
    sc->ok = false;
}

/* int over the ray X0 + i s t (s = +-1) of A e^(lam x), or, for s = 0, over the
 * real half line [X0, oo), at the working precision; adds into (re, im). */
static bool ni_rot_ray(NiBind* bind, const char* var, Expr* A, Expr* lam, int64_t X0, int s,
                       const NiOpts* o, mpfr_t re, mpfr_t im, double* err) {
    ni_bind_clear_temp(symtab_get_def(var));
    Expr* X0e = expr_new_integer(X0);
    Expr* z = s ? ni_ev2(SYM_Plus, expr_copy(X0e),
                         ni_ev2(SYM_Times, expr_new_integer(s),
                                ni_ev2(SYM_Times, expr_new_symbol(SYM_I), expr_new_symbol(var), var),
                                var), var)
                : ni_ev2(SYM_Plus, expr_copy(X0e), expr_new_symbol(var), var);
    expr_free(X0e);
    Expr* Az = ni_ev(SYM_ReplaceAll, (Expr*[]){ expr_copy(A),
                     expr_new_function(expr_new_symbol(SYM_Rule),
                                       (Expr*[]){ expr_new_symbol(var), expr_copy(z) }, 2) },
                     2, var);
    Expr* ph = ni_ev2(SYM_Power, expr_new_symbol(SYM_E), ni_ev2(SYM_Times, expr_copy(lam), z, var), var);
    Expr* G = ni_ev2(SYM_Times, Az, ph, var);
    if (s) G = ni_ev2(SYM_Times, ni_ev2(SYM_Times, expr_new_integer(s), expr_new_symbol(SYM_I), var),
                      G, var);
    if (!G) return false;
    NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = o->bits;
    NiCtx gc;
    ni_ctx_init(&gc, G, bind, sp);
    gc.ac_enabled = true;
    mpfr_t a0, rr, ri;
    mpfr_inits2(o->bits, a0, rr, ri, (mpfr_ptr)0);
    mpfr_set_zero(a0, 1);
    double e = INFINITY;
    bool conv = dequad_halfline_mpfr(ni_sample_mpfr, &gc, a0, o->bits, ni_mpfr_reltol(o) * 1e-2,
                                     ni_mpfr_levels(o), rr, ri, &e);
    bool ok = mpfr_number_p(rr) && mpfr_number_p(ri);
    if (ok) {
        mpfr_add(re, re, rr, MPFR_RNDN);
        mpfr_add(im, im, ri, MPFR_RNDN);
        *err += e;
    }
    mpfr_clears(a0, rr, ri, (mpfr_ptr)0);
    autocompiled_free(gc.ac);
    autocompiled_free(gc.ac_prec);
    expr_free(G);
    return ok && conv;
}

/* Returns true with (re, im, *abserr) set when the rotation applied and
 * converged; false (nothing usable written) otherwise.  ctx is the half-line
 * context (affine, unreflected), bind its variable frame. */
static bool ni_try_rotated_halfline_mpfr(NiCtx* ctx, NiBind* bind, Expr* body, const char* var,
                                         double a, const NiOpts* o,
                                         mpfr_t re, mpfr_t im, double* abserr) {
    if (!nc_has_phase(body, var)) return false;
    ni_bind_clear_temp(symtab_get_def(var));
    double nchk = (a > 0.0 ? a : 0.0) + 1.0;
    bool bad = false;
    Expr* hk = nc_bessel_to_hankel(body, var, nchk, &bad);
    if (bad || !hk) { expr_free(hk); return false; }
    Expr* ex = ni_ev(SYM_Expand, (Expr*[]){ ni_ev(SYM_TrigToExp, (Expr*[]){ hk }, 1, var) }, 1, var);
    if (!ex) return false;
    Expr* lams[NI_ROT_MAX_GROUPS];
    double _Complex lamd[NI_ROT_MAX_GROUPS];
    Expr* amps[NI_ROT_MAX_GROUPS];
    int ng = 0;
    bool ok = true;
    Expr* const* ts; Expr* one[1]; size_t nt;
    if (ex->type == EXPR_FUNCTION && ex->data.function.head->type == EXPR_SYMBOL
        && ex->data.function.head->data.symbol.name == SYM_Plus) {
        ts = ex->data.function.args; nt = ex->data.function.arg_count;
    } else { one[0] = ex; ts = one; nt = 1; }
    for (size_t i = 0; i < nt && ok; i++) {
        Expr* lam = NULL;
        Expr* f = nc_strip_phase(ts[i], var, &lam);
        double _Complex ld = 0.0;
        Expr* ln = lam ? numericalize(lam, numeric_machine_spec()) : NULL;
        bool lok = ln && ni_to_complex(ln, &ld);
        expr_free(ln);
        if (!f || !lok || creal(ld) > 1e-12) { expr_free(f); expr_free(lam); ok = false; break; }
        int j = 0;
        for (; j < ng; j++) if (cabs(lamd[j] - ld) <= 1e-12 * (1.0 + cabs(ld))) break;
        if (j == ng) {
            if (ng == NI_ROT_MAX_GROUPS) { expr_free(f); expr_free(lam); ok = false; break; }
            lams[ng] = lam; lamd[ng] = ld; amps[ng] = f; ng++;
        } else {
            expr_free(lam);
            amps[j] = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ amps[j], f }, 2);
        }
    }
    expr_free(ex);
    int nosc = 0;
    NiRotScan sc; sc.maxre = -INFINITY; sc.ok = ok; sc.var = var; sc.nroots = 0;
    for (int j = 0; j < ng && sc.ok; j++) {
        Expr* ae = nc_hankel_expand(amps[j]);
        expr_free(amps[j]);
        amps[j] = ni_ev(SYM_Together, (Expr*[]){ ae }, 1, var);  /* poles show as denominators */
        if (!amps[j]) { sc.ok = false; break; }
        ni_rot_scan(amps[j], &sc, a);
        if (fabs(cimag(lamd[j])) > 1e-12) nosc++;
    }
    bool done = false;
    if (sc.ok && nosc > 0) {
        double x0 = a + 1.0;
        if (isfinite(sc.maxre) && sc.maxre + 1.0 > x0) x0 = sc.maxre + 1.0;
        if (x0 < nchk) x0 = nchk;
        int64_t X0 = (int64_t)ceil(x0);
        long bits = o->bits;
        mpfr_t tr, ti, pr, pi;
        mpfr_inits2(bits, tr, ti, pr, pi, (mpfr_ptr)0);
        mpfr_set_zero(tr, 1); mpfr_set_zero(ti, 1);
        double err = 0.0;
        bool all = true;
        /* [a, X0]: the integrand itself */
        if ((double)X0 > a) {
            double ferr = 0.0;
            if (!ni_mpfr_finite(ctx, a, (double)X0, o, pr, pi, &ferr)) all = false;
            if (!mpfr_number_p(pr) || !mpfr_number_p(pi)) all = false;
            else { mpfr_add(tr, tr, pr, MPFR_RNDN); mpfr_add(ti, ti, pi, MPFR_RNDN); err += ferr; }
        }
        for (int j = 0; j < ng && all; j++) {
            int s = cimag(lamd[j]) > 1e-12 ? 1 : cimag(lamd[j]) < -1e-12 ? -1 : 0;
            if (!ni_rot_ray(bind, var, amps[j], lams[j], X0, s, o, tr, ti, &err)) all = false;
        }
        if (all) {
            mpfr_set(re, tr, MPFR_RNDN); mpfr_set(im, ti, MPFR_RNDN);
            *abserr = err;
            done = true;
        }
        mpfr_clears(tr, ti, pr, pi, (mpfr_ptr)0);
    }
    for (int j = 0; j < ng; j++) { expr_free(lams[j]); expr_free(amps[j]); }
    return done;
}
#endif /* USE_MPFR */

static Expr* ni_run_1d_halfline(Expr* body, const char* var, double a,
                                bool reflect, double sign, const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    ctx.ac_enabled = true;
    if (reflect) ctx.x_scale = -1.0;

#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        long bits = o->bits;
        mpfr_t am, re, im;
        mpfr_init2(am, bits); mpfr_init2(re, bits); mpfr_init2(im, bits);
        ni_endpoint_mpfr(o, a, am);
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
        /* Any other sum of linear-phase groups (several frequencies, Bessel
         * factors): rotate each group off the real axis. */
        if (!conv && o->method == NI_AUTO && !reflect)
            conv = ni_try_rotated_halfline_mpfr(&ctx, &bind, body, var, a, o, re, im, &abserr);

        /* A narrow or far-off-centre feature (exp(-x^2 - b/x^2), x^n e^-x,
         * a Lorentzian peak): the exp-sinh rule recentred on it (denint.h). */
        if (!conv && o->method == NI_AUTO)
            conv = ni_mpfr_rc(&ctx, DEN_RC_HALFLINE, am, NULL, o, 0, re, im, &abserr) == 1;
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
        double rmag = hypot(mpfr_get_d(re, MPFR_RNDN), mpfr_get_d(im, MPFR_RNDN));
        bool in_range = isfinite(rmag) && isfinite(abserr);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(am, re, im, (mpfr_ptr)0);
        if (out && !conv) ni_warn_goal_mpfr(o, in_range, abserr, rmag);
        autocompiled_free(ctx.ac);        /* MPFR branch returns here, before the */
        autocompiled_free(ctx.ac_prec);   /* machine-path free at the function end */
        return out;
    }
#endif

    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    double _Complex val; double abserr;
    bool conv = false;
    /* A well-aligned Fourier integrand amp(x) {Sin|Cos}[w x] on [0, oo): the
     * Ooura-Mori rule, as on the MPFR path.  Exp-sinh cannot converge on an
     * oscillatory tail, so the cascade below spent every one of its levels
     * (some 3 x 10^4 samples) before integrating between the zeros; this rule
     * needs a few hundred.  Declines (falls through) when it does not settle. */
    if ((o->method == NI_AUTO || o->method == NI_LEVIN)
        && !reflect && sign > 0.0 && a == 0.0) {
        OscDeKind fkind; double fomega;
        if (ni_fourier_detect(body, var, &fkind, &fomega))
            conv = oscde_fourier_machine(ni_sample_machine, &ctx, fomega, fkind,
                                         reltol, 12, 200000, &val, &abserr);
    }
    if (conv) {
        /* done */
    } else if (o->method == NI_LEVIN) {
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
    autocompiled_free(ctx.ac);
    autocompiled_free(ctx.ac_prec);

    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv) ni_warn_goal(o, abserr, cabs(val));
    return ni_from_complex_d(val * sign);
}

/* Doubly-infinite range via sinh-sinh. */
static Expr* ni_run_1d_wholeline(Expr* body, const char* var, double sign,
                                 const NiOpts* o) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    ctx.ac_enabled = true;

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
        double rmag = hypot(mpfr_get_d(re, MPFR_RNDN), mpfr_get_d(im, MPFR_RNDN));
        bool in_range = isfinite(rmag) && isfinite(abserr);
        Expr* out = (mpfr_number_p(re) && mpfr_number_p(im))
                  ? ni_mpfr_build(re, im, o->target_bits) : NULL;
        mpfr_clears(re, im, (mpfr_ptr)0);
        if (out && !conv) ni_warn_goal_mpfr(o, in_range, abserr, rmag);
        autocompiled_free(ctx.ac);        /* MPFR branch returns here, before the */
        autocompiled_free(ctx.ac_prec);   /* machine-path free at the function end */
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
    autocompiled_free(ctx.ac);
    autocompiled_free(ctx.ac_prec);

    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
    if (!conv) ni_warn_goal(o, abserr, cabs(val));
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

/* Box vertex cancellation guard (see ni_mc_vertex_guard). */
#define NI_VG_MAX_DIM   3
#define NI_VG_MAX_VERTS 8                /* 2^NI_VG_MAX_DIM */
#define NI_VG_ZONE_MIN  0x1p-17          /* smallest vertex zone (relative)      */
#define NI_VG_ZONE_MAX  0x1p-4           /* largest vertex zone (relative)       */
#define NI_VG_PROBE     0x1p-26          /* relative depth of the cancellation probe */
#define NI_VG_TRIGGER   0x1p-8           /* a sample this close probes the vertex */
#define NI_VG_BITS      128
/* vg_tau, the pointwise relative accuracy the zone is sized for, as a
 * multiple of the run's relative tolerance: the cubature's regions average a
 * sample's error over their volume, and at the default goal its 2^-17 zone
 * (rounding error ~4e-7 at the edge) is what this gives; the nested rule asks
 * every inner integral for the goal, so its samples need it pointwise. */
#define NI_VG_TAU_CUB   1e4
#define NI_VG_TAU_NEST  10.0

typedef struct {
    Expr* body;
    NiBind* binds;
    size_t d;
    NumericSpec spec;
    AutoCompiled* ac;   /* compiled body f(v_1..v_d); NULL = interpreter */
    AutoCompiled* ac_prec;  /* MPFR-compiled body (spec.mode MPFR), or NULL */
    /* Vertex cancellation guard of the machine box rules (ni_mc_vertex_guard):
     * the box (borrowed), and per vertex whether it is unprobed (0), found
     * clean (1) or cancelling (2), with the zone of a cancelling one (its
     * relative size), and vg_tau, the pointwise relative accuracy the zone is
     * sized for.  vg_on is false for every other caller. */
    bool          vg_on;
    const double* vg_a;
    const double* vg_b;
    double        vg_tau;
    signed char   vg_state[NI_VG_MAX_VERTS];
    double        vg_zone[NI_VG_MAX_VERTS];
    AutoCompiled* vg_ac;    /* body compiled at NI_VG_BITS, or NULL */
    bool          vg_ac_tried;
} NiMcCtx;

/* Compile `body` as a function of the d bound variables (machine reals).  NULL
 * if uncompilable.  Caller frees with autocompiled_free. */
static AutoCompiled* ni_mc_compile(Expr* body, NiBind* binds, size_t d) {
    const Expr** syms = malloc(d * sizeof(*syms));
    Expr** owned = malloc(d * sizeof(*owned));
    if (!syms || !owned) { free(syms); free(owned); return NULL; }
    for (size_t j = 0; j < d; j++) { owned[j] = expr_new_symbol(binds[j].name); syms[j] = owned[j]; }
    AutoCompiled* ac = autocompile_new(body, syms, d);
    for (size_t j = 0; j < d; j++) expr_free(owned[j]);
    free(syms); free(owned);
    return ac;
}

#ifdef USE_MPFR
/* MPFR-sampled cubature point: the abscissae are doubles (the rule's nodes),
 * but the integrand is evaluated at the working precision and only its VALUE is
 * rounded to a double.  That is what an integrand which cancels near a
 * singular vertex needs: 3 - Cos[x] - Cos[y] - Cos[z] at r = 1e-8 is 1.5e-16,
 * all rounding error in doubles and exact to the working precision here. */
static bool ni_mc_sample_mpfr_at(NiMcCtx* c, const double* x, size_t d,
                                 NumericSpec spec, AutoCompiled* prog,
                                 double _Complex* out);

static bool ni_mc_sample_mpfr(NiMcCtx* c, const double* x, size_t d,
                              double _Complex* out) {
    return ni_mc_sample_mpfr_at(c, x, d, c->spec, c->ac_prec, out);
}

static bool ni_mc_sample_mpfr_at(NiMcCtx* c, const double* x, size_t d,
                                 NumericSpec spec, AutoCompiled* prog,
                                 double _Complex* out) {
    Expr* xs[24];
    if (d > 24) return false;
    mpfr_t t;
    mpfr_init2(t, (mpfr_prec_t)spec.bits);
    for (size_t j = 0; j < d; j++) {
        mpfr_set_d(t, x[j], MPFR_RNDN);
        xs[j] = expr_new_mpfr_copy(t);
    }
    mpfr_clear(t);
    Expr* num = NULL;
    if (prog) num = autocompiled_eval_mpfr(prog, (const Expr* const*)xs);
    if (!num) {
        for (size_t j = 0; j < d; j++) ni_bind_set(&c->binds[j], xs[j]);
        eval_clock_bump();
        arith_warnings_mute_push();
        Expr* raw = eval_and_free(expr_copy(c->body));
        arith_warnings_mute_pop();
        num = raw ? numericalize(raw, spec) : NULL;
        expr_free(raw);
    }
    for (size_t j = 0; j < d; j++) expr_free(xs[j]);
    if (!num) return false;
    bool ok = ni_to_complex(num, out) && isfinite(creal(*out)) && isfinite(cimag(*out));
    expr_free(num);
    return ok;
}
#endif

static bool ni_mc_sample_double(NiMcCtx* c, const double* x, size_t d,
                                double _Complex* out);

#ifdef USE_MPFR
/* The sample at NI_VG_BITS (the guard's re-evaluation). */
static bool ni_mc_sample_guard(NiMcCtx* c, const double* x, size_t d,
                               double _Complex* out) {
    NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = NI_VG_BITS;
    if (!c->vg_ac_tried) {
        c->vg_ac_tried = true;
        const Expr* syms[NI_VG_MAX_DIM] = { NULL, NULL, NULL };
        Expr* owned[NI_VG_MAX_DIM] = { NULL, NULL, NULL };
        if (d > NI_VG_MAX_DIM) return false;
        for (size_t j = 0; j < d; j++) {
            owned[j] = expr_new_symbol(c->binds[j].name); syms[j] = owned[j];
        }
        c->vg_ac = autocompile_new_prec(c->body, syms, d, NI_VG_BITS);
        for (size_t j = 0; j < d; j++) expr_free(owned[j]);
    }
    return ni_mc_sample_mpfr_at(c, x, d, sp, c->vg_ac, out);
}

/* Relative discrepancy of the double sample at p against the NI_VG_BITS one
 * (infinite when the double sample fails); false when the MPFR one fails. */
static bool ni_vg_discrepancy(NiMcCtx* c, const double* p, size_t d, double* rel) {
    double _Complex fd = 0.0, fm = 0.0;
    bool okd = ni_mc_sample_double(c, p, d, &fd);
    if (!ni_mc_sample_guard(c, p, d, &fm)) return false;
    double m = cabs(fm);
    *rel = !okd ? INFINITY : (m > 0.0 ? cabs(fd - fm) / m : (cabs(fd) > 0.0 ? INFINITY : 0.0));
    return true;
}

/* The point of the box at relative depth t from vertex v, along the
 * direction `dir` (each coordinate's share of t). */
static void ni_vg_point(const NiMcCtx* c, int v, size_t d, double t,
                        const double* dir, double* p) {
    for (size_t j = 0; j < d; j++) {
        double span = c->vg_b[j] - c->vg_a[j];
        p[j] = ((v >> j) & 1) ? c->vg_b[j] - t * dir[j] * span
                              : c->vg_a[j] + t * dir[j] * span;
    }
}

/* Box vertex cancellation guard.  An integrable singularity at a vertex of
 * the box is where an integrand formula cancels: the lattice Green's function
 * (1 - c)/((1 - c)^2 + s^2), c = (Cos[k1] + Cos[k2])/2, has 1 - c ~ k^2/4
 * against an absolute rounding error of 1e-16, so within 1e-8 of the corner
 * it is exactly 0 in doubles instead of ~4/k^2, and the integral the rules
 * converge to is 1e-8 off (8 digits whatever the goal).  A sample within the
 * vertex's zone in every coordinate is re-evaluated at NI_VG_BITS when that
 * vertex cancels.  Whether it does is probed once per vertex, at depth
 * NI_VG_PROBE on the diagonal, by comparing the double sample with the MPFR
 * one; a vertex where they agree is never checked again, so a smooth corner
 * costs two samples.
 *
 * The zone of a cancelling vertex is sized to the accuracy asked for
 * (vg_tau): the double sample's relative error grows toward the vertex
 * (as 1e-16/k^2 here), and the zone is where it exceeds vg_tau, found by
 * probing outward-in at depths 2^-4 .. 2^-17 (two directions each, for a
 * rounding error that happens to be small at one point), and never smaller
 * than NI_VG_ZONE_MIN.  A 12-digit nested rule needs it: its inner integrals
 * at outer nodes near the vertex are asked for 12 digits each, and rounding
 * noise of 1e-9 just outside a 2^-17 zone cost the 3-D Watson integral its
 * convergence and the 2-D Green's function its 13th digit.
 * Returns true when it produced *out. */
static bool ni_mc_vertex_guard(NiMcCtx* c, const double* x, size_t d,
                               double _Complex* out, bool* ok) {
    int v = 0;
    double dist = 0.0;      /* relative distance to vertex v (max over coordinates) */
    for (size_t j = 0; j < d; j++) {
        double span = c->vg_b[j] - c->vg_a[j];
        if (!(span > 0.0)) return false;
        double lo = (x[j] - c->vg_a[j]) / span, hi = (c->vg_b[j] - x[j]) / span;
        double r = lo;
        if (hi < lo) { v |= 1 << j; r = hi; }
        if (!(r < NI_VG_ZONE_MAX)) return false;
        if (r > dist) dist = r;
    }
    if (c->vg_state[v] == 0) {
        /* probe once the samples come close: a vertex the rule never
         * approaches costs nothing (not even the MPFR compilation) */
        if (!(dist < NI_VG_TRIGGER)) return false;
        const double diag[NI_VG_MAX_DIM] = { 1.0, 1.0, 1.0 };
        const double skew[NI_VG_MAX_DIM] = { 1.0, 0.71, 0.53 };
        double p[NI_VG_MAX_DIM], rel = 0.0;
        ni_vg_point(c, v, d, NI_VG_PROBE, diag, p);
        c->vg_state[v] = (ni_vg_discrepancy(c, p, d, &rel) && rel > 1e-8) ? 2 : 1;
        c->vg_zone[v] = NI_VG_ZONE_MIN;
        if (c->vg_state[v] == 2) {
            for (double t = NI_VG_ZONE_MAX * 0.5; t > NI_VG_ZONE_MIN; t *= 0.5) {
                double r1 = 0.0, r2 = 0.0;
                ni_vg_point(c, v, d, t, diag, p);
                bool o1 = ni_vg_discrepancy(c, p, d, &r1);
                ni_vg_point(c, v, d, t, skew, p);
                bool o2 = ni_vg_discrepancy(c, p, d, &r2);
                if ((o1 && r1 > c->vg_tau) || (o2 && r2 > c->vg_tau)) {
                    c->vg_zone[v] = fmin(2.0 * t, NI_VG_ZONE_MAX);
                    break;
                }
            }
        }
    }
    if (c->vg_state[v] != 2 || !(dist < c->vg_zone[v])) return false;
    *ok = ni_mc_sample_guard(c, x, d, out);
    return true;
}
#endif

static bool ni_mc_sample(void* vctx, const double* x, size_t d, double _Complex* out) {
    NiMcCtx* c = (NiMcCtx*)vctx;
#ifdef USE_MPFR
    if (c->spec.mode == NUMERIC_MODE_MPFR) return ni_mc_sample_mpfr(c, x, d, out);
    if (c->vg_on && d <= NI_VG_MAX_DIM) {
        bool ok = false;
        if (ni_mc_vertex_guard(c, x, d, out, &ok)) return ok;
    }
#endif
    return ni_mc_sample_double(c, x, d, out);
}

static bool ni_mc_sample_double(NiMcCtx* c, const double* x, size_t d,
                                double _Complex* out) {
    /* Compiled fast path: evaluate over the d machine abscissae directly.  A
     * non-finite / non-real compiled result (e.g. the interpreter would return a
     * complex value here) falls through to the interpreter for this point. */
    if (c->ac) {
        double _Complex z;
        if (autocompiled_eval_complex(c->ac, x, &z)
            && isfinite(creal(z)) && isfinite(cimag(z))) {
            *out = z;
            return true;
        }
    }
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
    ctx.ac = ni_mc_compile(body, binds, d);
    ctx.ac_prec = NULL;
    ctx.vg_on = false; ctx.vg_a = ctx.vg_b = NULL; ctx.vg_tau = 0.0;
    memset(ctx.vg_state, 0, sizeof ctx.vg_state);
    ctx.vg_ac = NULL; ctx.vg_ac_tried = false;

    /* Monte-Carlo precision goal defaults low (≈2 digits), per Mathematica. */
    double reltol = (o->prec_goal > 0.0) ? pow(10.0, -o->prec_goal) : 1e-2;
    double abstol = (o->acc_goal > 0.0) ? pow(10.0, -o->acc_goal) : 0.0;
    long maxp = (o->max_points > 0) ? o->max_points : 2000000;

    double _Complex val; double abserr;
    bool conv = mc_integrate_machine(ni_mc_sample, &ctx, a, b, d, quasi,
                                     abstol, reltol, maxp, &val, &abserr);

    for (size_t j = 0; j < d; j++) ni_bind_restore(&binds[j]);
    autocompiled_free(ctx.ac);
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

/* Largest PrecisionGoal (digits) the MPFR-sampled cubature accepts: its rule
 * weights and running sums are doubles, and ~13 digits survive the summation
 * of up to 10^6 samples with room to spare. */
#define NI_CUB_MPFR_MAX_GOAL 12.0
static bool ni_cub_mpfr_goal_ok(const NiOpts* o) {
#ifdef USE_MPFR
    double rt = ni_mpfr_reltol(o);
    return rt >= pow(10.0, -NI_CUB_MPFR_MAX_GOAL);
#else
    (void)o;
    return false;
#endif
}

/* Iterated adaptive quadrature over a constant 2-D or 3-D box (machine
 * precision): the fallback of the cubature for an integrable point or edge
 * singularity that Genz-Malik bisection cannot resolve to the goal.  A
 * lattice Green's function (1 - c)/((1 - c)^2 + s^2), c = (Cos[k1] +
 * Cos[k2])/2, s = 2 e Sin[k1], is singular at k = 0 inside a parabolic cusp
 * |k1| < k2^2/(8 e); the cubature's regions shrink onto the corner and its
 * error falls only algebraically (x8 per doubling of the samples in 2-D, x4
 * in 3-D), so 12 digits would take millions of samples.  Integrated one axis
 * at a time, each inner integral is a peak at (or near) its endpoint and the
 * outer integrand is algebraically singular there (~ k1^-1/2): QUADPACK's
 * QAGS (gk_qags_machine) settles both with its epsilon extrapolation in tens
 * of QK21 panels, and its bisection stops at the extrapolated level instead
 * of going onto the singular point -- which keeps the vertex guard's MPFR
 * re-evaluations few.  This is what nested QUADPACK (scipy's dblquad /
 * tplquad) does.  Every level runs in C on the same compiled sampler, so the
 * cost is that of the samples alone (the iterated path of the evaluator
 * re-enters NIntegrate per outer node).
 *
 * Tolerances: every level is asked for NI_NEST_HEADROOM of the goal (the
 * outer QAGS estimate does not see the inner errors, and its truncation error
 * at the goal left the 2-D Green's function at the edge of 12 digits; the
 * extrapolation makes the extra digit cheap); an inner level is also given an
 * absolute floor, that goal times the integral's scale (the cubature's
 * estimate) over the outer ranges, so an inner integral that cancels to ~0
 * does not chase relative digits the total cannot see.  An inner integral
 * that misses the requested goal (not the headroom) counts against the
 * result by its excess error times the outer volume it stands for (an upper
 * bound on its weight): inner integrals at outer nodes deep in the singular
 * corner need not resolve to their own relative goal.  The result converged
 * when that total error meets the requested goal. */
#define NI_NEST_LIMIT 100      /* QAGS subintervals per 1-D integral */
#define NI_NEST_HEADROOM 0.1   /* the 1-D rules' goal, relative to the box's */
#define NI_NEST_ADAPT_SUBDIV 200  /* bisections of the plain rule behind a QAGS miss */
#define NI_NEST_MAX_DIM NI_VG_MAX_DIM

typedef struct NiNest NiNest;
typedef struct { NiNest* n; size_t axis; } NiNestLevel;
struct NiNest {
    NiMcCtx*    mc;
    size_t      d;
    double      x[NI_NEST_MAX_DIM];
    const double* a;
    const double* b;
    double      abstol[NI_NEST_MAX_DIM], reltol;
    double      miss;            /* largest inner error x outer range over the tolerance */
    double      miss_cap;        /* the box's tolerance: a larger miss ends the run */
    bool        aborted;
    NiNestLevel lev[NI_NEST_MAX_DIM];
};

/* The integrand of level `axis`: the sample itself at the last axis, else the
 * integral over the next axis at the current coordinates. */
static bool ni_nest_level(void* v, double t, double _Complex* out) {
    NiNestLevel* L = (NiNestLevel*)v;
    NiNest* n = L->n;
    if (n->aborted) return false;
    n->x[L->axis] = t;
    if (L->axis + 1 == n->d) return ni_mc_sample(n->mc, n->x, n->d, out);
    size_t j = L->axis + 1;
    GkResult R = gk_qags_machine(ni_nest_level, &n->lev[j], n->a[j], n->b[j],
                                 n->abstol[j], n->reltol, NI_NEST_LIMIT);
    if (R.status == GK_NONNUMERIC) return false;
    if (R.status != GK_OK) {
        /* QAGS reads bisections that change neither value nor error as
         * roundoff and stops (ier = 2) -- which is also what an interior peak
         * narrower than its nodes looks like (the ridge s = 0 leaving the
         * singular point into the box: a Lorentzian of width ~k1^2 at
         * k2 = -k1/2).  The plain adaptive rule keeps bisecting until it
         * finds it; it is asked only for an integral QAGS missed. */
        GkResult R2 = gk_integrate_machine(ni_nest_level, &n->lev[j], n->a[j], n->b[j],
                                           n->abstol[j], n->reltol, NI_NEST_ADAPT_SUBDIV,
                                           0, true);
        if (R2.status == GK_OK || (R2.status == GK_NOCONV && R2.abs_err < R.abs_err)) R = R2;
    }
    if (R.status != GK_OK) {
        /* how much this miss could move the total: its excess error spread
         * over the outer box (the product of the outer ranges) */
        double tol = fmax(n->abstol[j], n->reltol * cabs(R.value)) / NI_NEST_HEADROOM;
        double w = 1.0;
        for (size_t q = 0; q < j; q++) w *= n->b[q] - n->a[q];
        double m = (R.abs_err - tol) * w;
        if (m > n->miss) n->miss = m;
        /* a miss the goal cannot absorb: this rule cannot meet it, so stop
         * paying for it */
        if (n->miss > n->miss_cap) { n->aborted = true; return false; }
    }
    *out = R.value;
    return isfinite(creal(R.value)) && isfinite(cimag(R.value));
}

static CubStatus ni_nest_machine(NiMcCtx* mc, const double* a, const double* b,
                                 size_t d, double abstol, double reltol, double scale,
                                 double _Complex* val, double* abserr) {
    if (d < 2 || d > NI_NEST_MAX_DIM) return CUB_NONNUMERIC;
    NiNest n;
    n.mc = mc; n.d = d; n.a = a; n.b = b;
    /* the 1-D rules are asked for NI_NEST_HEADROOM of the goal (see above) */
    double rq = NI_NEST_HEADROOM * reltol, aq = NI_NEST_HEADROOM * abstol;
    n.reltol = rq;
    n.miss = 0.0;
    n.miss_cap = fmax(abstol, reltol * scale);
    n.aborted = false;
    /* the guard re-probes this box, its zones sized for the inner goals */
    mc->vg_a = a; mc->vg_b = b;
    mc->vg_tau = NI_VG_TAU_NEST * reltol;
    memset(mc->vg_state, 0, sizeof mc->vg_state);
    double outer = 1.0;
    for (size_t j = 0; j < d; j++) {
        n.x[j] = a[j];
        n.lev[j].n = &n; n.lev[j].axis = j;
        /* the absolute floor of level j: the goal over the ranges outside it */
        n.abstol[j] = (j == 0) ? aq : 0.1 * fmax(aq, rq * scale) / outer;
        outer *= b[j] - a[j];
    }
    NiNestLevel top = { &n, 0 };
    GkResult R = gk_qags_machine(ni_nest_level, &top, a[0], b[0], aq, rq,
                                 NI_NEST_LIMIT);
    if (n.aborted) { *val = R.value; *abserr = INFINITY; return CUB_NOCONV; }
    if (R.status == GK_NONNUMERIC) return CUB_NONNUMERIC;
    *val = R.value;
    double tol = fmax(abstol, reltol * cabs(R.value));
    *abserr = R.abs_err + (n.miss > 0.0 ? n.miss : 0.0);
    return (*abserr <= tol) ? CUB_OK : CUB_NOCONV;
}


/* One box of the singular retry: the cubature (with the vertex guard on this
 * box), then for d <= 3 the iterated rule when it does not converge.  With
 * `forecast` the cubature may stop early for the iterated rule; when that
 * does not converge either (a ridge leaving the singular point into the box,
 * which the inner integrals of the iterated rule cannot see, where the
 * cubature's regions follow it) the cubature is run again to its budget. */
static CubStatus ni_cub_box(NiMcCtx* ctx, const double* a, const double* b, size_t d,
                            double abstol, double reltol, long max_eval, bool forecast,
                            double _Complex* val, double* abserr) {
    ctx->vg_a = a; ctx->vg_b = b;
    ctx->vg_tau = NI_VG_TAU_CUB * reltol;
    memset(ctx->vg_state, 0, sizeof ctx->vg_state);
    bool cut = forecast && d <= NI_NEST_MAX_DIM;
    CubStatus st = cub_integrate_machine_hot(ni_mc_sample, ctx, a, b, d, abstol, reltol,
                                             reltol, max_eval, val, abserr, NULL, NULL, NULL,
                                             &cut);
    if (st == CUB_NOCONV && d <= NI_NEST_MAX_DIM) {
        double _Complex v2; double e2;
        CubStatus s2 = ni_nest_machine(ctx, a, b, d, abstol, reltol, cabs(*val), &v2, &e2);
        if (s2 == CUB_OK || (s2 == CUB_NOCONV && e2 < *abserr)) {
            st = s2; *val = v2; *abserr = e2;
        }
    }
    if (st == CUB_NOCONV && cut) {
        double _Complex v3; double e3;
        ctx->vg_a = a; ctx->vg_b = b;
        ctx->vg_tau = NI_VG_TAU_CUB * reltol;
        memset(ctx->vg_state, 0, sizeof ctx->vg_state);
        CubStatus s3 = cub_integrate_machine(ni_mc_sample, ctx, a, b, d, abstol, reltol,
                                             max_eval, &v3, &e3);
        if (s3 == CUB_OK || (s3 == CUB_NOCONV && e3 < *abserr)) {
            st = s3; *val = v3; *abserr = e3;
        }
    }
    return st;
}

/* A machine-precision box the cubature did not resolve.  The refinement piles
 * up at an integrable point singularity; when that point lies inside the box
 * (on the k1 = 0 face of [0, Pi] x [-Pi, Pi], where the lattice Green's
 * function of a walk biased across both axes is singular) the box is split
 * there, so that it is a VERTEX of every piece: the vertex guard then
 * re-evaluates its cancelling samples, and the iterated rule's inner and
 * outer integrals see it at an endpoint, where their bisection goes.  A point
 * already at a vertex of the box needs no split (the guard has it).  Each
 * piece is asked for the relative goal and for its share of the absolute one
 * (a small piece must not chase digits the total does not need); the result
 * replaces the first run's when it converged or has the smaller error. */
static void ni_cub_singular_retry(NiMcCtx* ctx, const double* a, const double* b,
                                  size_t d, double abstol, double reltol,
                                  long max_eval, bool forecast,
                                  const double* hot, const double* hot_w,
                                  CubStatus* st, double _Complex* val, double* abserr) {
    int axes[NI_VG_MAX_DIM]; int k = 0;
    double sp[NI_VG_MAX_DIM];
    if (hot)
        for (size_t j = 0; j < d; j++) {
            double span = b[j] - a[j];
            /* The pile-up sits NEXT to the singular point, not on it (along
             * the ridge s = 0 that leaves it, in the walk biased on both
             * axes: 1e-5 off the 0 of [-Pi, Pi]), so the point is snapped to
             * the coarsest dyadic point of the range (the ends and the middle
             * first) near it -- where a singularity of a symmetric range sits.
             * "Near" is 1e-4, four widths of the smallest region when that is
             * wider, and 2^-(L+8) for a point of level L: a run the forecast
             * stopped early has piled up farther from the point (2.4e-4 off
             * the k1 = 0 face, next to the regions frozen at its cancelling
             * samples), and a coarse dyadic point is the likelier site. */
            double t = (hot[j] - a[j]) / span;
            double near = hot_w ? fmax(1e-4, 4.0 * hot_w[j] / span) : 1e-4;
            sp[j] = hot[j];
            for (int L = 1; L <= 30; L++) {
                double m = floor(t * ldexp(1.0, L) + 0.5);
                if (fabs(t - ldexp(m, -L)) <= fmax(near, ldexp(1.0, -(L + 8)))) {
                    sp[j] = a[j] + span * ldexp(m, -L);
                    break;
                }
            }
            if (sp[j] - a[j] > 1e-9 * span && b[j] - sp[j] > 1e-9 * span) axes[k++] = (int)j;
        }
    hot = sp;
    if (k == 0) {
        /* the point is a vertex (or unknown): the iterated rule on the box */
        if (d <= NI_NEST_MAX_DIM) {
            double _Complex v2; double e2;
            CubStatus s2 = ni_nest_machine(ctx, a, b, d, abstol, reltol, cabs(*val), &v2, &e2);
            if (s2 == CUB_OK || (s2 == CUB_NOCONV && e2 < *abserr)) {
                *st = s2; *val = v2; *abserr = e2;
            }
        }
        return;
    }
    int npieces = 1 << k;
    double share = fmax(abstol, reltol * cabs(*val)) / (double)npieces;
    double _Complex sum = 0.0; double err = 0.0;
    bool all_ok = true;
    for (int m = 0; m < npieces; m++) {
        double pa[NI_VG_MAX_DIM], pb[NI_VG_MAX_DIM];
        for (size_t j = 0; j < d; j++) { pa[j] = a[j]; pb[j] = b[j]; }
        for (int q = 0; q < k; q++) {
            int j = axes[q];
            if ((m >> q) & 1) pa[j] = hot[j]; else pb[j] = hot[j];
        }
        double _Complex v; double e;
        CubStatus sp = ni_cub_box(ctx, pa, pb, d, fmax(abstol, share), reltol,
                                  max_eval, forecast, &v, &e);
        if (sp == CUB_NONNUMERIC) { ctx->vg_a = a; ctx->vg_b = b; return; }
        if (sp != CUB_OK) all_ok = false;
        sum += v; err += e;
    }
    ctx->vg_a = a; ctx->vg_b = b;
    if (all_ok || err < *abserr) {
        *st = all_ok ? CUB_OK : CUB_NOCONV;
        *val = sum; *abserr = err;
    }
}

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
    ctx.ac = NULL; ctx.ac_prec = NULL;
    /* the vertex guard: machine samples only (an MPFR run computes at the
     * working precision already) */
    ctx.vg_on = !o->prec_mpfr && d <= NI_VG_MAX_DIM;
    ctx.vg_a = a; ctx.vg_b = b;
    memset(ctx.vg_state, 0, sizeof ctx.vg_state);
    ctx.vg_ac = NULL; ctx.vg_ac_tried = false;
    double reltol, abstol;
    ni_machine_tols(o, &reltol, &abstol);
    double reltol_deep = -1.0;          /* set below; < 0: as reltol */
#ifdef USE_MPFR
    if (o->prec_mpfr) {
        ctx.spec.mode = NUMERIC_MODE_MPFR;
        ctx.spec.bits = o->bits;
        const Expr** syms = malloc(d * sizeof(*syms));
        Expr** owned = malloc(d * sizeof(*owned));
        if (syms && owned) {
            for (size_t j = 0; j < d; j++) {
                owned[j] = expr_new_symbol(binds[j].name); syms[j] = owned[j];
            }
            ctx.ac_prec = autocompile_new_prec(body, syms, d, o->bits);
            for (size_t j = 0; j < d; j++) expr_free(owned[j]);
        }
        free(syms); free(owned);
        reltol = ni_mpfr_reltol(o);
        abstol = 0.0;
    } else
#endif
    {
        ctx.ac = ni_mc_compile(body, binds, d);
        /* The Automatic goal, 10^-8, is the goal of the RESULT.  The 1-D rules
         * overshoot it by far (Gauss-Kronrod converges geometrically), and so
         * does the cubature on a smooth integrand, but at a singularity the
         * Genz-Malik estimate tracks the true error closely: asked for 10^-8
         * it delivered 10^-8 and no more, 7-8 digits of a lattice Green's
         * function where nested QUADPACK gives 10-13.  Two digits of headroom
         * once the refinement goes deep (reltol_deep, see cubature.h), and
         * for the singular retries, which only run there. */
        if (!(o->prec_goal > 0.0)) reltol_deep = reltol * 1e-2;
    }
    long max_eval = (o->max_points > 0) ? o->max_points : 0;

    double _Complex val; double abserr;
    double hot[NI_VG_MAX_DIM] = { 0.0, 0.0, 0.0 };
    double hot_w[NI_VG_MAX_DIM] = { 0.0, 0.0, 0.0 };
    bool hot_ok = false;
    bool mach = !o->prec_mpfr && o->max_points <= 0 && d <= NI_VG_MAX_DIM;
    if (!(reltol_deep >= 0.0)) reltol_deep = reltol;
    ctx.vg_tau = NI_VG_TAU_CUB * reltol_deep;
    /* mach: the singular retry (face split, iterated rule) follows a run that
     * does not converge, so the cubature may stop early when it cannot
     * finish (cubature.h); when the retry does not converge either, the
     * whole unhurried path -- the cubature to its budget, whose pile-up then
     * locates the singular point closely, and the retry -- runs after it,
     * and the better of the two results is kept. */
    bool cut = mach;
    CubStatus st = cub_integrate_machine_hot(ni_mc_sample, &ctx, a, b, d,
                                             abstol, reltol, reltol_deep, max_eval,
                                             &val, &abserr, mach ? hot : NULL, hot_w, &hot_ok,
                                             &cut);
    if (st == CUB_NOCONV && mach)
        ni_cub_singular_retry(&ctx, a, b, d, abstol, reltol_deep, max_eval, true,
                              hot_ok ? hot : NULL, hot_w, &st, &val, &abserr);
    if (st == CUB_NOCONV && cut) {
        double _Complex v2; double e2;
        ctx.vg_a = a; ctx.vg_b = b;
        ctx.vg_tau = NI_VG_TAU_CUB * reltol_deep;
        memset(ctx.vg_state, 0, sizeof ctx.vg_state);
        CubStatus s2 = cub_integrate_machine_hot(ni_mc_sample, &ctx, a, b, d,
                                                 abstol, reltol, reltol_deep, max_eval,
                                                 &v2, &e2, hot, hot_w, &hot_ok, NULL);
        if (s2 == CUB_NOCONV)
            ni_cub_singular_retry(&ctx, a, b, d, abstol, reltol_deep, max_eval, false,
                                  hot_ok ? hot : NULL, hot_w, &s2, &v2, &e2);
        if (s2 == CUB_NONNUMERIC) st = s2;
        else if (s2 == CUB_OK || e2 < abserr) { st = s2; val = v2; abserr = e2; }
    }

    for (size_t j = 0; j < d; j++) ni_bind_restore(&binds[j]);
    autocompiled_free(ctx.ac);
    autocompiled_free(ctx.ac_prec);
    autocompiled_free(ctx.vg_ac);
    free(a); free(b); free(binds);

    if (st == CUB_NONNUMERIC) return NULL;   /* fall back to iterated quadrature */
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return NULL;
#ifdef USE_MPFR
    if (o->prec_mpfr) {
        /* The sums were accumulated in doubles, which holds the goal (at most
         * NI_CUB_MPFR_MAX_GOAL digits, see the caller); the result carries the
         * requested precision as every WorkingPrecision result does. */
        double rmag = cabs(val);
        if (st == CUB_NOCONV) ni_warn_goal(o, abserr, rmag);
        mpfr_t re, im;
        mpfr_init2(re, (mpfr_prec_t)o->target_bits);
        mpfr_init2(im, (mpfr_prec_t)o->target_bits);
        mpfr_set_d(re, creal(val), MPFR_RNDN);
        mpfr_set_d(im, cimag(val), MPFR_RNDN);
        Expr* out = ni_mpfr_build(re, im, o->target_bits);
        mpfr_clears(re, im, (mpfr_ptr)0);
        return out;
    }
#endif
    if (st == CUB_NOCONV) ni_warn_goal(o, abserr, cabs(val));
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
        if (!conv) ni_warn_goal(o, err, cabs(val));
        return ni_from_complex_d(val);
    }
    return NULL;   /* no reducible axis — fall back to iterated quadrature */
}

/* Route a single {x, xmin, xmax} integral by the kinds of its bounds. */
#ifdef USE_MPFR
/* True when `e` mentions the symbol `var` (head positions included). */
static bool ni_mentions(const Expr* e, const char* var) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == var;
    if (e->type != EXPR_FUNCTION) return false;
    if (ni_mentions(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ni_mentions(e->data.function.args[i], var)) return true;
    return false;
}

/* The integrand with every exact numeric constant that is a factor of a
 * product or a term of a sum depending on `var` replaced by its value at the
 * working precision, once, instead of at every sample.  At a sample the
 * other operand of such a Times or Plus is an MPFR number, so contagion
 * would numericalize the constant at that same precision anyway; doing it
 * per sample re-ran it every time.  (99/100 E^(I t))^-10000, a contour
 * integral's Cauchy weight, evaluates to (100/99)^10000 E^(-10000 I t), and
 * the 66000-bit rational cost 2 ms a sample -- a gcd to scan it and a long
 * division to round it -- over the 32769 samples of a trapezoidal rule.
 * Only numeric, non-atomic or Rational / big constants are replaced (a small
 * integer coefficient is left as it is: exponents and arguments are never
 * touched, since only Times and Plus operands are), and nothing under a
 * holding head is entered.  Returns a new reference (the input itself when
 * nothing changes). */
static Expr* ni_prenumericalize(Expr* e, const char* var, NumericSpec spec) {
    if (!e || e->type != EXPR_FUNCTION || !ni_mentions(e, var)) return expr_copy(e);
    Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        uint32_t at = symtab_get_def(h->data.symbol.name)->attributes;
        if (at & (ATTR_HOLDALL | ATTR_HOLDFIRST | ATTR_HOLDREST | ATTR_HOLDALLCOMPLETE))
            return expr_copy(e);
    }
    bool arith = h->type == EXPR_SYMBOL &&
                 (h->data.symbol.name == SYM_Times || h->data.symbol.name == SYM_Plus);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    bool changed = false;
    for (size_t i = 0; i < n; i++) {
        Expr* x = e->data.function.args[i];
        Expr* y = NULL;
        if (arith && !ni_mentions(x, var) && x->type != EXPR_INTEGER &&
            x->type != EXPR_REAL && x->type != EXPR_MPFR &&
            expr_is_numeric_quantity(x)) {
            CommonInexactInfo xi = common_scan_inexact(x);
            if (!xi.has_inexact) {
                y = numericalize(x, spec);
                if (y && !(y->type == EXPR_MPFR || is_complex(y, NULL, NULL))) {
                    expr_free(y);
                    y = NULL;
                }
            }
        }
        if (!y) y = ni_prenumericalize(x, var, spec);
        if (y != x) changed = true;
        a[i] = y;
    }
    if (!changed) {
        for (size_t i = 0; i < n; i++) expr_free(a[i]);
        free(a);
        return expr_copy(e);
    }
    Expr* out = expr_new_function(expr_copy(h), a, n);
    free(a);
    return out;
}
#endif

static Expr* ni_dispatch_1d_impl(Expr* body, const char* var,
                                 Expr* amin, Expr* amax, const NiOpts* o);

static Expr* ni_dispatch_1d(Expr* body, const char* var,
                            Expr* amin, Expr* amax, const NiOpts* o) {
#ifdef USE_MPFR
    if (o->prec_mpfr && o->bits > 0) {
        NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = o->bits;
        Expr* pb = ni_prenumericalize(body, var, sp);
        NiOpts po = *o;
        if (pb != body) { po.pre_body = pb; po.exact_body = body; }
        Expr* r = ni_dispatch_1d_impl(pb, var, amin, amax, &po);
        expr_free(pb);
        return r;
    }
#endif
    return ni_dispatch_1d_impl(body, var, amin, amax, o);
}

static Expr* ni_dispatch_1d_impl(Expr* body, const char* var,
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

/* False, with NIntegrate::inumr, when the integrand has no numeric value at
 * two interior points of [a, b] (an unknown function, a broken interpolant). */
static bool ni_check_numeric(Expr* body, const char* var, double a, double b) {
    NiBind bind;
    ni_bind_snapshot(&bind, var);
    NiCtx ctx;
    ni_ctx_init(&ctx, body, &bind, numeric_machine_spec());
    const double fr[2] = { 0.3819660112501051, 0.6180339887498949 };
    bool numeric = false;
    for (int i = 0; i < 2 && !numeric; i++) {
        Expr* num = ni_eval_at(&ctx, expr_new_real(a + fr[i] * (b - a)));
        double _Complex v;
        if (num && (ni_to_complex(num, &v) || ni_beyond_double(num))) numeric = true;
        expr_free(num);
    }
    ni_bind_restore(&bind);
    if (!numeric)
        ni_warn("inumr", "the integrand has evaluated to non-numerical values "
                "for all sampling points in the region with boundaries {{%.17g, %.17g}}.",
                a, b);
    return numeric;
}

/* ------------------------------------------------------------------ *
 *  Interior points of non-analyticity at WorkingPrecision             *
 * ------------------------------------------------------------------ *
 *
 * Every arbitrary-precision rule here converges geometrically only for an
 * integrand analytic inside the interval; tanh-sinh and Clenshaw-Curtis absorb
 * singularities at the ENDPOINTS, not inside.  Log[Abs[x^2 - 2]] on [0, 4] has
 * a logarithmic singularity at Sqrt[2], so at 30 digits every rule stalled near
 * 5 correct digits and NIntegrate returned 3.6013587 for 3.6013302 (after 7 s
 * of levels); split at Sqrt[2], the same rules take 5 ms.  (The machine path
 * finds such a point by adaptive bisection, which the doubled cost of each
 * MPFR panel does not allow.)
 *
 * The points are read off the integrand: Log[u], Abs[u] and u^p (p not a
 * non-negative integer: poles, roots, fractional powers) are non-analytic
 * where u = 0, so for each such u that is a polynomial in the variable with
 * numeric coefficients, its real roots inside (a, b) become breakpoints.  The
 * roots are kept exact (Sqrt[2], Root objects), so the MPFR endpoint of each
 * piece is the singular point to the working precision (ni_endpoint_mpfr). */
#define NI_SING_MAX 32

#ifdef USE_MPFR


static void ni_sing_collect(const Expr* e, const char* var, Expr** us, int* nu) {
    if (!e || e->type != EXPR_FUNCTION || *nu >= NI_SING_MAX) return;
    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    const Expr* u = NULL;
    if (hd->type == EXPR_SYMBOL) {
        const char* h = hd->data.symbol.name;
        if ((h == SYM_Log || h == SYM_Abs) && n == 1) u = e->data.function.args[0];
        else if (h == SYM_Power && n == 2) {
            const Expr* p = e->data.function.args[1];
            bool nonneg_int = (p->type == EXPR_INTEGER && p->data.integer >= 0);
            double pd;
            if (!nonneg_int && ni_free_of_var(p, var) && ni_to_double_real((Expr*)p, &pd))
                u = e->data.function.args[0];
        }
    }
    if (u && !ni_free_of_var(u, var)) {
        bool dup = false;
        for (int i = 0; i < *nu && !dup; i++) dup = expr_eq(us[i], (Expr*)u);
        if (!dup) us[(*nu)++] = expr_copy((Expr*)u);
    }
    ni_sing_collect(hd, var, us, nu);
    for (size_t i = 0; i < n; i++) ni_sing_collect(e->data.function.args[i], var, us, nu);
}

/* Screen before the exact Solve: false only when the linear or quadratic u
 * certainly has no real root in [a, b] (the common 1/(c + x^2), Log[x - c]
 * with c outside), read off its coefficients from three machine values; any
 * doubt, or a higher degree, answers true. */
static bool ni_poly_may_vanish_on(Expr* u, const char* var, int deg, double a, double b) {
    if (deg > 2) return true;
    double v[3];
    const double xs[3] = { 0.0, 1.0, -1.0 };
    for (int i = 0; i < 3; i++) {
        Expr* r = ni_ev(SYM_ReplaceAll, (Expr*[]){ expr_copy(u),
                        expr_new_function(expr_new_symbol(SYM_Rule),
                                          (Expr*[]){ expr_new_symbol(var),
                                                     expr_new_integer((int64_t)xs[i]) }, 2) },
                        2, var);
        Expr* num = r ? numericalize(r, numeric_machine_spec()) : NULL;
        expr_free(r);
        bool ok = num && ni_to_double_real(num, &v[i]) && isfinite(v[i]);
        expr_free(num);
        if (!ok) return true;                 /* complex coefficients etc. */
    }
    double c0 = v[0], c1 = 0.5 * (v[1] - v[2]), c2 = 0.5 * (v[1] + v[2]) - v[0];
    double scale = fabs(c0) + fabs(c1) + fabs(c2);
    if (!(scale > 0.0)) return true;
    double tol = 1e-8 * (1.0 + fabs(a) + fabs(b));
    if (deg == 1 || fabs(c2) <= 1e-12 * scale) {
        if (fabs(c1) <= 1e-12 * scale) return true;
        double r = -c0 / c1;
        return r >= a - tol && r <= b + tol;
    }
    double disc = c1 * c1 - 4.0 * c2 * c0;
    if (disc < -1e-8 * (c1 * c1 + fabs(4.0 * c2 * c0))) return false;   /* complex pair */
    double sq = sqrt(disc > 0.0 ? disc : 0.0);
    double r1 = (-c1 - sq) / (2.0 * c2), r2 = (-c1 + sq) / (2.0 * c2);
    return (r1 >= a - tol && r1 <= b + tol) || (r2 >= a - tol && r2 <= b + tol);
}

/* Exact real roots in (a, b) of the polynomial factors of the integrand's
 * singular arguments, sorted, appended to roots[] (owned) with their doubles. */
static int ni_interior_singular_points(Expr* body, const char* var, double a, double b,
                                       Expr** roots, double* vals, int maxn) {
    Expr* us[NI_SING_MAX];
    int nu = 0, nr = 0;
    ni_sing_collect(body, var, us, &nu);
    for (int i = 0; i < nu; i++) {
        Expr* u = us[i];
        /* numeric coefficients only: the roots must be numbers */
        int deg = ni_numeric_poly_degree(u, var);
        if (deg < 1 || deg > NI_SING_MAX_DEGREE) continue;
        if (!ni_poly_may_vanish_on(u, var, deg, a, b)) continue;
        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                                     (Expr*[]){ expr_copy(u), expr_new_integer(0) }, 2);
        Expr* sol = ni_ev(SYM_Solve, (Expr*[]){ eq, expr_new_symbol(var),
                                                expr_new_symbol(SYM_Reals) }, 3, var);
        if (sol && sol->type == EXPR_FUNCTION && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol.name == SYM_List) {
            for (size_t j = 0; j < sol->data.function.arg_count && nr < maxn; j++) {
                Expr* rs = sol->data.function.args[j];
                if (rs->type != EXPR_FUNCTION || rs->data.function.arg_count != 1) continue;
                Expr* r = rs->data.function.args[0];
                if (r->type != EXPR_FUNCTION || r->data.function.arg_count != 2
                    || r->data.function.head->type != EXPR_SYMBOL
                    || r->data.function.head->data.symbol.name != SYM_Rule)
                    continue;
                Expr* val = r->data.function.args[1];
                /* the double of the correctly rounded value: ni_endpoint_mpfr
                 * matches the piece endpoint to its exact node by this double */
                NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = 160;
                Expr* num = numericalize(val, sp);
                double rv = NAN;
                bool in = false;
                if (num) {
                    mpfr_t re, im;
                    mpfr_inits2(160, re, im, (mpfr_ptr)0);
                    bool inexact;
                    if (get_approx_mpfr(num, re, im, &inexact) && mpfr_zero_p(im)
                        && mpfr_number_p(re)) {
                        rv = mpfr_get_d(re, MPFR_RNDN);
                        in = rv > a + 1e-12 * (b - a) && rv < b - 1e-12 * (b - a);
                    }
                    mpfr_clears(re, im, (mpfr_ptr)0);
                }
                expr_free(num);
                for (int k = 0; in && k < nr; k++)
                    if (fabs(vals[k] - rv) <= 1e-14 * (1.0 + fabs(rv))) in = false;
                if (!in) continue;
                /* insertion sort by value */
                int k = nr;
                while (k > 0 && vals[k - 1] > rv) {
                    vals[k] = vals[k - 1]; roots[k] = roots[k - 1]; k--;
                }
                vals[k] = rv; roots[k] = expr_copy(val);
                nr++;
            }
        }
        expr_free(sol);
    }
    for (int i = 0; i < nu; i++) expr_free(us[i]);
    return nr;
}
#endif /* USE_MPFR */


/* ------------------------------------------------------------------ *
 *  Phase substitution for a power-law phase                           *
 * ------------------------------------------------------------------ *
 *
 * An integrand A(x) K(c x^n), K = Sin or Cos, oscillates ever faster toward
 * an endpoint -- at x = 0 for n < 0 (Sin[1/x] on [0, 1]: infinitely many
 * oscillations; the adaptive rule reported 0.504166 "converged" for
 * 0.5040671) or at infinity for n > 1 (Sin[x^3] on [0, oo): the lobes shrink
 * like k^(-2/3), each a full panel at WorkingPrecision -> 40, which never
 * finished).  In the phase t = c x^n the oscillation is uniform:
 *     int A(x) K(c x^n) dx = int_{t0}^oo A((t/c)^(1/n)) (t/c)^(1/n)/(|n| t) K(t) dt,
 * a semi-infinite Fourier integral with frequency 1, which the half-line
 * Fourier rules (Ooura-Mori at WorkingPrecision, between the zeros at machine
 * precision) resolve.  t0 = 0 (n > 0 on [0, oo)) or c b^n (n < 0 on [0, b]);
 * for t0 > 0, K(t0 + s) is expanded into Cos[s] and Sin[s] so each part keeps
 * the aligned phase those rules need.  Returns the evaluated sum, or NULL
 * when the integrand or range is not of this form. */

/* evaluate f[args...] (args consumed) with `var` detached from any binding */
static Expr* ni_ev(const char* head, Expr** args, size_t n, const char* var) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    int vtok = mth_bind_push_name(var, 0);
    eval_clock_bump();
    Expr* r = eval_and_free(call);
    mth_bind_pop(vtok);
    eval_clock_bump();
    return r;
}
static Expr* ni_ev2(const char* head, Expr* a, Expr* b, const char* var) {
    Expr* args[2] = { a, b };
    return ni_ev(head, args, 2, var);
}
static Expr* ni_ev1(const char* head, Expr* a, const char* var) {
    Expr* args[1] = { a };
    return ni_ev(head, args, 1, var);
}

static Expr* ni_try_phase_substitution(Expr* body, const char* var, Expr* amin,
                                       Expr* amax, Expr* const* opts, size_t nopts) {
    LevinKernel k; Expr* g = NULL; Expr* amp = NULL;
    Expr* fb = expr_copy(body);
    if (levin_is_head(fb, SYM_Times)) {
        fb = expr_unshare(fb);
        eval_flatten_args(fb, SYM_Times);
    }
    bool ok = levin_detect_kernel(fb, var, &k, &g, &amp);
    expr_free(fb);
    if (!ok) return NULL;
    if (k != LEVIN_KERNEL_SIN && k != LEVIN_KERNEL_COS) { expr_free(g); expr_free(amp); return NULL; }
    Expr* out = NULL;
    /* n = x g'/g must be a number, c = g/x^n a positive number */
    Expr* gp = levin_phase_derivative(g, var);
    Expr* nE = gp ? ni_ev1(SYM_Together,
                           ni_ev2(SYM_Times, expr_new_symbol(var),
                                  ni_ev2(SYM_Times, gp, ni_ev2(SYM_Power, expr_copy(g),
                                                              expr_new_integer(-1), var), var),
                                  var), var) : NULL;
    double nd = 0.0, cd = 0.0, ad = 0.0, bd = 0.0;
    bool b_inf = false;
    Expr* cE = NULL;
    if (!nE || levin_occurs(nE, var) || !ni_to_double_real(nE, &nd) || nd == 0.0 || nd == 1.0)
        goto done;
    cE = ni_ev1(SYM_Together, ni_ev2(SYM_Times, expr_copy(g),
                ni_ev2(SYM_Power, expr_new_symbol(var),
                       ni_ev2(SYM_Times, expr_new_integer(-1), expr_copy(nE), var), var), var), var);
    if (!cE || levin_occurs(cE, var)) goto done;
    {
        Expr* cN = numericalize(cE, numeric_machine_spec());
        bool okc = cN && ni_to_double_real(cN, &cd);
        expr_free(cN);
        if (!okc || !(cd != 0.0)) goto done;
    }
    {
        Expr* aN = ni_num_endpoint(amin);
        bool oka = aN && ni_to_double_real(aN, &ad);
        expr_free(aN);
        double _Complex dir;
        if (ni_infinity_dir(amax, &dir) && creal(dir) > 0 && cimag(dir) == 0) b_inf = true;
        else {
            Expr* bN = ni_num_endpoint(amax);
            bool okb = bN && ni_to_double_real(bN, &bd);
            expr_free(bN);
            if (!okb) goto done;
        }
        if (!oka || ad != 0.0) goto done;
        if (nd > 1.0 && !b_inf) goto done;            /* phase grows toward oo */
        if (nd < 0.0 && (b_inf || !(bd > 0.0))) goto done;   /* toward x = 0 */
        if (nd > 0.0 && nd < 1.0) goto done;
    }
    {
        /* sign of c: K(-t) = -Sin t / Cos t */
        double ksign = 1.0;
        if (cd < 0.0) {
            cE = ni_ev2(SYM_Times, expr_new_integer(-1), cE, var);
            if (k == LEVIN_KERNEL_SIN) ksign = -1.0;
        }
        /* x(t) = (t/c)^(1/n);  A(x(t)) x(t)/(|n| t) */
        Expr* xt = ni_ev2(SYM_Power, ni_ev2(SYM_Times, expr_new_symbol(var),
                          ni_ev2(SYM_Power, expr_copy(cE), expr_new_integer(-1), var), var),
                          ni_ev2(SYM_Power, expr_copy(nE), expr_new_integer(-1), var), var);
        Expr* rule = ni_ev2(SYM_Rule, expr_new_symbol(var), expr_copy(xt), var);
        Expr* ampt = ni_ev2(SYM_ReplaceAll, expr_copy(amp), rule, var);
        Expr* jac = ni_ev2(SYM_Times, xt,
                           ni_ev2(SYM_Power, ni_ev2(SYM_Times,
                                  ni_ev1(SYM_Abs, expr_copy(nE), var), expr_new_symbol(var), var),
                                  expr_new_integer(-1), var), var);
        Expr* B = ni_ev2(SYM_Times, ampt, jac, var);       /* B(t) */
        /* t0: c b^n for n < 0.  For n > 1 the substitution would create an
         * algebraic singularity t^(1/n - 1) at t = 0 where the original is
         * smooth (Sin[x^3] ~ x^3), so the first lobe x in [0, (Pi/c)^(1/n)] is
         * integrated as it stands and the phase variable starts at t0 = Pi. */
        Expr* head_part = NULL;
        Expr* t0 = (nd > 0.0) ? expr_new_symbol(SYM_Pi)
                 : ni_ev2(SYM_Times, expr_copy(cE), ni_ev2(SYM_Power, expr_copy(amax),
                                                           expr_copy(nE), var), var);
        if (nd > 0.0) {
            Expr* x1 = ni_ev2(SYM_Power, ni_ev2(SYM_Times, expr_new_symbol(SYM_Pi),
                              ni_ev2(SYM_Power, expr_copy(cE), expr_new_integer(-1), var), var),
                              ni_ev2(SYM_Power, expr_copy(nE), expr_new_integer(-1), var), var);
            Expr** args = malloc((2 + nopts) * sizeof(Expr*));
            args[0] = expr_copy(body);
            Expr* sp[3] = { expr_new_symbol(var), expr_new_integer(0), x1 };
            args[1] = expr_new_function(expr_new_symbol(SYM_List), sp, 3);
            for (size_t j = 0; j < nopts; j++) args[2 + j] = expr_copy(opts[j]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_NIntegrate), args, 2 + nopts);
            free(args);
            head_part = eval_and_free(call);
            if (!head_part || levin_is_head(head_part, SYM_NIntegrate)) {
                expr_free(head_part); expr_free(t0); expr_free(B);
                goto done;
            }
        }
        double t0d = 0.0;
        { Expr* tN = numericalize(expr_copy(t0), numeric_machine_spec());
          ni_to_double_real(tN, &t0d); expr_free(tN); }
        /* parts: (coefficient, kernel head) over s in [0, oo), B(t0 + s) */
        Expr* Bs = B;
        if (t0d != 0.0) {
            Expr* r2 = ni_ev2(SYM_Rule, expr_new_symbol(var),
                              ni_ev2(SYM_Plus, expr_copy(t0), expr_new_symbol(var), var), var);
            Bs = ni_ev2(SYM_ReplaceAll, B, r2, var);
        }
        const char* kh = (k == LEVIN_KERNEL_SIN) ? SYM_Sin : SYM_Cos;
        Expr* coef[2]; const char* kern[2]; int nparts;
        if (t0d == 0.0) {
            coef[0] = expr_new_integer(1); kern[0] = kh; nparts = 1;
        } else if (k == LEVIN_KERNEL_SIN) {     /* Sin[t0+s] = Sin t0 Cos s + Cos t0 Sin s */
            coef[0] = ni_ev1(SYM_Sin, expr_copy(t0), var); kern[0] = SYM_Cos;
            coef[1] = ni_ev1(SYM_Cos, expr_copy(t0), var); kern[1] = SYM_Sin; nparts = 2;
        } else {                                /* Cos[t0+s] = Cos t0 Cos s - Sin t0 Sin s */
            coef[0] = ni_ev1(SYM_Cos, expr_copy(t0), var); kern[0] = SYM_Cos;
            coef[1] = ni_ev2(SYM_Times, expr_new_integer(-1),
                             ni_ev1(SYM_Sin, expr_copy(t0), var), var); kern[1] = SYM_Sin; nparts = 2;
        }
        Expr* total = NULL;
        bool good = true;
        for (int i = 0; i < nparts; i++) {
            /* an exactly vanishing coefficient (Sin[t0] for the lobe start
             * t0 = Pi) needs no integral: half the work for every n > 1 */
            if (coef[i]->type == EXPR_INTEGER && coef[i]->data.integer == 0) continue;
            Expr* integrand = ni_ev2(SYM_Times, expr_copy(Bs),
                                     ni_ev1(kern[i], expr_new_symbol(var), var), var);
            Expr** args = malloc((2 + nopts) * sizeof(Expr*));
            args[0] = integrand;
            Expr* sp[3] = { expr_new_symbol(var), expr_new_integer(0), expr_new_symbol(SYM_Infinity) };
            args[1] = expr_new_function(expr_new_symbol(SYM_List), sp, 3);
            for (size_t j = 0; j < nopts; j++) args[2 + j] = expr_copy(opts[j]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_NIntegrate), args, 2 + nopts);
            free(args);
            Expr* r = eval_and_free(call);
            if (!r || levin_is_head(r, SYM_NIntegrate)) { expr_free(r); good = false; }
            else {
                Expr* term = ni_ev2(SYM_Times, coef[i], r, var);
                coef[i] = NULL;
                total = total ? ni_ev2(SYM_Plus, total, term, var) : term;
            }
            if (!good) break;
        }
        for (int i = 0; i < nparts; i++) expr_free(coef[i]);
        expr_free(Bs);
        expr_free(t0);
        if (good && total) {
            out = (ksign < 0) ? ni_ev2(SYM_Times, expr_new_integer(-1), total, var) : total;
            total = NULL;
            if (head_part) { out = ni_ev2(SYM_Plus, head_part, out, var); head_part = NULL; }
        }
        expr_free(total);
        expr_free(head_part);
    }
done:
    expr_free(cE);
    expr_free(nE);
    expr_free(g);
    expr_free(amp);
    return out;
}

/* Integrate res = NIntegrate[f, specs..., opts...] with the parsed options:
 * the dispatch onto Monte-Carlo, cubature, iterated, contour and 1-D paths. */
static Expr* ni_integrate(Expr* res, size_t argc, size_t pos_end, const NiOpts* op) {
    NiOpts o = *op;
    Expr* body = res->data.function.args[0];

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

    /* At WorkingPrecision the cubature runs with MPFR-evaluated samples but
     * double accumulation, so it is taken only for a precision goal a double
     * sum still holds; a higher goal keeps the iterated MPFR rules. */
    if (o.method == NI_AUTO && nspecs >= 2 && !o.exclusions
        && (!o.prec_mpfr || ni_cub_mpfr_goal_ok(&o))) {
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
    o.spec_nodes = &spec0->data.function.args[1];
    o.spec_nnodes = nnodes;
    const char* var = spec0->data.function.args[0]->data.symbol.name;

    Expr* out = NULL;
    if (nnodes == 2) {
        /* {x, xmin, xmax}: real finite / infinite / complex segment / ray. */
        Expr* amin = ni_num_endpoint(spec0->data.function.args[1]);
        Expr* amax = ni_num_endpoint(spec0->data.function.args[2]);
        double a, b;
        bool real_finite = ni_to_double_real(amin, &a) && ni_to_double_real(amax, &b);
        double pts[64]; int npts = 0;
        /* An integrand with no numeric value at interior points has no
         * integral to approximate: unevaluated with NIntegrate::inumr (the
         * machine cascade's last rules used to sum nothing and report 0.). */
        if (real_finite && a != b && !ni_check_numeric(eff_body, var, a, b)) {
            expr_free(amin);
            expr_free(amax);
            expr_free(inner);
            return NULL;
        }
        /* a power-law phase: integrate in the phase variable */
        if (o.method == NI_AUTO && nspecs == 1 && !o.exclusions) {
            Expr* ps = ni_try_phase_substitution(eff_body, var, spec0->data.function.args[1],
                                                 spec0->data.function.args[2],
                                                 &res->data.function.args[pos_end],
                                                 argc - pos_end);
            if (ps) {
                expr_free(amin);
                expr_free(amax);
                expr_free(inner);
                return ps;
            }
        }
#ifdef USE_MPFR
        /* interior points of non-analyticity, split exactly at WorkingPrecision */
        if (real_finite && o.prec_mpfr && o.method == NI_AUTO && nspecs == 1
            && !o.exclusions && a != b) {
            Expr* roots[NI_SING_MAX];
            double vals[NI_SING_MAX];
            double lo = a < b ? a : b, hi = a < b ? b : a;
            int nr = ni_interior_singular_points(eff_body, var, lo, hi, roots, vals,
                                                 NI_SING_MAX);
            if (nr > 0) {
                double _Complex nodes[NI_SING_MAX + 2];
                Expr* enodes[NI_SING_MAX + 2];
                int nn = 0;
                bool rev = a > b;
                nodes[nn] = a; enodes[nn++] = spec0->data.function.args[1];
                for (int i = 0; i < nr; i++) {
                    int k = rev ? nr - 1 - i : i;
                    nodes[nn] = vals[k]; enodes[nn++] = roots[k];
                }
                nodes[nn] = b; enodes[nn++] = spec0->data.function.args[2];
                NiOpts os = o;
                os.spec_nodes = enodes;
                os.spec_nnodes = (size_t)nn;
                out = ni_run_contour(eff_body, var, nodes, nn, &os);
                for (int i = 0; i < nr; i++) expr_free(roots[i]);
                expr_free(amin);
                expr_free(amax);
                expr_free(inner);
                return out;
            }
        }
#endif
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

/* Symbolic processing must not make the integrand costlier to sample than the
 * integrand as written.  Evaluating it with the variables unbound can expand
 * a compact closed form into a huge one: LaguerreL[1000, r E^(I phi)] becomes
 * its degree-1000 polynomial in E^(I phi), which costs a thousand terms per
 * sample and, worse, cancels catastrophically (|L_1000| ~ 1e1018 against
 * terms of 1e2500 at r = 3940), where the held form is evaluated pointwise by
 * the stable three-term recurrence.  So when the evaluated integrand is many
 * times larger than the held one, the held one is kept -- provided it has a
 * numeric value with the variables bound to interior points of their ranges
 * (it would not if, say, it were a held D[f[x], x], which only means something
 * symbolically).  The probe is silent: its messages, if any, are muted and
 * kept out of Check[]. */
#define NI_SWELL_RATIO  16
#define NI_SWELL_MIN    256
static bool ni_symbolic_swell(Expr* held, Expr* evald) {
    int64_t lh = leaf_count_internal(held, true);
    int64_t le = leaf_count_internal(evald, true);
    return le > NI_SWELL_MIN && le > NI_SWELL_RATIO * lh;
}

static bool ni_held_body_numeric(Expr* body, Expr** specs, size_t nspecs) {
    if (nspecs == 0 || nspecs > 16) return false;
    NiBind binds[16]; size_t nb = 0;
    bool ok = true;
    unsigned long fired = mth_msg_fired_count();
    mth_msg_suppress_push();
    for (size_t i = 0; i < nspecs && ok; i++) {
        Expr* sp = specs[i];
        if (sp->type != EXPR_FUNCTION || sp->data.function.arg_count < 3 ||
            sp->data.function.args[0]->type != EXPR_SYMBOL) { ok = false; break; }
        Expr* lo = ni_num_endpoint(sp->data.function.args[1]);
        Expr* hi = ni_num_endpoint(sp->data.function.args[2]);
        NiBound bl, bh;
        bl.kind = bh.kind = NI_B_BAD;
        if (lo) ni_classify(lo, &bl);
        if (hi) ni_classify(hi, &bh);
        expr_free(lo); expr_free(hi);
        /* an interior point: the golden section of a finite range, one unit
         * inside a half line, near the origin on the whole line */
        double _Complex pt;
        bool lf = (bl.kind == NI_B_FIN || bl.kind == NI_B_CPLX);
        bool hf = (bh.kind == NI_B_FIN || bh.kind == NI_B_CPLX);
        if (lf && hf)                                   pt = bl.z + 0.3819660112501051 * (bh.z - bl.z);
        else if (lf && bh.kind == NI_B_PINF)            pt = bl.z + 0.6180339887498949;
        else if (hf && bl.kind == NI_B_NINF)            pt = bh.z - 0.6180339887498949;
        else if (bl.kind == NI_B_NINF && bh.kind == NI_B_PINF) pt = 0.3819660112501051;
        else { ok = false; break; }
        ni_bind_snapshot(&binds[nb], sp->data.function.args[0]->data.symbol.name);
        nb++;
        ni_bind_set(&binds[nb - 1], ni_from_complex_d(pt));
    }
    if (ok) {
        /* every variable is bound; the context binding only names one of them,
         * so evaluate the body directly at the bound values */
        eval_clock_bump();
        arith_warnings_mute_push();
        Expr* raw = eval_and_free(expr_copy(body));
        arith_warnings_mute_pop();
        Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
        if (raw) expr_free(raw);
        double _Complex v;
        ok = num && (ni_to_complex(num, &v) || ni_beyond_double(num));
        if (num) expr_free(num);
    }
    /* the binding frames are a stack: pop innermost first */
    while (nb > 0) ni_bind_restore(&binds[--nb]);
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired);
    return ok;
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

    /* {x, y} \[Element] reg is accepted only for the finite element path
     * below (a mesh function integrated over its own region) */
    bool region_spec = pos_end == 2 && fem_is_region_spec(res->data.function.args[1]);
    for (size_t i = 1; i < pos_end && !region_spec; i++)
        if (!ni_is_spec(res->data.function.args[i])) return NULL;

    /* Localise the integration variables and evaluate the integrand
     * symbolically once, then integrate the result, which is the reference's
     * documented order. NIntegrate is HoldAll, so without this a held wrapper
     * such as f /. b -> 1 reached the arbitrary-precision sampler
     * unevaluated. NIntegrate[y^(b-1)/(1+y)^2 /. b -> 1, {y, 0, Infinity},
     * WorkingPrecision -> 20] returned 0. When the evaluated integrand
     * differs, the call is re-dispatched with it. The second entry finds it
     * unchanged, so this does not recurse. */
    {
        Expr* vars[16]; size_t nv = 0;
        for (size_t i = 1; i < pos_end && nv < 16; i++) {
            Expr* sp = res->data.function.args[i];
            if (sp->type == EXPR_FUNCTION && sp->data.function.arg_count >= 1 &&
                sp->data.function.args[0]->type == EXPR_SYMBOL)
                vars[nv++] = expr_copy(sp->data.function.args[0]);
            else if (region_spec)
                for (size_t k = 0; k < sp->data.function.args[0]->data.function.arg_count && nv < 16; k++)
                    if (sp->data.function.args[0]->data.function.args[k]->type == EXPR_SYMBOL)
                        vars[nv++] = expr_copy(sp->data.function.args[0]->data.function.args[k]);
        }
        /* A Block frame without Block's re-evaluation on exit: the integrand
         * is evaluated with the variables unbound, and the result is used as
         * is. (Block[{y}, y^2] would come back as 25 under an outer y = 5.)
         * The frame is drained by a TimeConstrained unwind. */
        void* frame = mth_block_enter(vars, nv);
        Expr* eb = evaluate(body);
        mth_block_leave(frame);
        for (size_t i = 0; i < nv; i++) expr_free(vars[i]);
        /* a function on a finite element mesh integrated over its own
         * region: element by element (fem.c) */
        if (eb && pos_end <= 3) {
            Expr* r = fem_nintegrate(eb, res->data.function.args + 1, pos_end - 1);
            if (r) { expr_free(eb); return r; }
        }
        if (region_spec) { expr_free(eb); return NULL; }
        /* an expansion that swells the integrand is not used when the held
         * integrand samples numerically (ni_symbolic_swell) */
        if (eb && ni_symbolic_swell(body, eb) &&
            ni_held_body_numeric(body, res->data.function.args + 1, pos_end - 1)) {
            expr_free(eb);
            eb = NULL;
        }
        if (eb && !expr_eq(eb, body)) {
            Expr** args = malloc(argc * sizeof(Expr*));
            args[0] = eb;
            for (size_t i = 1; i < argc; i++) args[i] = expr_copy(res->data.function.args[i]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_NIntegrate), args, argc);
            free(args);
            Expr* r = evaluate(call);
            expr_free(call);
            /* still unevaluated: hand back NULL so the original call stays */
            if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
                r->data.function.head->data.symbol.name == SYM_NIntegrate) {
                expr_free(r);
                return NULL;
            }
            return r;
        }
        if (eb) expr_free(eb);
    }

    NiOpts o;
    o.method = NI_AUTO;
    o.method_name = NULL;
    o.prec_mpfr = false; o.bits = 0; o.target_bits = 0;
    /* Default AccuracyGoal -> MachinePrecision (a ~15.95-digit absolute-error
     * floor); PrecisionGoal stays Automatic (-> WorkingPrecision/2).  The floor
     * only relaxes the combined tolerance 10^-a + |x| 10^-p, since a finite a
     * adds a positive term where the old Automatic contributed none. */
    o.acc_goal = NUMERIC_MACHINE_PRECISION_DIGITS; o.prec_goal = -1.0;
    o.acc_goal_user = false;
    o.max_recursion = -1; o.min_recursion = 0; o.max_points = -1;
    o.exclusions = NULL;
    o.rule_type = 0; o.romberg = true; o.nc_points = 3;
    o.spec_nodes = NULL; o.spec_nnodes = 0;
    o.pre_body = NULL; o.exact_body = NULL;
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

    /* Machine pass; when a sample left the double range, the pass is replaced
     * by an MPFR pass at machine-equivalent precision (see g_ni_beyond). */
    bool saved_beyond = g_ni_beyond, saved_escalated = g_ni_escalated;
    long saved_sbits = g_ni_sample_bits;
    long saved_ebits = g_ni_endpoint_bits;
    double saved_gnoise = g_ni_guard_noise_l2;
    g_ni_endpoint_bits = 0;
    g_ni_beyond = false;
    g_ni_escalated = false;
    g_ni_sample_bits = 0;
    g_ni_guard_noise_l2 = -INFINITY;
    Expr* out = ni_integrate(res, argc, pos_end, &o);
#ifdef USE_MPFR
    if (!o.prec_mpfr && g_ni_beyond) {
        NiOpts m = o;
        m.prec_mpfr = true;
        m.target_bits = numeric_digits_to_bits(NUMERIC_MACHINE_PRECISION_DIGITS);
        m.bits = m.target_bits + 64;
        /* PrecisionGoal 12 unless given: the result carries 53 bits, and the
         * nested rules make the last few digits nearly free */
        if (!(m.prec_goal > 0.0)) m.prec_goal = 12.0;
        g_ni_beyond = false;
        g_ni_escalated = true;
        Expr* mo = ni_integrate(res, argc, pos_end, &m);
        /* the machine pass dropped its out-of-range samples: never keep it */
        expr_free(out);
        out = mo;
    }
#endif
    /* WorkingPrecision asked of an integrand whose values carry fewer digits:
     * the rules stop at the noise floor of those values (they no longer spin
     * through every level), and the shortfall is reported. */
    if (out && o.prec_mpfr && g_ni_sample_bits > 0
        && g_ni_sample_bits < o.target_bits && !arith_warnings_muted())
        ni_warn("precw", "the precision of the integrand values (%.3g digits) is "
                "less than the WorkingPrecision (%.3g digits).",
                numeric_bits_to_digits(g_ni_sample_bits),
                numeric_bits_to_digits(o.target_bits));
#ifdef USE_MPFR
    /* An endpoint cancellation the guard could not resolve (see
     * g_ni_guard_noise_l2): the result carries only the digits that error
     * leaves, never the WorkingPrecision it did not reach. */
    if (out && o.prec_mpfr && isfinite(g_ni_guard_noise_l2))
        out = ni_cap_precision(out, g_ni_guard_noise_l2, o.target_bits, false);
    /* An exact endpoint that numericalised short of the working precision:
     * the result carries the endpoint's relative precision at most. */
    if (out && o.prec_mpfr && g_ni_endpoint_bits > 0
        && g_ni_endpoint_bits < o.target_bits) {
        Expr *cre = NULL, *cim = NULL;
        Expr* mag = is_complex(out, &cre, &cim) ? cre : out;
        if (mag->type == EXPR_MPFR && mpfr_number_p(mag->data.mpfr)
            && !mpfr_zero_p(mag->data.mpfr))
            out = ni_cap_precision(out,
                                   de_mpfr_log2_abs(mag->data.mpfr)
                                       - (double)g_ni_endpoint_bits,
                                   o.target_bits, true);
    }
#endif
    g_ni_beyond = saved_beyond;
    g_ni_escalated = saved_escalated;
    g_ni_sample_bits = saved_sbits;
    g_ni_endpoint_bits = saved_ebits;
    g_ni_guard_noise_l2 = saved_gnoise;
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
        "(\"GlobalAdaptive\", \"DoubleExponential\", \"Trapezoidal\" (plain "
        "doubling, for periodic integrands), \"TrapezoidalRule\", "
        "\"RiemannRule\", \"NewtonCotesRule\", \"LevinRule\", "
        "\"OscillatorySingularity\", \"MonteCarlo\", \"PrincipalValue\", ...), "
        "WorkingPrecision, PrecisionGoal, AccuracyGoal, MaxRecursion, MaxPoints, "
        "Exclusions. Fixed rules take Method sub-options \"Type\" "
        "(\"Left\"/\"Right\"/\"Midpoint\"), \"RombergQuadrature\" and \"Points\". "
        "A machine-precision integrand that leaves the double range is "
        "re-integrated at machine-equivalent arbitrary precision.");
}
