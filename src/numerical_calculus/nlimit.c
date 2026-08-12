/*
 * nlimit.c — NLimit[expr, z -> z0, opts]
 *
 * Numerically estimates  lim_{z -> z0} expr  by sampling expr on a geometric
 * sequence of points approaching z0 and accelerating the resulting sequence.
 *
 *   Sample points
 *   -------------
 *   Finite z0:    z_k = z0 - d * Scale * 2^-k          (k = 0 .. Terms-1)
 *                 d is the (unit) Direction vector; the points sit on the -d
 *                 side of z0, i.e. one moves *along* d to reach z0.  The
 *                 default Direction -> Automatic == -1 places the points at
 *                 z0 + Scale*2^-k, approaching "from larger values".
 *   Infinite z0:  z_k = u * Scale * 2^k                (k = 0 .. Terms-1)
 *                 u is the unit ray of the infinite limit point (the direction
 *                 of Infinity / I Infinity / DirectedInfinity[..]); the points
 *                 march outward to infinity along that ray.
 *
 *   Acceleration
 *   ------------
 *   Method -> Automatic (default): run BOTH Richardson and Wynn's epsilon (at
 *     every admissible degree) and return the estimate whose internal
 *     convergence residual is smallest.  Richardson models an integer-power
 *     (analytic) error tail; its fixed 2^j-1 denominators cannot annihilate a
 *     geometric or fractional-power (branch-point) tail — e.g. the sqrt(step)
 *     imaginary part of 2 ArcTan[Sqrt[(1+x)/(1-x)]] approaching x -> 1 from the
 *     larger-values side.  Wynn's epsilon captures exactly those tails.
 *     Selecting by best self-consistency (smallest step) picks the right tool
 *     per problem, so smooth limits keep Richardson's accuracy while
 *     algebraic/branch approaches gain Wynn's.
 *
 *   Method -> EulerSum: Richardson / Romberg extrapolation of the sequence S_k,
 *     treated as a function of the geometric step.  Following the same
 *     convention as ND's "EulerSum" (see nderiv.c) the tableau uses the
 *     all-powers denominator 2^j - 1:
 *         T(i,0) = S_i,
 *         T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
 *         result = T(Terms-1, Terms-1).
 *
 *   Method -> SequenceLimit: Wynn's epsilon algorithm (iterated Shanks /
 *     Aitken).  With ε_{-1}=0, ε_0 = S_n,
 *         ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1/(ε_k^{(n+1)} - ε_k^{(n)});
 *     the even columns ε_{2d} are the limit estimates.  WynnDegree -> d selects
 *     column 2d and requires at least 2(d+1) terms for a convergence check.
 *
 *   Robustness
 *   ----------
 *   Two independent gates, both returning the form unevaluated when they fire.
 *
 *   Oscillatory divergence (NLimit::osc).  An extrapolator returns a number for
 *   any input, including a sequence with no limit, so before trusting one we
 *   check that the oscillation envelope actually decays: samples whose
 *   increments keep reversing direction get a wide auxiliary ladder, and a
 *   non-decaying envelope is refused.  See the NL_OSC_* block below.
 *
 *   Noise (NLimit::noise).  The last two extrapolates are compared; if they
 *   fail to agree to a loose tolerance relative to the result magnitude (or the
 *   result is non-finite), the form is returned unevaluated.  As with
 *   Mathematica, spurious tiny residuals are *not* recognised as zero — apply
 *   Chop when needed.
 *
 * Options: Method (Automatic | EulerSum | SequenceLimit), WorkingPrecision (MachinePrecision
 * | digits -> MPFR), Direction (Automatic == -1, or a complex approach vector),
 * Scale (initial step / distance, default 1), Terms (default 13), WynnDegree
 * (default 1).  The default of 13 (not the historical 7) is deliberate: the
 * accuracy of a branch-point / fractional-power approach is limited by the
 * *depth of the extrapolation tableau*, not by the arithmetic precision, so a
 * short sequence starves the extrapolators regardless of WorkingPrecision.
 * Thirteen samples resolve such tails to ~12 machine digits while leaving
 * smooth limits at their existing roundoff floor; the extra six evaluations are
 * sub-millisecond, and the Automatic best-of selector discards any deep tableau
 * that only adds roundoff, so the larger default never degrades an easy case.
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.  Every temporary OwnValue
 * created for the sampler is removed on all return paths.
 */

#include "nlimit.h"

#include <complex.h>
#include <gmp.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#  include "numeric_complex.h"   /* numeric_mpfr_make_complex, mpfr_complex_div */
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "common.h"       /* common_method_alias -- the NLimit`m heads */
#include "eval.h"
#include "nc_accuracy.h"  /* shared AccuracyGoal/PrecisionGoal handling */
#include "numeric.h"
#include "seqaccel.h"     /* shared Richardson + Wynn-epsilon kernels */
#include "sym_names.h"
#include "symtab.h"

/* Hard cap on Terms: 2^(Terms-1) must stay an exact double for the machine
 * Richardson denominator, and any larger tableau is numerical nonsense. */
#define NL_MAX_TERMS 50

/* Adaptive-refinement guards.  Growing the sample count past the point where
 * roundoff/cancellation dominates does not improve a machine-precision limit --
 * it starts producing extrapolants with a spuriously small residual but a
 * garbage value.  A refinement is therefore accepted only if it moves the
 * estimate by no more than NL_REFINE_JUMP times the current error estimate (a
 * genuine improvement cannot jump further than that), and growth stops after
 * NL_REFINE_STALL consecutive samples that yield no accepted improvement. */
#define NL_REFINE_JUMP  8.0
#define NL_REFINE_STALL 3

/* Convergence / noise gate.  Measured against the *sample scale* (the largest
 * |S_k|), not the result, so that genuine near-zero limits are not mistaken for
 * noise.  The gate is deliberately lenient: its job is to reject clearly
 * divergent or non-convergent sequences (e.g. 1/x as x->0, whose extrapolant
 * escapes far beyond the sampled data), not to certify accuracy.  A bounded but
 * only roughly-resolved limit (e.g. Sin[x]/x as x->Infinity, ~0) is accepted —
 * matching NLimit, which returns such values rather than refusing.
 *
 *   - reject if the result is non-finite;
 *   - reject if |result| exceeds NL_DIVERGE_FACTOR * (max |S_k|)  (divergence);
 *   - reject if the last two extrapolates differ by more than
 *     NL_NOISE_RELTOL * (max |S_k|)                               (no settling). */
#define NL_DIVERGE_FACTOR 4.0
#define NL_NOISE_RELTOL   0.5
#define NL_NOISE_ABSTOL   1e-6

static bool nl_accept(double result_mag, double step, double maxsample) {
    if (!isfinite(result_mag) || !isfinite(step)) return false;
    if (maxsample < NL_NOISE_ABSTOL) maxsample = NL_NOISE_ABSTOL;
    if (result_mag > NL_DIVERGE_FACTOR * maxsample) return false;
    if (step > NL_NOISE_RELTOL * maxsample + NL_NOISE_ABSTOL) return false;
    return true;
}

/* Oscillatory-divergence gate.
 * ---------------------------
 * nl_accept only asks whether the *extrapolant* is plausible next to the sample
 * scale.  It cannot ask whether the sample sequence converges at all, and the
 * extrapolators happily return a number for a sequence that has no limit: e.g.
 * (Cos[x^2]/x^2 - Cos[(x+1)^2]/(x+1)^2) x^3 at x -> Infinity, an amplitude-
 * modulated oscillation of growing envelope, used to yield -5.28256.
 *
 * A limit exists only if the oscillation envelope decays, so that is what we
 * measure.  Two stages:
 *
 *   1. Screen (free).  Count direction reversals of the increments
 *      D_k = S_k - S_{k-1}; a reversal is Re(D_k conj(D_{k-1})) < 0.  Monotone
 *      and smoothly-converging samples score 0 and skip stage 2 entirely, so
 *      every ordinary limit costs exactly what it did before.
 *
 *   2. Envelope (only when the screen fires).  Over the default sampling window
 *      of Terms=13 octaves a decaying and a non-decaying envelope are still not
 *      distinguishable, so the diagnosis gets its own, much wider ladder:
 *      NL_OSC_OCTAVES octaves sampled twice, at Scale*2^k and Scale*phi*2^k.
 *      Because 1/2 < phi < 1 the merged sequence is a deterministic alternation
 *      (no sort needed), and the phi offset breaks the power-of-two aliasing
 *      that hides e.g. Sin[Pi x].  With
 *          env = (max |f| over the half nearest the limit point)
 *              / (max |f| over the half furthest from it)
 *      the measured separation over the regression battery is env <= 0.22 for
 *      every limit that exists (worst: Sin[x]/x^(1/4) -> 0) against env >= 0.92
 *      for every non-limit (worst: Sin[x] + Sin[Pi x]).  NL_OSC_ENV_MAX = 0.6
 *      sits between them with ~2.7x and ~1.5x margin.
 *
 * The gate abstains -- never rejects -- when the ladder yields too few finite
 * points, so overflow or underflow far from the limit point cannot manufacture
 * a refusal.  Known blind spots: envelopes decaying slower than about x^(-1/10)
 * across the window are indistinguishable from non-decaying ones (refusing is
 * the safer error), and a sequence sampled exactly on its own zeros (Sin[Pi x])
 * shows no reversal for the screen to see. */
#define NL_OSC_OCTAVES        20
#define NL_OSC_ENV_MAX        0.6
#define NL_OSC_PHI            0.6180339887498949  /* (Sqrt[5]-1)/2, in (1/2, 1) */
#define NL_OSC_MIN_REVERSALS  1
/* Both ladders over at least half the window: a truncated ladder must still
 * span enough octaves for the envelope ratio to mean anything. */
#define NL_OSC_MIN_POINTS     (NL_OSC_OCTAVES)

#define NL_OSC_MESSAGE                                                        \
    "The sampled values oscillate with a non-decaying envelope; no limiting "  \
    "value exists as the limit point is approached."

static void nl_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "NLimit::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf <-> C value helpers (mirrors nderiv.c / nresidue.c)   *
 * ------------------------------------------------------------------ */

static bool nl_to_double_real(Expr* e, double* out) {
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

static bool nl_to_complex(Expr* e, double _Complex* out) {
    if (!e) return false;
    double rv;
    if (nl_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double r, i;
        if (nl_to_double_real(re, &r) && nl_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

static Expr* nl_from_complex_d(double _Complex c) {
    double r = creal(c), i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
static Expr* nl_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}

/* (out_re, out_im) = (a_re + a_im i)(b_re + b_im i). Output may alias inputs. */
static void nl_cmul(mpfr_t out_re, mpfr_t out_im,
                    const mpfr_t a_re, const mpfr_t a_im,
                    const mpfr_t b_re, const mpfr_t b_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t t1, t2, t3, t4;
    mpfr_inits2(p, t1, t2, t3, t4, (mpfr_ptr)0);
    mpfr_mul(t1, a_re, b_re, MPFR_RNDN);
    mpfr_mul(t2, a_im, b_im, MPFR_RNDN);
    mpfr_mul(t3, a_re, b_im, MPFR_RNDN);
    mpfr_mul(t4, a_im, b_re, MPFR_RNDN);
    mpfr_sub(out_re, t1, t2, MPFR_RNDN);
    mpfr_add(out_im, t3, t4, MPFR_RNDN);
    mpfr_clears(t1, t2, t3, t4, (mpfr_ptr)0);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors nderiv.c)                    *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* name;
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} NlBind;

static void nl_bind_snapshot(NlBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void nl_bind_clear_temp(SymbolDef* def) {
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

static void nl_bind_set(NlBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    nl_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

static void nl_bind_restore(NlBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    nl_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Sampler: evaluate expr with the variable bound to a complex point  *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*       f;     /* borrowed expr (free in the variable) */
    NlBind*     bind;
    NumericSpec spec;
} NlCtx;

/* Evaluate f after binding the variable to `value` (consumed). Returns the
 * numericalised result (caller frees) or NULL. */
static Expr* nl_eval_at(NlCtx* c, Expr* value) {
    nl_bind_set(c->bind, value);   /* takes a copy internally */
    expr_free(value);
    eval_clock_bump();
    Expr* raw = eval_and_free(expr_copy(c->f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, c->spec);
    expr_free(raw);
    return num;
}

static bool nl_sample_machine(NlCtx* c, double _Complex z, double _Complex* out) {
    Expr* num = nl_eval_at(c, nl_from_complex_d(z));
    if (!num) return false;
    bool ok = nl_to_complex(num, out);
    if (ok && (!isfinite(creal(*out)) || !isfinite(cimag(*out)))) ok = false;
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nl_sample_mpfr(NlCtx* c, const mpfr_t z_re, const mpfr_t z_im,
                           mpfr_t out_re, mpfr_t out_im) {
    Expr* num = nl_eval_at(c, nl_from_complex_mpfr(z_re, z_im));
    if (!num) return false;
    bool inexact;
    bool ok = get_approx_mpfr(num, out_re, out_im, &inexact);
    if (ok && (!mpfr_number_p(out_re) || !mpfr_number_p(out_im))) ok = false;
    expr_free(num);
    return ok;
}
#endif

/* Numericalise a constant (var-free) expression to a machine complex value. */
static bool nl_const_complex_machine(Expr* e, double _Complex* out) {
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
    expr_free(raw);
    bool ok = num && nl_to_complex(num, out);
    expr_free(num);
    return ok;
}

#ifdef USE_MPFR
static bool nl_const_complex_mpfr(Expr* e, long bits, mpfr_t re, mpfr_t im) {
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = bits;
    spec.preserve_inexact = false;
    Expr* raw = eval_and_free(expr_copy(e));
    Expr* num = raw ? numericalize(raw, spec) : NULL;
    expr_free(raw);
    bool inexact;
    bool ok = num && get_approx_mpfr(num, re, im, &inexact);
    expr_free(num);
    return ok;
}
#endif

/* ------------------------------------------------------------------ *
 *  Limit-point classification: finite vs. an infinite ray            *
 * ------------------------------------------------------------------ */

/* If z0 denotes an infinite limit point, return true and store in *ray an
 * owned expression for the (un-normalised) direction of the ray from the
 * origin.  Otherwise return false and leave *ray untouched.
 *
 * Recognised: Infinity (+1), ComplexInfinity / DirectedInfinity[] (+1),
 * DirectedInfinity[d] (d), and any Times[...] containing an Infinity or a
 * DirectedInfinity factor (e.g. -Infinity, I Infinity) whose ray is the
 * product of the remaining factors with the directed-infinity direction. */
static bool nl_infinite_ray(Expr* z0, Expr** ray) {
    if (!z0) return false;
    if (z0->type == EXPR_SYMBOL) {
        if (z0->data.symbol.name == SYM_Infinity
            || z0->data.symbol.name == SYM_ComplexInfinity) {
            *ray = expr_new_integer(1);
            return true;
        }
        return false;
    }
    if (z0->type != EXPR_FUNCTION) return false;
    Expr* head = z0->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;

    if (head->data.symbol.name == SYM_DirectedInfinity) {
        if (z0->data.function.arg_count == 1) {
            *ray = expr_copy(z0->data.function.args[0]);
        } else {
            *ray = expr_new_integer(1);
        }
        return true;
    }

    if (head->data.symbol.name == SYM_Times) {
        size_t n = z0->data.function.arg_count;
        bool found = false;
        Expr* prod = NULL;          /* accumulating Times of the remaining factors */
        for (size_t i = 0; i < n; i++) {
            Expr* f = z0->data.function.args[i];
            Expr* contrib = NULL;
            if (f->type == EXPR_SYMBOL
                && (f->data.symbol.name == SYM_Infinity
                    || f->data.symbol.name == SYM_ComplexInfinity)) {
                found = true;
                continue;           /* Infinity contributes unit magnitude */
            }
            if (f->type == EXPR_FUNCTION
                && f->data.function.head->type == EXPR_SYMBOL
                && f->data.function.head->data.symbol.name == SYM_DirectedInfinity) {
                found = true;
                contrib = (f->data.function.arg_count == 1)
                            ? expr_copy(f->data.function.args[0])
                            : expr_new_integer(1);
            } else {
                contrib = expr_copy(f);
            }
            if (!prod) {
                prod = contrib;
            } else {
                Expr** v = malloc(sizeof(Expr*) * 2);
                v[0] = prod; v[1] = contrib;
                prod = expr_new_function(expr_new_symbol(SYM_Times), v, 2);
                free(v);
            }
        }
        if (!found) { expr_free(prod); return false; }
        *ray = prod ? prod : expr_new_integer(1);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Extrapolation kernels live in seqaccel.{c,h} (shared with NSum).   *
 *  NLimit keeps only its own acceptance/noise policy (nl_accept) and  *
 *  the Expr construction around the raw extrapolated value.           *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* L1 magnitude |re| + |im| of (re, im) as a double, for the gauges. */
static double nl_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Oscillatory-divergence gate (see the NL_OSC_* block above)          *
 * ------------------------------------------------------------------ */

/* Number of times the sample-to-sample increment turns by more than 90 degrees.
 * A monotone or smoothly-converging sequence scores 0. */
static int nl_count_reversals(const double _Complex* S, int n) {
    int rev = 0;
    for (int k = 2; k < n; k++) {
        double _Complex a = S[k] - S[k - 1];
        double _Complex b = S[k - 1] - S[k - 2];
        if (creal(a * conj(b)) < 0.0) rev++;
    }
    return rev;
}

/* mags[] holds |f| along the diagnostic ladder, ordered furthest-from ->
 * nearest-to the limit point.  Stores the near/far envelope ratio and returns
 * true, or returns false to abstain (too little usable data). */
static bool nl_envelope_ratio(const double* mags, int n, double* out) {
    if (n < NL_OSC_MIN_POINTS) return false;
    int h = n / 2;
    double far = 0.0, near = 0.0;
    for (int i = 0; i < h; i++)     if (mags[i] > far)  far  = mags[i];
    for (int i = n - h; i < n; i++) if (mags[i] > near) near = mags[i];
    if (!isfinite(far) || !isfinite(near) || far <= 0.0) return false;
    *out = near / far;
    return true;
}

/* The k-th ladder point carries the factor phi on one of its two visits; taking
 * phi first when marching outward and second when closing in on a finite point
 * makes the merged sequence run furthest -> nearest without a sort. */
static bool nl_osc_use_phi(bool infinite, int visit) {
    return infinite ? (visit == 0) : (visit == 1);
}

/* True if the sampled function oscillates with a non-decaying envelope. */
static bool nl_oscillates_machine(NlCtx* c, double _Complex z0,
                                  double _Complex dir, double _Complex scale,
                                  bool infinite) {
    double mags[2 * NL_OSC_OCTAVES];
    int n = 0;
    for (int k = 0; k < NL_OSC_OCTAVES; k++) {
        for (int visit = 0; visit < 2; visit++) {
            double h = ldexp(1.0, infinite ? k : -k);
            if (nl_osc_use_phi(infinite, visit)) h *= NL_OSC_PHI;
            double _Complex z = infinite ? dir * scale * h
                                         : z0 - dir * scale * h;
            double _Complex v;
            if (nl_sample_machine(c, z, &v)) mags[n++] = cabs(v);
        }
    }
    double env;
    if (!nl_envelope_ratio(mags, n, &env)) return false;
    return env > NL_OSC_ENV_MAX;
}

#ifdef USE_MPFR
/* MPFR analogue.  The ladder is evaluated at the caller's working precision --
 * measuring the envelope of a high-precision integrand in machine-precision
 * noise would be its own false-positive source -- but the ratio itself only
 * needs a double.  (z0r, z0i) is the finite limit point and (dsr, dsi) the
 * unit-direction * Scale product, exactly as in nl_run_mpfr's sampling loop. */
static bool nl_oscillates_mpfr(NlCtx* c, const mpfr_t z0r, const mpfr_t z0i,
                               const mpfr_t dsr, const mpfr_t dsi,
                               bool infinite, long bits) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t zr, zi, hk, vr, vi, phi;
    mpfr_inits2(p, zr, zi, hk, vr, vi, phi, (mpfr_ptr)0);
    mpfr_set_d(phi, NL_OSC_PHI, MPFR_RNDN);

    double mags[2 * NL_OSC_OCTAVES];
    int n = 0;
    for (int k = 0; k < NL_OSC_OCTAVES; k++) {
        for (int visit = 0; visit < 2; visit++) {
            mpfr_set_ui(hk, 1, MPFR_RNDN);
            if (infinite) mpfr_mul_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);
            else          mpfr_div_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);
            if (nl_osc_use_phi(infinite, visit))
                mpfr_mul(hk, hk, phi, MPFR_RNDN);
            mpfr_mul(zr, dsr, hk, MPFR_RNDN);
            mpfr_mul(zi, dsi, hk, MPFR_RNDN);
            if (!infinite) {
                mpfr_sub(zr, z0r, zr, MPFR_RNDN);
                mpfr_sub(zi, z0i, zi, MPFR_RNDN);
            }
            if (nl_sample_mpfr(c, zr, zi, vr, vi))
                mags[n++] = hypot(mpfr_get_d(vr, MPFR_RNDN),
                                  mpfr_get_d(vi, MPFR_RNDN));
        }
    }
    mpfr_clears(zr, zi, hk, vr, vi, phi, (mpfr_ptr)0);

    double env;
    if (!nl_envelope_ratio(mags, n, &env)) return false;
    return env > NL_OSC_ENV_MAX;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */

typedef struct {
    const char* method;    /* SYM_Automatic (default) / SYM_EulerSum / SYM_SequenceLimit / SYM_Levin */
    int         levin_variant; /* SeqaccelLevinVariant when method == SYM_Levin (default u) */
    Expr*       direction; /* borrowed Direction value (NULL => Automatic) */
    Expr*       scale;     /* borrowed Scale value (NULL => 1)             */
    int         terms;     /* starting sample count / tableau depth        */
    int         wynn;      /* WynnDegree (default 1)                       */
    bool        prec_mpfr; /* WorkingPrecision selects MPFR                */
    long        bits;      /* MPFR working precision in bits               */
    double      acc_goal;  /* AccuracyGoal digits (nc_accuracy.h sentinels) */
    double      prec_goal; /* PrecisionGoal digits (nc_accuracy.h sentinels)*/
} NlOpts;

static bool nl_is_known_option(const char* s) {
    return s == SYM_Method
        || s == SYM_Direction
        || s == SYM_Scale
        || s == SYM_Terms
        || s == SYM_WynnDegree
        || s == SYM_WorkingPrecision
        || s == SYM_AccuracyGoal
        || s == SYM_PrecisionGoal;
}

static bool nl_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && nl_is_known_option(lhs->data.symbol.name);
}

static bool nl_parse_working_precision(Expr* val, bool* mpfr, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
        *mpfr = false; *bits = 0; return true;
    }
    double digits;
    if (!nl_to_double_real(val, &digits) || digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) { *mpfr = false; *bits = 0; }
    else { *mpfr = true; *bits = numeric_digits_to_bits(digits); }
#else
    *mpfr = false; *bits = 0;
#endif
    return true;
}

static bool nl_apply_option(Expr* rule, NlOpts* o) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL
            && (rhs->data.symbol.name == SYM_EulerSum
                || rhs->data.symbol.name == SYM_SequenceLimit)) {
            o->method = rhs->data.symbol.name; return true;
        }
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
            o->method = SYM_Automatic; return true;
        }
        /* Levin transform: bare symbol Levin (= u), or string "Levin"/"LevinU"/
         * "LevinT"/"LevinV" selecting the remainder-estimate variant. */
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Levin) {
            o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_U; return true;
        }
        if (rhs->type == EXPR_STRING) {
            const char* s = rhs->data.string;
            if (!strcmp(s, "Levin") || !strcmp(s, "LevinU")) {
                o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_U; return true;
            }
            if (!strcmp(s, "LevinT")) {
                o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_T; return true;
            }
            if (!strcmp(s, "LevinV")) {
                o->method = SYM_Levin; o->levin_variant = SEQACCEL_LEVIN_V; return true;
            }
            if (!strcmp(s, "EulerSum"))      { o->method = SYM_EulerSum;      return true; }
            if (!strcmp(s, "SequenceLimit")) { o->method = SYM_SequenceLimit; return true; }
            if (!strcmp(s, "Automatic"))     { o->method = SYM_Automatic;     return true; }
        }
        nl_warn("badmeth", "Method must be Automatic, EulerSum, SequenceLimit or "
                           "\"Levin\"; using Automatic");
        o->method = SYM_Automatic;
        return true;
    }
    if (name == SYM_Direction) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
            o->direction = NULL;       /* Automatic == -1 */
        else
            o->direction = rhs;
        return true;
    }
    if (name == SYM_Scale) { o->scale = rhs; return true; }
    if (name == SYM_Terms) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 2) {
            o->terms = (int)rhs->data.integer;
            if (o->terms > NL_MAX_TERMS) o->terms = NL_MAX_TERMS;
            return true;
        }
        nl_warn("badopt", "Terms must be an integer >= 2");
        return false;
    }
    if (name == SYM_WynnDegree) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer >= 1) {
            o->wynn = (int)rhs->data.integer;
            return true;
        }
        nl_warn("badopt", "WynnDegree must be a positive integer");
        return false;
    }
    if (name == SYM_WorkingPrecision) {
        if (!nl_parse_working_precision(rhs, &o->prec_mpfr, &o->bits)) {
            nl_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_AccuracyGoal) {
        if (!nc_parse_goal(rhs, &o->acc_goal)) {
            nl_warn("badopt", "AccuracyGoal must be a positive number, "
                              "Automatic, Infinity or MachinePrecision");
            return false;
        }
        return true;
    }
    if (name == SYM_PrecisionGoal) {
        if (!nc_parse_goal(rhs, &o->prec_goal)) {
            nl_warn("badopt", "PrecisionGoal must be a positive number, "
                              "Automatic, Infinity or MachinePrecision");
            return false;
        }
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Extrapolation dispatch (Automatic = best-of Richardson + Wynn)      *
 * ------------------------------------------------------------------ */

/* Automatic method: try Richardson and Wynn's epsilon at every admissible
 * degree, and return the estimate whose internal convergence residual (*step)
 * is smallest, i.e. the extrapolant that is most self-consistent.  Richardson's
 * diagonal step of exactly zero means genuine exact convergence (as for an
 * integer-power tail such as 1 - 2^-k) and is fully trusted; a Wynn residual of
 * exactly zero is instead treated as unreliable — it arises from a degenerate
 * single-entry column or a singular-bridge coincidence, not from certified
 * convergence — and is used only as a last resort.  Returns false if no
 * candidate could be built. */
static bool nl_auto_machine(const double _Complex* S, int terms,
                            double _Complex* result, double* step) {
    double _Complex best_r = 0.0;   double best_s = INFINITY; bool have = false;
    double _Complex first_r = 0.0;  bool have_first = false;
    double _Complex r; double s;

    /* Richardson: eligible for any finite step, including exact zero. */
    if (seqaccel_richardson_machine(S, terms, &r, &s)
        && isfinite(creal(r)) && isfinite(cimag(r))) {
        if (!have_first) { first_r = r; have_first = true; }
        if (isfinite(s) && s < best_s) { best_r = r; best_s = s; have = true; }
    }
    /* Wynn: require a strictly positive residual (a certified convergence gap). */
    int maxdeg = (terms - 1) / 2;
    for (int d = 1; d <= maxdeg; d++) {
        if (!seqaccel_wynn_machine(S, terms, d, &r, &s)) continue;
        if (!isfinite(creal(r)) || !isfinite(cimag(r))) continue;
        if (!have_first) { first_r = r; have_first = true; }
        if (isfinite(s) && s > 0.0 && s < best_s) {
            best_r = r; best_s = s; have = true;
        }
    }
    /* Levin u: same "positive residual = certified gap" policy as Wynn.  Strong
     * on logarithmically/algebraically convergent samples the tableau misses.
     * Only admitted when the sample sequence is actually settling (its last
     * increment is smaller than its first): on a divergent sequence Levin can
     * collapse to a spurious finite value with a deceptively small residual, so
     * this gate keeps it from hijacking the best-of. */
    if (terms >= 3 && cabs(S[terms - 1] - S[terms - 2]) < cabs(S[1] - S[0])
        && seqaccel_levin_machine(S, terms, SEQACCEL_LEVIN_U, 1.0, &r, &s)
        && isfinite(creal(r)) && isfinite(cimag(r))) {
        if (!have_first) { first_r = r; have_first = true; }
        if (isfinite(s) && s > 0.0 && s < best_s) { best_r = r; best_s = s; have = true; }
    }
    if (have)       { *result = best_r;  *step = best_s; return true; }
    if (have_first) { *result = first_r; *step = 0.0;    return true; }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Sequence construction + dispatch                                   *
 * ------------------------------------------------------------------ */

/* Dispatch the chosen extrapolation Method over S[0..terms-1]. */
static bool nl_extrapolate(const double _Complex* S, int terms, const NlOpts* o,
                           double _Complex* result, double* step) {
    if (o->method == SYM_SequenceLimit)
        return seqaccel_wynn_machine(S, terms, o->wynn, result, step);
    if (o->method == SYM_EulerSum)
        return seqaccel_richardson_machine(S, terms, result, step);
    if (o->method == SYM_Levin)
        return seqaccel_levin_machine(S, terms, o->levin_variant, 1.0, result, step);
    return nl_auto_machine(S, terms, result, step);        /* Automatic */
}

/* True unless BOTH goals are Infinity, i.e. accuracy is not used as a
 * termination criterion at all (in which case Terms is honoured verbatim). */
static bool nl_goal_active(const NlOpts* o) {
    return !(isinf(o->acc_goal) && isinf(o->prec_goal));
}

/* Machine path.  Samples an octave-spaced sequence and accelerates it, growing
 * the sample count from o->terms up to NL_MAX_TERMS until the combined
 * AccuracyGoal/PrecisionGoal tolerance is met, keeping the estimate with the
 * smallest internal residual.  On a genuine non-limit (divergent / noise-
 * dominated) it returns NULL as before; when a limit is recognised but the goal
 * is not reached at the cap, it warns and returns the best approximation. */
static Expr* nl_run_machine(Expr* expr, const char* var, double _Complex z0,
                            double _Complex dir, double _Complex scale,
                            bool infinite, NlOpts* o) {
    const double wp_digits = NUMERIC_MACHINE_PRECISION_DIGITS;
    const bool goal_active = nl_goal_active(o);
    int start = o->terms; if (start < 2) start = 2;
    int cap = goal_active ? NL_MAX_TERMS : start;
    if (cap < start) cap = start;

    double _Complex* S = malloc(sizeof(double _Complex) * (size_t)cap);
    if (!S) return NULL;

    NlBind bind; nl_bind_snapshot(&bind, var);
    NlCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = numeric_machine_spec();

    /* Sample the starting prefix. */
    int filled = 0;
    bool ok = true;
    for (; filled < start; filled++) {
        double _Complex z;
        if (infinite) z = dir * scale * ldexp(1.0, filled);     /* u * Scale * 2^k */
        else          z = z0 - dir * scale * ldexp(1.0, -filled); /* z0 - d*Scale*2^-k */
        if (!nl_sample_machine(&ctx, z, &S[filled])) { ok = false; break; }
    }
    /* The oscillation diagnosis samples expr too, so it must run while the
     * variable is still bound -- before nl_bind_restore, not after. */
    bool oscillatory = ok
        && nl_count_reversals(S, filled) >= NL_OSC_MIN_REVERSALS
        && nl_oscillates_machine(&ctx, z0, dir, scale, infinite);

    /* Adaptive refinement: extrapolate at the current prefix, keep the estimate
     * with the smallest residual, grow the sample count until the goal is met
     * or the cap is reached. */
    double _Complex best_r = 0.0; double best_step = INFINITY;
    double best_maxsample = 0.0; bool have_best = false;
    if (ok && !oscillatory) {
        int terms = filled;
        int stall = 0;
        for (;;) {
            double maxsample = 0.0;
            for (int k = 0; k < terms; k++) {
                double v = cabs(S[k]);
                if (v > maxsample) maxsample = v;
            }
            double _Complex result; double step = 0.0;
            bool improved = false;
            if (nl_extrapolate(S, terms, o, &result, &step)
                && isfinite(creal(result)) && isfinite(cimag(result))
                && step < best_step
                && (!have_best
                    || cabs(result - best_r) <= NL_REFINE_JUMP * best_step + 1e-300)) {
                best_r = result; best_step = step;
                best_maxsample = maxsample; have_best = true; improved = true;
            }
            if (!goal_active) break;             /* honour Terms verbatim */
            if (have_best) {
                double tol = nc_combined_tol(o->acc_goal, o->prec_goal,
                                             cabs(best_r), wp_digits);
                if (best_step <= tol
                    && nl_accept(cabs(best_r), best_step, best_maxsample))
                    break;                       /* goal met and accepted */
            }
            stall = improved ? 0 : (stall + 1);
            if (stall >= NL_REFINE_STALL) break; /* refinement no longer helping */
            if (terms >= cap) break;             /* out of refinement budget */
            double _Complex z;                   /* grow by one octave */
            if (infinite) z = dir * scale * ldexp(1.0, terms);
            else          z = z0 - dir * scale * ldexp(1.0, -terms);
            if (!nl_sample_machine(&ctx, z, &S[terms])) break;  /* can't extend */
            terms++; filled = terms;
        }
    }
    nl_bind_restore(&bind);
    free(S);

    if (!ok) {
        nl_warn("notnum", "the expression is not numerical at a sample point");
        return NULL;
    }
    if (oscillatory) { nl_warn("osc", NL_OSC_MESSAGE); return NULL; }
    if (!have_best) { nl_warn("ndterm", "not enough Terms for the chosen Method"); return NULL; }

    /* Divergent / noise-dominated extrapolant: no limiting value recognised. */
    if (!nl_accept(cabs(best_r), best_step, best_maxsample)) {
        nl_warn("noise", "Cannot recognize a limiting value. This may be due to "
                         "noise from roundoff; try higher WorkingPrecision, fewer "
                         "Terms, or a different Scale.");
        return NULL;
    }
    /* Accepted, but possibly short of the requested accuracy. */
    if (goal_active) {
        double tol = nc_combined_tol(o->acc_goal, o->prec_goal, cabs(best_r), wp_digits);
        if (best_step > tol) nc_warn_goal("NLimit", best_step, tol);
    }
    return nl_from_complex_d(best_r);
}

#ifdef USE_MPFR
/* MPFR analogue of nl_auto_machine: best-of Richardson + Wynn (all degrees) by
 * smallest internal residual.  (out_re, out_im) are pre-initialised by the
 * caller.  Same zero-residual-is-unreliable policy as the machine path. */
static bool nl_auto_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                         long bits, mpfr_t out_re, mpfr_t out_im,
                         double* step, bool* finite) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t tr, ti;
    mpfr_init2(tr, p); mpfr_init2(ti, p);
    double best_s = INFINITY; bool have = false, have_first = false;
    double s; bool fin;

    if (seqaccel_richardson_mpfr(Sr, Si, terms, bits, tr, ti, &s, &fin) && fin) {
        if (!have_first) {
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; *step = s; have_first = true;
        }
        if (isfinite(s) && s < best_s) {   /* exact-zero step = exact convergence */
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; best_s = s; *step = s; have = true;
        }
    }
    int maxdeg = (terms - 1) / 2;
    for (int d = 1; d <= maxdeg; d++) {
        if (!seqaccel_wynn_mpfr(Sr, Si, terms, d, bits, tr, ti, &s, &fin)) continue;
        if (!fin) continue;
        if (!have_first) {
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; *step = s; have_first = true;
        }
        if (isfinite(s) && s > 0.0 && s < best_s) {
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; best_s = s; *step = s; have = true;
        }
    }
    /* Levin u: same positive-residual policy as Wynn, gated on the sample
     * sequence actually settling (see nl_auto_machine for the rationale). */
    bool contracting = false;
    if (terms >= 3) {
        mpfr_t d0, dn; mpfr_init2(d0, p); mpfr_init2(dn, p);
        mpfr_sub(d0, Sr[1], Sr[0], MPFR_RNDN);
        double first_inc = fabs(mpfr_get_d(d0, MPFR_RNDN));
        mpfr_sub(d0, Si[1], Si[0], MPFR_RNDN);
        first_inc += fabs(mpfr_get_d(d0, MPFR_RNDN));
        mpfr_sub(dn, Sr[terms - 1], Sr[terms - 2], MPFR_RNDN);
        double last_inc = fabs(mpfr_get_d(dn, MPFR_RNDN));
        mpfr_sub(dn, Si[terms - 1], Si[terms - 2], MPFR_RNDN);
        last_inc += fabs(mpfr_get_d(dn, MPFR_RNDN));
        contracting = last_inc < first_inc;
        mpfr_clear(d0); mpfr_clear(dn);
    }
    if (contracting
        && seqaccel_levin_mpfr(Sr, Si, terms, SEQACCEL_LEVIN_U, 1.0, bits, tr, ti, &s, &fin)
        && fin) {
        if (!have_first) {
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; *step = s; have_first = true;
        }
        if (isfinite(s) && s > 0.0 && s < best_s) {
            mpfr_set(out_re, tr, MPFR_RNDN); mpfr_set(out_im, ti, MPFR_RNDN);
            *finite = true; best_s = s; *step = s; have = true;
        }
    }
    mpfr_clear(tr); mpfr_clear(ti);
    return have || have_first;
}

/* Dispatch the chosen extrapolation Method over the MPFR prefix Sr/Si[0..terms). */
static bool nl_extrapolate_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                                long bits, const NlOpts* o,
                                mpfr_t rr, mpfr_t ri, double* step, bool* finite) {
    if (o->method == SYM_SequenceLimit)
        return seqaccel_wynn_mpfr(Sr, Si, terms, o->wynn, bits, rr, ri, step, finite);
    if (o->method == SYM_EulerSum)
        return seqaccel_richardson_mpfr(Sr, Si, terms, bits, rr, ri, step, finite);
    if (o->method == SYM_Levin)
        return seqaccel_levin_mpfr(Sr, Si, terms, o->levin_variant, 1.0, bits, rr, ri, step, finite);
    return nl_auto_mpfr(Sr, Si, terms, bits, rr, ri, step, finite);   /* Automatic */
}

/* Sample expr at octave index k into (Srk, Sik); mirrors the machine path's
 * z_k formula.  zr/zi/hk are caller-owned scratch at working precision. */
static bool nl_sample_index_mpfr(NlCtx* ctx, int k, bool infinite,
                                 const mpfr_t z0r, const mpfr_t z0i,
                                 const mpfr_t dsr, const mpfr_t dsi,
                                 mpfr_t zr, mpfr_t zi, mpfr_t hk,
                                 mpfr_t Srk, mpfr_t Sik) {
    if (infinite) {
        mpfr_set_ui(hk, 1, MPFR_RNDN);
        mpfr_mul_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);   /* 2^k */
        mpfr_mul(zr, dsr, hk, MPFR_RNDN);
        mpfr_mul(zi, dsi, hk, MPFR_RNDN);
    } else {
        mpfr_set_ui(hk, 1, MPFR_RNDN);
        mpfr_div_2ui(hk, hk, (unsigned long)k, MPFR_RNDN);   /* 2^-k */
        mpfr_mul(zr, dsr, hk, MPFR_RNDN);
        mpfr_mul(zi, dsi, hk, MPFR_RNDN);
        mpfr_sub(zr, z0r, zr, MPFR_RNDN);                    /* z0 - d*Scale*2^-k */
        mpfr_sub(zi, z0i, zi, MPFR_RNDN);
    }
    return nl_sample_mpfr(ctx, zr, zi, Srk, Sik);
}

static Expr* nl_run_mpfr(Expr* expr, const char* var, Expr* z0_expr,
                         Expr* dir_expr, Expr* scale_expr, bool infinite,
                         Expr* ray_expr, NlOpts* o) {
    long bits = o->bits;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; spec.preserve_inexact = false;

    /* Resolve z0 (finite) or ray (infinite), Direction and Scale to complex. */
    mpfr_t z0r, z0i, dr, di, sr, si, mag;
    mpfr_init2(z0r, p); mpfr_init2(z0i, p);
    mpfr_init2(dr, p);  mpfr_init2(di, p);
    mpfr_init2(sr, p);  mpfr_init2(si, p);
    mpfr_init2(mag, p);
    bool prep = true;

    if (infinite) {
        prep = nl_const_complex_mpfr(ray_expr, bits, dr, di);
        mpfr_set_ui(z0r, 0, MPFR_RNDN); mpfr_set_ui(z0i, 0, MPFR_RNDN);
    } else {
        prep = nl_const_complex_mpfr(z0_expr, bits, z0r, z0i);
        if (prep) {
            if (dir_expr) prep = nl_const_complex_mpfr(dir_expr, bits, dr, di);
            else { mpfr_set_si(dr, -1, MPFR_RNDN); mpfr_set_ui(di, 0, MPFR_RNDN); }
        }
    }
    /* Normalise the direction/ray to a unit vector. */
    if (prep) {
        mpfr_hypot(mag, dr, di, MPFR_RNDN);
        if (mpfr_zero_p(mag)) prep = false;
        else { mpfr_div(dr, dr, mag, MPFR_RNDN); mpfr_div(di, di, mag, MPFR_RNDN); }
    }
    if (prep) {
        if (scale_expr) prep = nl_const_complex_mpfr(scale_expr, bits, sr, si);
        else { mpfr_set_ui(sr, 1, MPFR_RNDN); mpfr_set_ui(si, 0, MPFR_RNDN); }
    }
    if (prep && mpfr_zero_p(sr) && mpfr_zero_p(si)) prep = false;
    if (!prep) {
        mpfr_clears(z0r, z0i, dr, di, sr, si, mag, (mpfr_ptr)0);
        nl_warn("notnum", "the limit point, Direction or Scale is not numerical");
        return NULL;
    }

    const bool goal_active = nl_goal_active(o);
    const double wp_digits = numeric_bits_to_digits(bits);
    int start = o->terms; if (start < 2) start = 2;
    int cap = goal_active ? NL_MAX_TERMS : start;
    if (cap < start) cap = start;

    mpfr_t* Sr = malloc(sizeof(mpfr_t) * (size_t)cap);
    mpfr_t* Si = malloc(sizeof(mpfr_t) * (size_t)cap);
    for (int k = 0; k < cap; k++) { mpfr_init2(Sr[k], p); mpfr_init2(Si[k], p); }
    /* scratch for the point and the d*Scale product */
    mpfr_t dsr, dsi, zr, zi, hk;
    mpfr_inits2(p, dsr, dsi, zr, zi, hk, (mpfr_ptr)0);
    nl_cmul(dsr, dsi, dr, di, sr, si);          /* d * Scale (unit dir * scale) */

    NlBind bind; nl_bind_snapshot(&bind, var);
    NlCtx ctx; ctx.f = expr; ctx.bind = &bind; ctx.spec = spec;

    /* Sample the starting prefix. */
    int filled = 0;
    bool ok = true;
    while (filled < start) {
        if (!nl_sample_index_mpfr(&ctx, filled, infinite, z0r, z0i, dsr, dsi,
                                  zr, zi, hk, Sr[filled], Si[filled])) { ok = false; break; }
        filled++;
    }
    /* Same ordering constraint as the machine path: diagnose before unbinding.
     * The reversal screen is a sign test, so it reads the samples as doubles. */
    bool oscillatory = false;
    if (ok) {
        double _Complex Sd[NL_MAX_TERMS];
        for (int k = 0; k < filled; k++)
            Sd[k] = mpfr_get_d(Sr[k], MPFR_RNDN)
                  + mpfr_get_d(Si[k], MPFR_RNDN) * I;
        oscillatory = nl_count_reversals(Sd, filled) >= NL_OSC_MIN_REVERSALS
                   && nl_oscillates_mpfr(&ctx, z0r, z0i, dsr, dsi, infinite, bits);
    }

    /* Adaptive refinement: keep the smallest-residual estimate, grow the sample
     * count until the combined goal tolerance is met or the cap is reached. */
    mpfr_t best_re, best_im, rr, ri;
    mpfr_init2(best_re, p); mpfr_init2(best_im, p);
    mpfr_init2(rr, p);      mpfr_init2(ri, p);
    double best_step = INFINITY, best_maxsample = 0.0;
    bool have_best = false;
    if (ok && !oscillatory) {
        int terms = filled;
        int stall = 0;
        for (;;) {
            double maxsample = 0.0;
            for (int k = 0; k < terms; k++) {
                double v = nl_l1_d(Sr[k], Si[k]);
                if (v > maxsample) maxsample = v;
            }
            double step = 0.0; bool finite = false;
            bool improved = false;
            if (nl_extrapolate_mpfr(Sr, Si, terms, bits, o, rr, ri, &step, &finite)
                && finite && isfinite(step) && step < best_step) {
                double jump = 0.0;               /* value-consistency guard */
                if (have_best)
                    jump = hypot(mpfr_get_d(rr, MPFR_RNDN) - mpfr_get_d(best_re, MPFR_RNDN),
                                 mpfr_get_d(ri, MPFR_RNDN) - mpfr_get_d(best_im, MPFR_RNDN));
                if (!have_best || jump <= NL_REFINE_JUMP * best_step + 1e-300) {
                    mpfr_set(best_re, rr, MPFR_RNDN); mpfr_set(best_im, ri, MPFR_RNDN);
                    best_step = step; best_maxsample = maxsample;
                    have_best = true; improved = true;
                }
            }
            if (!goal_active) break;
            if (have_best) {
                double bmag = nl_l1_d(best_re, best_im);
                double tol = nc_combined_tol(o->acc_goal, o->prec_goal, bmag, wp_digits);
                if (best_step <= tol && nl_accept(bmag, best_step, best_maxsample))
                    break;
            }
            stall = improved ? 0 : (stall + 1);
            if (stall >= NL_REFINE_STALL) break;
            if (terms >= cap) break;
            if (!nl_sample_index_mpfr(&ctx, terms, infinite, z0r, z0i, dsr, dsi,
                                      zr, zi, hk, Sr[terms], Si[terms])) break;
            terms++; filled = terms;
        }
    }
    nl_bind_restore(&bind);

    Expr* out = NULL;
    if (!ok) {
        nl_warn("notnum", "the expression is not numerical at a sample point");
    } else if (oscillatory) {
        nl_warn("osc", NL_OSC_MESSAGE);
    } else if (have_best
               && nl_accept(nl_l1_d(best_re, best_im), best_step, best_maxsample)) {
        if (goal_active) {
            double bmag = nl_l1_d(best_re, best_im);
            double tol = nc_combined_tol(o->acc_goal, o->prec_goal, bmag, wp_digits);
            if (best_step > tol) nc_warn_goal("NLimit", best_step, tol);
        }
        out = nl_from_complex_mpfr(best_re, best_im);
    } else {
        nl_warn("noise", "Cannot recognize a limiting value. This may be due "
                         "to noise from roundoff; try higher WorkingPrecision, "
                         "fewer Terms, or a different Scale.");
    }

    mpfr_clear(best_re); mpfr_clear(best_im); mpfr_clear(rr); mpfr_clear(ri);
    for (int k = 0; k < cap; k++) { mpfr_clear(Sr[k]); mpfr_clear(Si[k]); }
    free(Sr); free(Si);
    mpfr_clears(z0r, z0i, dr, di, sr, si, mag, dsr, dsi, zr, zi, hk, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_nlimit(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;     /* expr, z -> z0  (+ options) */

    /* Peel trailing options; exactly two positional args must remain. */
    size_t pos_end = argc;
    while (pos_end > 0 && nl_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nl_is_option_arg(res->data.function.args[i])) {
            nl_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) return NULL;

    /* Parse the z -> z0 rule. */
    Expr* rule = res->data.function.args[1];
    if (rule->type != EXPR_FUNCTION
        || rule->data.function.head->type != EXPR_SYMBOL
        || (rule->data.function.head->data.symbol.name != SYM_Rule
            && rule->data.function.head->data.symbol.name != SYM_RuleDelayed)
        || rule->data.function.arg_count != 2
        || rule->data.function.args[0]->type != EXPR_SYMBOL) {
        return NULL;               /* not a var -> point spec: stay unevaluated */
    }
    const char* var = rule->data.function.args[0]->data.symbol.name;
    Expr* z0_expr = rule->data.function.args[1];
    Expr* expr = res->data.function.args[0];

    NlOpts o;
    o.method = SYM_Automatic;
    o.levin_variant = SEQACCEL_LEVIN_U;
    o.direction = NULL;
    o.scale = NULL;
    o.terms = 13;                                  /* starting sample count */
    o.wynn = 1;
    o.prec_mpfr = false;
    o.bits = 0;
    o.acc_goal = NC_GOAL_MACHINE; /* AccuracyGoal -> MachinePrecision */
    o.prec_goal = NC_GOAL_AUTO;                    /* PrecisionGoal -> Automatic */
    for (size_t i = pos_end; i < argc; i++) {
        if (!nl_apply_option(res->data.function.args[i], &o)) return NULL;
    }

    /* Classify the limit point. */
    Expr* ray = NULL;
    bool infinite = nl_infinite_ray(z0_expr, &ray);

#ifdef USE_MPFR
    if (o.prec_mpfr) {
        Expr* out = nl_run_mpfr(expr, var, z0_expr, o.direction, o.scale,
                                infinite, ray, &o);
        expr_free(ray);
        return out;
    }
#endif

    /* Machine path: resolve z0/ray, Direction and Scale to complex up front. */
    double _Complex z0 = 0.0, dir, scale;
    bool prep;
    if (infinite) {
        prep = nl_const_complex_machine(ray, &dir);
    } else {
        prep = nl_const_complex_machine(z0_expr, &z0);
        if (prep) {
            if (o.direction) prep = nl_const_complex_machine(o.direction, &dir);
            else dir = -1.0 + 0.0 * I;          /* Automatic == -1 */
        }
    }
    if (prep) {
        double m = cabs(dir);
        if (m == 0.0) prep = false; else dir /= m;   /* unit direction/ray */
    }
    if (prep) {
        if (o.scale) prep = nl_const_complex_machine(o.scale, &scale);
        else scale = 1.0 + 0.0 * I;
        if (prep && scale == 0.0) prep = false;
    }
    expr_free(ray);
    if (!prep) {
        /* Non-numeric limit point / Direction / Scale: stay unevaluated unless
         * it was an explicitly bad option. */
        return NULL;
    }

    return nl_run_machine(expr, var, z0, dir, scale, infinite, &o);
}

/* ------------------------------------------------------------------ *
 *  Per-method entry points                                            *
 *                                                                     *
 *  Every Method setting is also callable as a head of its own:        *
 *  NLimit`EulerSum[f, z -> z0] is exactly                             *
 *  NLimit[f, z -> z0, Method -> "EulerSum"].  The two positional      *
 *  arguments and every other option (Direction, Scale, Terms,         *
 *  WynnDegree, WorkingPrecision) are forwarded untouched; see         *
 *  common_method_alias.  The Levin variants differ only in the        *
 *  remainder estimate, so each gets its own head rather than hiding   *
 *  behind a string setting.                                           *
 * ------------------------------------------------------------------ */

#define NLIMIT_METHOD_ENTRY(fn, name)                                    \
    static Expr* fn(Expr* res) {                                         \
        return common_method_alias(res, "NLimit", 2, name, builtin_nlimit); \
    }

NLIMIT_METHOD_ENTRY(builtin_nlimit_automatic,     "Automatic")
NLIMIT_METHOD_ENTRY(builtin_nlimit_eulersum,      "EulerSum")
NLIMIT_METHOD_ENTRY(builtin_nlimit_sequencelimit, "SequenceLimit")
NLIMIT_METHOD_ENTRY(builtin_nlimit_levin,         "Levin")
NLIMIT_METHOD_ENTRY(builtin_nlimit_levinu,        "LevinU")
NLIMIT_METHOD_ENTRY(builtin_nlimit_levint,        "LevinT")
NLIMIT_METHOD_ENTRY(builtin_nlimit_levinv,        "LevinV")

#undef NLIMIT_METHOD_ENTRY

static void nl_register_method(const char* sym, Expr* (*fn)(Expr*),
                               const char* doc) {
    symtab_add_builtin(sym, fn);
    symtab_get_def(sym)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring(sym, doc);
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void nlimit_init(void) {
    symtab_add_builtin("NLimit", builtin_nlimit);
    /* Protected only. NLimit is deliberately NOT ATTR_LISTABLE: threading would
     * split the z -> z0 spec across bogus single-argument calls. */
    symtab_get_def("NLimit")->attributes |= ATTR_PROTECTED;
    /* NLimit's own docstring lives in info.c (registered later in
     * core_init, so it wins regardless); only the per-method heads are
     * documented here, next to the code they describe. */

    /* --- The Method settings, exposed as heads ---------------------- */

    nl_register_method("NLimit`Automatic", builtin_nlimit_automatic,
        "NLimit`Automatic[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "Automatic]: run Richardson extrapolation, Wynn's epsilon algorithm "
        "(at every admissible degree) and Levin's u-transform on the sampled "
        "sequence, and return the estimate with the smallest internal "
        "convergence residual.  Richardson models an integer-power (analytic) "
        "error tail while Wynn captures the geometric and fractional-power "
        "(branch-point) tails its fixed 2^j-1 denominators cannot annihilate, "
        "so selecting by best self-consistency picks the right tool per "
        "problem.  Levin is admitted only when the sample increments are "
        "contracting, since on a divergent sequence it can collapse to a "
        "spurious value with a deceptively small residual.  This is what "
        "plain NLimit does.  (The MPFR path runs Richardson and Wynn only.)");

    nl_register_method("NLimit`EulerSum", builtin_nlimit_eulersum,
        "NLimit`EulerSum[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "EulerSum]: Richardson / Romberg extrapolation of the sampled "
        "sequence S_k treated as a function of the geometric step.  Following "
        "the same convention as ND, the tableau uses the all-powers "
        "denominator 2^j - 1: T(i,0) = S_i, T(i,j) = T(i,j-1) + "
        "(T(i,j-1) - T(i-1,j-1))/(2^j - 1), and the result is T(n-1,n-1) for "
        "n = Terms.  Best on a smooth approach with an analytic error tail; "
        "it cannot annihilate a geometric or fractional-power tail.");

    nl_register_method("NLimit`SequenceLimit", builtin_nlimit_sequencelimit,
        "NLimit`SequenceLimit[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "SequenceLimit]: Wynn's epsilon algorithm (iterated Shanks / Aitken) "
        "on the sampled sequence.  With eps(-1) = 0 and eps(0) = S_n, "
        "eps(k+1,n) = eps(k-1,n+1) + 1/(eps(k,n+1) - eps(k,n)); the even "
        "columns are the limit estimates.  WynnDegree -> d selects column 2d "
        "and needs at least 2(d+1) terms for a convergence check.  Handles "
        "the geometric and branch-point tails that defeat Richardson.");

    nl_register_method("NLimit`Levin", builtin_nlimit_levin,
        "NLimit`Levin[f, z -> z0] is NLimit[f, z -> z0, Method -> \"Levin\"]: "
        "the Levin sequence transformation with the u-type remainder "
        "estimate, identical to NLimit`LevinU.  Levin's transformation "
        "divides out an explicit model of the remainder before extrapolating, "
        "which makes it markedly stronger than Richardson on slowly "
        "convergent and alternating sequences.  Use NLimit`LevinT or "
        "NLimit`LevinV for the other remainder estimates.");

    nl_register_method("NLimit`LevinU", builtin_nlimit_levinu,
        "NLimit`LevinU[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "\"LevinU\"]: the Levin u-transform, whose remainder estimate is "
        "k a_k (a_k the k-th term difference).  The general-purpose variant, "
        "and the one NLimit`Levin selects; effective on both logarithmic and "
        "linear convergence.");

    nl_register_method("NLimit`LevinT", builtin_nlimit_levint,
        "NLimit`LevinT[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "\"LevinT\"]: the Levin t-transform, whose remainder estimate is the "
        "term difference a_k itself.  Suited to strictly alternating "
        "sequences, where the remainder is bounded by the first omitted "
        "term; weaker than the u-variant on monotone sequences.");

    nl_register_method("NLimit`LevinV", builtin_nlimit_levinv,
        "NLimit`LevinV[f, z -> z0] is NLimit[f, z -> z0, Method -> "
        "\"LevinV\"]: the Levin v-transform, whose remainder estimate "
        "a_k a_(k+1) / (a_k - a_(k+1)) interpolates between the t- and "
        "u-variants.  Often the most robust of the three when the "
        "convergence type of the sampled sequence is not known in advance.");
}
