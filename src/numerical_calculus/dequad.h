/*
 * dequad.h — double-exponential (exp-sinh) quadrature for a half line
 *
 * Evaluates  I = ∫_a^∞ f(x) dx  for an integrand that decays at +∞ and is
 * analytic on (a, ∞), via the exp-sinh transformation
 *
 *     x(t) = a + exp((π/2) sinh t),     t ∈ (-∞, ∞),
 *     I    = ∫ f(x(t)) x'(t) dt,   x'(t) = (π/2) cosh t · exp((π/2) sinh t),
 *
 * to which the trapezoidal rule converges double-exponentially.  The step h is
 * halved level by level until two successive estimates agree to the requested
 * tolerance; the tails are truncated where the transformed integrand underflows
 * relative to the running maximum.  This is the quadrature Mathematica's NSum
 * uses for the Euler–Maclaurin tail integral (Method -> DoubleExponential), and
 * is reusable groundwork for a future real-axis NIntegrate.
 *
 * The integrand is supplied through a sample callback (machine `double _Complex`
 * or an MPFR real/imag pair), mirroring quadrature.h.  A sample that returns
 * false (non-numeric / non-finite) truncates that tail rather than aborting.
 */
#ifndef MATHILDA_DEQUAD_H
#define MATHILDA_DEQUAD_H

#include <complex.h>
#include <math.h>
#include <stdbool.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

/* Machine sample: write f(x) into *out; return false to truncate the tail. */
typedef bool (*DeQuadSampleMachine)(void* ctx, double x, double _Complex* out);

/* Compute ∫_a^∞ f dx.  Writes *result and, via *abserr, the |last - previous|
 * level difference (an error proxy).  Returns true if the levels settled to
 * `reltol` within `max_levels` refinements, false otherwise (the best estimate
 * is still written). */
bool dequad_halfline_machine(DeQuadSampleMachine f, void* ctx, double a,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr);

#ifdef USE_MPFR
/* MPFR sample: write Re/Im of f(x) into (out_re,out_im) at the working
 * precision; return false to truncate the tail. */
typedef bool (*DeQuadSampleMPFR)(void* ctx, const mpfr_t x,
                                 mpfr_t out_re, mpfr_t out_im);

/* MPFR ∫_a^∞ f dx at `bits` precision.  (out_re,out_im) must be pre-initialised
 * by the caller.  *abserr receives the final level difference as a double. */
bool dequad_halfline_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr);
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Tail-truncation guard shared by the DE drivers (dequad.c, denint.c) *
 * ------------------------------------------------------------------ *
 *
 * Each DE level sums transformed terms outward from t = 0 and may stop a side
 * early on two signals: the terms have fallen below a relative floor (the
 * double-exponential tail is negligible), or, having decayed, they start to
 * *rise* again, the signature of roundoff noise amplified by a growing
 * transform weight (a cancelling summand in the NSum Euler-Maclaurin tail, a
 * noisy integrand at a singular end).
 *
 * The drivers used to read both signals off single samples, and a single sample
 * cannot tell a tail from an interior zero of the integrand.  Exp[-x]
 * BesselJ[0, x]^2 has a double zero at every Bessel zero: a term landing near
 * one is tiny and the next is large, which read as "decayed below 0.1 max,
 * then rose 4x, so noise" and cut away everything beyond the zero, usually the
 * bulk of the integral.  Successive levels cut at nearby zeros and agreed, so
 * the truncated value was reported as converged (0.636 for 0.6426 on [0, oo);
 * 4.3e-6 for 0.6426 at WorkingPrecision -> 20 on [0, 20]).
 *
 * The guard now asks for the evidence each signal actually rests on:
 *   - floor cut: the terms stay <= floor * maxterm over a t-window of
 *     DE_TAIL_WINDOW (a genuine DE tail decays over a stretch of t; an
 *     interior zero is a dip of vanishing t-width);
 *   - noise cut: a rise after the terms have stayed <= 0.1 * maxterm over a
 *     t-window only arms a *checkpoint* (the caller saves its partial sum).
 *     Summing continues, and the rise is confirmed as noise only if the terms
 *     go on to exceed DE_TAIL_NOISE_GROWTH times the envelope of that decayed
 *     window: amplified noise grows without bound, while the lobe after an
 *     interior zero only climbs back to the decaying envelope.  On
 *     confirmation the caller restores the checkpoint and stops the side, which
 *     is exactly the old cut; an unconfirmed rise keeps every term.
 * The window is a length in t, so it costs one extra sample at level 0 and
 * about window/h at step h. */
#define DE_TAIL_WINDOW        0.5
#define DE_TAIL_DECAYED       0.1     /* "decayed": term <= this * maxterm      */
#define DE_TAIL_NOISE_GROWTH  1e4     /* confirm noise: term > this * envelope  */

typedef struct {
    double floor_t;   /* |t| at which the current floor-level run began (<0: none) */
    double dec_t;     /* |t| at which the current decayed run began (<0: none)     */
    double dec_env;   /* largest term in the current decayed run                  */
    double prev_m;    /* magnitude of the previous term on this side              */
    bool   armed;     /* a noise checkpoint is held by the caller                 */
    double arm_env;   /* envelope of the decayed run that armed it                */
} DeTailGuard;

enum {
    DE_TAIL_ADD = 0,        /* keep the term, continue                           */
    DE_TAIL_ADD_STOP,       /* keep the term, stop this side (floor reached)     */
    DE_TAIL_CHECKPOINT,     /* save the partial sum, then keep the term          */
    DE_TAIL_ROLLBACK_STOP   /* noise confirmed: restore the checkpoint and stop  */
};

static inline void de_tail_guard_init(DeTailGuard* g) {
    g->floor_t = -1.0; g->dec_t = -1.0; g->dec_env = 0.0;
    g->prev_m = INFINITY; g->armed = false; g->arm_env = 0.0;
}

/* Classify the term of magnitude m at |t| = at, given the largest term seen so
 * far this level (*maxterm, updated here).  floor_rel is the relative floor. */
static inline int de_tail_step(DeTailGuard* g, double at, double m,
                               double* maxterm, double floor_rel) {
    if (g->armed && m > DE_TAIL_NOISE_GROWTH * g->arm_env)
        return DE_TAIL_ROLLBACK_STOP;
    int act = DE_TAIL_ADD;
    if (!g->armed && g->dec_t >= 0.0 && at - g->dec_t >= DE_TAIL_WINDOW
        && m > 4.0 * g->prev_m) {
        g->armed = true;
        g->arm_env = g->dec_env;
        act = DE_TAIL_CHECKPOINT;
    }
    if (m > *maxterm) *maxterm = m;
    double mx = *maxterm;
    g->prev_m = m;
    if (mx > 0.0 && m <= DE_TAIL_DECAYED * mx) {
        if (g->dec_t < 0.0) { g->dec_t = at; g->dec_env = m; }
        else if (m > g->dec_env) g->dec_env = m;
    } else {
        g->dec_t = -1.0;
    }
    if (mx > 0.0 && m <= floor_rel * mx) {
        if (g->floor_t < 0.0) g->floor_t = at;
        if (at - g->floor_t >= DE_TAIL_WINDOW && act == DE_TAIL_ADD)
            return DE_TAIL_ADD_STOP;
    } else {
        g->floor_t = -1.0;
    }
    return act;
}

#endif /* MATHILDA_DEQUAD_H */
