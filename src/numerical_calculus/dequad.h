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

/* The guard works on log2 magnitudes, so the same classification serves the
 * machine drivers and the MPFR drivers, whose terms can lie anywhere in MPFR's
 * exponent range.  Read through doubles, a term beyond 1.8e308 became +Inf and
 * one below 2.2e-308 became 0: an overflowing term "rose without bound" and was
 * cut as noise, and Inf <= reltol * Inf then passed the convergence test, so
 * NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 30] stopped at the first
 * term past 1e308 and returned 6.6e25591 for 1.10e36300.  Every field holding
 * a magnitude below holds its log2 (-Inf for zero). */
typedef struct {
    double floor_t;   /* |t| at which the current floor-level run began (<0: none) */
    double dec_t;     /* |t| at which the current decayed run began (<0: none)     */
    double dec_env;   /* log2 of the largest term in the current decayed run      */
    double prev_m;    /* log2 magnitude of the previous term on this side         */
    double side_max;  /* log2 of the largest term so far on this side, seeded with
                       * the shared t = 0 term                                    */
    bool   armed;     /* a noise checkpoint is held by the caller                 */
    double arm_env;   /* log2 envelope of the decayed run that armed it           */
} DeTailGuard;

enum {
    DE_TAIL_ADD = 0,        /* keep the term, continue                           */
    DE_TAIL_ADD_STOP,       /* keep the term, stop this side (floor reached)     */
    DE_TAIL_CHECKPOINT,     /* save the partial sum, then keep the term          */
    DE_TAIL_ROLLBACK_STOP   /* noise confirmed: restore the checkpoint and stop  */
};

/* lcenter: log2 magnitude of the level's t = 0 term (-Inf if none). */
static inline void de_tail_guard_init(DeTailGuard* g, double lcenter) {
    g->floor_t = -1.0; g->dec_t = -1.0; g->dec_env = -INFINITY;
    g->prev_m = INFINITY; g->side_max = lcenter;
    g->armed = false; g->arm_env = -INFINITY;
}

/* log2 of a nonnegative double magnitude (-Inf for zero). */
static inline double de_log2_mag(double m) {
    return (m > 0.0) ? log2(m) : -INFINITY;
}

/* Classify the term of log2 magnitude lm at |t| = at, given the log2 of the
 * largest term seen so far this level (*lmax, updated here).  lfloor is the
 * log2 of the relative floor. */
static inline int de_tail_step_log2(DeTailGuard* g, double at, double lm,
                                    double* lmax, double lfloor) {
    static const double L_GROWTH  = 13.287712379549449;   /* log2(1e4)  */
    static const double L_DECAYED = -3.3219280948873622;  /* log2(0.1)  */
    if (g->armed && lm > L_GROWTH + g->arm_env)
        return DE_TAIL_ROLLBACK_STOP;
    int act = DE_TAIL_ADD;
    if (!g->armed && g->dec_t >= 0.0 && at - g->dec_t >= DE_TAIL_WINDOW
        && lm > 2.0 + g->prev_m) {
        g->armed = true;
        g->arm_env = g->dec_env;
        act = DE_TAIL_CHECKPOINT;
    }
    if (lm > *lmax) *lmax = lm;
    double mx = *lmax;
    g->prev_m = lm;
    /* Both cuts read "the tail of this side", so they wait until the side has
     * been significant: its largest term (seeded with the shared t = 0 term)
     * has reached the decayed band of the level maximum.  The level maximum
     * alone is shared with the other side, so a side whose integrand starts at
     * a deep zero and climbs (the order-2000 zero of (1 + e^(i phi))^2000 at
     * the centre of [0, 2 Pi]) sat below the other side's peak for a whole
     * window and was cut before its own bulk: half the integral was lost.  A
     * side that starts at or rises to a significant term is unaffected. */
    if (lm > g->side_max) g->side_max = lm;
    bool have_max = mx > -INFINITY && g->side_max > L_DECAYED + mx;
    if (have_max && lm <= L_DECAYED + mx) {
        if (g->dec_t < 0.0) { g->dec_t = at; g->dec_env = lm; }
        else if (lm > g->dec_env) g->dec_env = lm;
    } else {
        g->dec_t = -1.0;
    }
    if (have_max && lm <= lfloor + mx) {
        if (g->floor_t < 0.0) g->floor_t = at;
        if (at - g->floor_t >= DE_TAIL_WINDOW && act == DE_TAIL_ADD)
            return DE_TAIL_ADD_STOP;
    } else {
        g->floor_t = -1.0;
    }
    return act;
}

/* Linear-magnitude front end for the machine drivers: m and *maxterm are plain
 * magnitudes, floor_rel the relative floor. */
static inline int de_tail_step(DeTailGuard* g, double at, double m,
                               double* maxterm, double floor_rel) {
    double lmax = de_log2_mag(*maxterm);
    int act = de_tail_step_log2(g, at, de_log2_mag(m), &lmax, de_log2_mag(floor_rel));
    if (m > *maxterm) *maxterm = m;
    return act;
}

#ifdef USE_MPFR
/* log2 |v| over MPFR's whole exponent range (-Inf for zero). */
static inline double de_mpfr_log2_abs(const mpfr_t v) {
    if (mpfr_zero_p(v) || !mpfr_number_p(v)) return -INFINITY;
    long e;
    double d = mpfr_get_d_2exp(&e, v, MPFR_RNDN);       /* 0.5 <= |d| < 1 */
    return log2(fabs(d)) + (double)e;
}

/* log2(|re| + |im|), the gauge the MPFR drivers use for terms and level
 * differences: no overflow to Inf, no underflow to 0. */
static inline double de_mpfr_log2_l1(const mpfr_t re, const mpfr_t im) {
    double lr = de_mpfr_log2_abs(re), li = de_mpfr_log2_abs(im);
    double hi = lr > li ? lr : li, lo = lr > li ? li : lr;
    if (hi == -INFINITY) return -INFINITY;
    return hi + log2(1.0 + exp2(lo - hi));
}

/* Level-to-level convergence test of the MPFR drivers, scale-free: the
 * difference d of two successive estimates I is accepted when
 * |d| <= reltol |I|, or when |d| is at the rounding noise of the level's
 * largest term (log2 lmax_term, already multiplied by the step) at `bits`
 * precision, which is all an integral that cancels to ~0 can show. */
static inline bool de_mpfr_level_converged(const mpfr_t d_re, const mpfr_t d_im,
                                           const mpfr_t I_re, const mpfr_t I_im,
                                           double reltol, double lterm_scaled,
                                           long bits) {
    double ld = de_mpfr_log2_l1(d_re, d_im);
    if (ld == -INFINITY) return true;
    if (ld <= log2(reltol) + de_mpfr_log2_l1(I_re, I_im)) return true;
    return lterm_scaled > -INFINITY && ld <= lterm_scaled - (double)bits + 8.0;
}

/* The level difference as a double magnitude for the drivers' *abserr out
 * parameter (saturates at the double range; the convergence decision itself
 * never goes through it). */
static inline double de_mpfr_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}
#endif /* USE_MPFR */

#endif /* MATHILDA_DEQUAD_H */
