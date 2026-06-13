/*
 * dequad.c — double-exponential (exp-sinh) half-line quadrature (see dequad.h)
 *
 * Transformation:  x(t) = a + exp((π/2) sinh t),  covering (a, ∞) as t runs
 * over (-∞, ∞), with weight x'(t) = (π/2) cosh t · exp((π/2) sinh t).  The
 * trapezoidal rule on t converges double-exponentially for integrands analytic
 * on the half line and decaying at +∞.  Each refinement halves the step h and
 * compares against the previous estimate; the tails are truncated where the
 * transformed integrand falls below a relative floor.
 */

#include "dequad.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEQ_TMAX_MACHINE 7.0      /* |t| beyond this overflows the transform   */
#define DEQ_REL_FLOOR    1e-17    /* tail cut: term < floor * running max       */

/* ------------------------------------------------------------------ *
 *  Machine precision                                                  *
 * ------------------------------------------------------------------ */

/* One transformed sample at parameter t. Returns false (skip/truncate) if the
 * transform overflows, x is non-finite, or the user sample fails. */
static bool deq_term_machine(DeQuadSampleMachine f, void* ctx, double a,
                             double t, double _Complex* term) {
    double w = (M_PI / 2.0) * sinh(t);
    if (w > 700.0) return false;               /* exp(w) would overflow        */
    double ex = exp(w);
    double x = a + ex;
    if (!isfinite(x)) return false;
    double weight = (M_PI / 2.0) * cosh(t) * ex;
    if (!isfinite(weight)) return false;
    double _Complex fv;
    if (!f(ctx, x, &fv)) return false;
    double _Complex tm = fv * weight;
    if (!isfinite(creal(tm)) || !isfinite(cimag(tm))) return false;
    *term = tm;
    return true;
}

bool dequad_halfline_machine(DeQuadSampleMachine f, void* ctx, double a,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr) {
    if (reltol <= 0.0) reltol = 1e-13;
    if (max_levels < 1) max_levels = 1;
    double h = 1.0;                            /* level-0 step                  */
    double _Complex prev = 0.0, best = 0.0;
    bool have_prev = false;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++, h *= 0.5) {
        double _Complex sum = 0.0, c0;
        if (deq_term_machine(f, ctx, a, 0.0, &c0)) sum = c0;
        double maxterm = cabs(sum);

        /* Each tail: add terms until they fall below the relative floor, or —
         * crucially for integrands that hit a roundoff noise floor while the
         * exp-sinh weight keeps growing — until the term, having decayed into
         * the tail, starts *rising* again (the noise-amplification signature). */
        for (int side = -1; side <= 1; side += 2) {
            double prev_m = INFINITY;
            bool decaying = false;
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > DEQ_TMAX_MACHINE) break;
                double _Complex tm;
                if (!deq_term_machine(f, ctx, a, t, &tm)) break;
                double m = cabs(tm);
                if (decaying && k > 3 && m > 4.0 * prev_m) break; /* tail noise blow-up */
                sum += tm;
                if (m > maxterm) maxterm = m;
                if (m < 0.1 * maxterm) decaying = true;          /* past the peak    */
                if (k > 2 && m <= DEQ_REL_FLOOR * maxterm) break;
                prev_m = m;
            }
        }

        double _Complex Ival = h * sum;
        best = Ival;
        if (have_prev) {
            double diff = cabs(Ival - prev);
            *abserr = diff;
            if (diff <= reltol * cabs(Ival) + 1e-300) { *result = Ival; return true; }
        }
        prev = Ival;
        have_prev = true;
    }
    *result = best;
    return false;
}

/* ------------------------------------------------------------------ *
 *  MPFR precision                                                     *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR

/* Scratch carried across samples so each transformed point reuses buffers. */
typedef struct {
    mpfr_t halfpi, w, ex, cosh_t, weight, x, fr, fi, t;
} DeqMpfrScratch;

static void deq_scratch_init(DeqMpfrScratch* s, mpfr_prec_t p) {
    mpfr_inits2(p, s->halfpi, s->w, s->ex, s->cosh_t, s->weight,
                s->x, s->fr, s->fi, s->t, (mpfr_ptr)0);
    mpfr_const_pi(s->halfpi, MPFR_RNDN);
    mpfr_div_2ui(s->halfpi, s->halfpi, 1, MPFR_RNDN);   /* π/2 */
}

static void deq_scratch_clear(DeqMpfrScratch* s) {
    mpfr_clears(s->halfpi, s->w, s->ex, s->cosh_t, s->weight,
                s->x, s->fr, s->fi, s->t, (mpfr_ptr)0);
}

/* Transformed sample at parameter `tv` (a double), accumulating Re/Im of the
 * weighted value into (out_re,out_im) via add. Returns false to truncate. */
static bool deq_term_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a,
                          DeqMpfrScratch* s, double tv,
                          mpfr_t term_re, mpfr_t term_im) {
    mpfr_set_d(s->t, tv, MPFR_RNDN);
    mpfr_sinh(s->w, s->t, MPFR_RNDN);
    mpfr_mul(s->w, s->w, s->halfpi, MPFR_RNDN);        /* w = (π/2) sinh t      */
    mpfr_exp(s->ex, s->w, MPFR_RNDN);                  /* exp(w)                */
    if (!mpfr_number_p(s->ex)) return false;
    mpfr_add(s->x, a, s->ex, MPFR_RNDN);               /* x = a + exp(w)        */
    if (!mpfr_number_p(s->x)) return false;
    mpfr_cosh(s->cosh_t, s->t, MPFR_RNDN);
    mpfr_mul(s->weight, s->halfpi, s->cosh_t, MPFR_RNDN);
    mpfr_mul(s->weight, s->weight, s->ex, MPFR_RNDN);  /* (π/2) cosh t exp(w)   */
    if (!mpfr_number_p(s->weight)) return false;
    if (!f(ctx, s->x, s->fr, s->fi)) return false;
    if (!mpfr_number_p(s->fr) || !mpfr_number_p(s->fi)) return false;
    mpfr_mul(term_re, s->fr, s->weight, MPFR_RNDN);
    mpfr_mul(term_im, s->fi, s->weight, MPFR_RNDN);
    return true;
}

/* |re| + |im| as a double, for the relative tail/convergence gauges. */
static double deq_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}

bool dequad_halfline_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr) {
    if (reltol <= 0.0) reltol = 1e-13;
    if (max_levels < 1) max_levels = 1;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    DeqMpfrScratch sc; deq_scratch_init(&sc, p);

    /* t range grows with precision; the relative-floor break does the real work
     * but a slowly-decaying algebraic tail (∝ x^{-1-ε}, ε→0, as in ζ(s) for s
     * near 1) only crosses the floor at large t, so the cap must be generous.
     * MPFR's wide exponent range tolerates the huge exp((π/2) sinh t). */
    double tmax = 6.5 + 1.1 * log((double)bits);
    double floor_rel = ldexp(1.0, -(int)bits - 4);

    mpfr_t sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, diff_re, diff_im;
    mpfr_inits2(p, sum_re, sum_im, I_re, I_im, prev_re, prev_im,
                tr, ti, diff_re, diff_im, (mpfr_ptr)0);

    double h = 1.0;
    bool have_prev = false, converged = false;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++, h *= 0.5) {
        mpfr_set_ui(sum_re, 0, MPFR_RNDN);
        mpfr_set_ui(sum_im, 0, MPFR_RNDN);
        double maxterm = 0.0;
        if (deq_term_mpfr(f, ctx, a, &sc, 0.0, tr, ti)) {
            mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
            mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
            maxterm = deq_l1_d(tr, ti);
        }
        for (int side = -1; side <= 1; side += 2) {
            double prev_m = INFINITY;
            bool decaying = false;
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > tmax) break;
                if (!deq_term_mpfr(f, ctx, a, &sc, t, tr, ti)) break;
                double m = deq_l1_d(tr, ti);
                if (decaying && k > 3 && m > 4.0 * prev_m) break; /* tail noise blow-up */
                mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
                mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
                if (m > maxterm) maxterm = m;
                if (m < 0.1 * maxterm) decaying = true;          /* past the peak    */
                if (k > 2 && m <= floor_rel * maxterm) break;
                prev_m = m;
            }
        }
        mpfr_mul_d(I_re, sum_re, h, MPFR_RNDN);
        mpfr_mul_d(I_im, sum_im, h, MPFR_RNDN);

        mpfr_set(out_re, I_re, MPFR_RNDN);
        mpfr_set(out_im, I_im, MPFR_RNDN);
        if (have_prev) {
            mpfr_sub(diff_re, I_re, prev_re, MPFR_RNDN);
            mpfr_sub(diff_im, I_im, prev_im, MPFR_RNDN);
            double diff = deq_l1_d(diff_re, diff_im);
            double mag = deq_l1_d(I_re, I_im);
            *abserr = diff;
            if (diff <= reltol * mag + 1e-300) { converged = true; break; }
        }
        mpfr_set(prev_re, I_re, MPFR_RNDN);
        mpfr_set(prev_im, I_im, MPFR_RNDN);
        have_prev = true;
    }

    deq_scratch_clear(&sc);
    mpfr_clears(sum_re, sum_im, I_re, I_im, prev_re, prev_im,
                tr, ti, diff_re, diff_im, (mpfr_ptr)0);
    return converged;
}
#endif /* USE_MPFR */
