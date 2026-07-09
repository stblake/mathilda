/*
 * ncrule.c — equally-spaced composite quadrature rules (see ncrule.h)
 *
 * Each rule is a composite Newton–Cotes-family estimate over n equal panels,
 * refined by repeatedly doubling n.  Trapezoidal and Newton–Cotes optionally
 * drive a Richardson (Romberg) table over the halvings; the Riemann rectangle
 * rules are returned raw.  Only two rows of the Romberg table are kept at a time
 * (level j needs only level j−1), which bounds the MPFR variant's allocation.
 */
#include "ncrule.h"

#include <math.h>

/* Closed Newton–Cotes weights w[p][i] for one panel of p subintervals (p+1
 * points): ∫ ≈ h · Σ w[p][i] f_i, with h the subinterval width.  p = 1..4.  The
 * weights are rationals NC_NUM[p][i] / NC_DEN[p]; the double table is for the
 * machine kernel, while the MPFR kernel uses the exact integer numerator and a
 * single final division by the denominator so a high-precision result is not
 * capped at machine epsilon by a double-rounded 4/3, 1/45, … weight. */
static const double NC_W[5][5] = {
    { 0.0,      0.0,      0.0,      0.0,      0.0      },  /* unused        */
    { 1.0/2.0,  1.0/2.0,  0.0,      0.0,      0.0      },  /* trapezoid     */
    { 1.0/3.0,  4.0/3.0,  1.0/3.0,  0.0,      0.0      },  /* Simpson       */
    { 3.0/8.0,  9.0/8.0,  9.0/8.0,  3.0/8.0,  0.0      },  /* Simpson 3/8   */
    {14.0/45.0,64.0/45.0,24.0/45.0,64.0/45.0,14.0/45.0}   /* Boole         */
};

static const long NC_NUM[5][5] = {
    { 0,  0,  0,  0,  0},
    { 1,  1,  0,  0,  0},   /* trapezoid : [1,1]/2          */
    { 1,  4,  1,  0,  0},   /* Simpson   : [1,4,1]/3        */
    { 3,  9,  9,  3,  0},   /* Simpson3/8: [3,9,9,3]/8      */
    {14, 64, 24, 64, 14}    /* Boole     : [14,64,24,64,14]/45 */
};
static const long NC_DEN[5] = { 1, 2, 3, 8, 45 };

#define NCR_MAX_LEVELS 60

/* Leading error order q for Richardson extrapolation: the symmetric rules have
 * error expansions in even powers of h beyond the leading term, so column m
 * cancels the h^{q+2(m-1)} term.  Closed Newton–Cotes of p subintervals is
 * order p+2 for even p (Simpson, Boole) and p+1 for odd p (trapezoid, 3/8). */
static int ncr_order(NcrRuleKind kind, int p) {
    switch (kind) {
        case NCR_TRAPEZOIDAL:      return 2;
        case NCR_RIEMANN_MIDPOINT: return 2;
        case NCR_RIEMANN_LEFT:
        case NCR_RIEMANN_RIGHT:    return 1;
        case NCR_NEWTONCOTES:      return (p % 2 == 0) ? p + 2 : p + 1;
    }
    return 2;
}

/* Number of panels (subintervals) at refinement level j. */
static long ncr_panels(NcrRuleKind kind, int p, int j) {
    return (kind == NCR_NEWTONCOTES) ? ((long)p << j) : ((long)1 << j);
}

/* Integrand evaluations a level's estimate will perform (for budgeting). */
static long ncr_level_evals(NcrRuleKind kind, int p, long n) {
    if (kind == NCR_NEWTONCOTES) return (n / p) * (p + 1);
    if (kind == NCR_TRAPEZOIDAL) return n + 1;
    return n;                                   /* Riemann */
}

/* ------------------------------------------------------------------ *
 *  Machine precision                                                  *
 * ------------------------------------------------------------------ */

/* Composite estimate over n panels.  Endpoints are formed from the panel
 * fraction so a + 0 = a and a + 1·(b−a) = b exactly. */
static bool ncr_estimate(GkSampleMachine f, void* ctx, double a, double b,
                         NcrRuleKind kind, int p, long n, double _Complex* out) {
    double span = b - a;
    double _Complex acc = 0.0, v;

    if (kind == NCR_NEWTONCOTES) {
        long blocks = n / p;
        for (long blk = 0; blk < blocks; blk++) {
            for (int i = 0; i <= p; i++) {
                long k = blk * p + i;
                double x = a + ((double)k / (double)n) * span;
                if (!f(ctx, x, &v)) return false;
                acc += NC_W[p][i] * v;
            }
        }
        *out = (span / (double)n) * acc;
        return true;
    }

    if (kind == NCR_TRAPEZOIDAL) {
        if (!f(ctx, a, &v)) return false;
        acc += 0.5 * v;
        if (!f(ctx, b, &v)) return false;
        acc += 0.5 * v;
        for (long k = 1; k < n; k++) {
            double x = a + ((double)k / (double)n) * span;
            if (!f(ctx, x, &v)) return false;
            acc += v;
        }
        *out = (span / (double)n) * acc;
        return true;
    }

    /* Riemann left / right / midpoint. */
    for (long k = 0; k < n; k++) {
        double frac;
        if      (kind == NCR_RIEMANN_LEFT)  frac = (double)k / (double)n;
        else if (kind == NCR_RIEMANN_RIGHT) frac = (double)(k + 1) / (double)n;
        else                                frac = ((double)k + 0.5) / (double)n;
        double x = a + frac * span;
        if (!f(ctx, x, &v)) return false;
        acc += v;
    }
    *out = (span / (double)n) * acc;
    return true;
}

bool ncr_integrate_machine(GkSampleMachine f, void* ctx, double a, double b,
                           NcrRuleKind kind, int nc_points, bool romberg,
                           double reltol, double abstol,
                           int max_levels, long max_eval,
                           double _Complex* result, double* abserr) {
    *abserr = 0.0;
    if (a == b) { *result = 0.0; return true; }

    int p = nc_points - 1;
    if (kind == NCR_NEWTONCOTES) { if (p < 1) p = 1; if (p > 4) p = 4; }
    if (kind != NCR_TRAPEZOIDAL && kind != NCR_NEWTONCOTES) romberg = false;
    int q = ncr_order(kind, p);
    if (max_levels < 1) max_levels = 1;
    if (max_levels > NCR_MAX_LEVELS) max_levels = NCR_MAX_LEVELS;

    double _Complex prev[NCR_MAX_LEVELS + 1], cur[NCR_MAX_LEVELS + 1];
    double _Complex best = 0.0, best_prev = 0.0;
    bool have_prev = false;
    long evals = 0;

    for (int j = 0; j <= max_levels; j++) {
        long n = ncr_panels(kind, p, j);
        long lvl_evals = ncr_level_evals(kind, p, n);
        if (have_prev && max_eval > 0 && evals + lvl_evals > max_eval) break;

        double _Complex E;
        if (!ncr_estimate(f, ctx, a, b, kind, p, n, &E)) {
            if (!have_prev) { *result = NAN + NAN * I; return false; }  /* none */
            break;                              /* keep best so far    */
        }
        evals += lvl_evals;

        cur[0] = E;
        for (int m = 1; m <= j; m++) {
            double factor = pow(2.0, (double)(q + 2 * (m - 1)));
            cur[m] = cur[m - 1] + (cur[m - 1] - prev[m - 1]) / (factor - 1.0);
        }
        best = romberg ? cur[j] : cur[0];

        if (have_prev) {
            double err = cabs(best - best_prev);
            *abserr = err;
            if (err <= abstol || err <= reltol * cabs(best)) {
                *result = best;
                return true;
            }
        }
        for (int m = 0; m <= j; m++) prev[m] = cur[m];
        best_prev = best;
        have_prev = true;
    }

    *result = best;
    return false;                               /* best estimate, not converged */
}

/* ------------------------------------------------------------------ *
 *  Arbitrary precision (MPFR)                                         *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR

/* Abscissa for node fraction num/den (in lowest practical terms): x = a + span·
 * num/den, with the exact endpoints substituted to avoid drift. */
static void ncr_node_mpfr(mpfr_t x, const mpfr_t a, const mpfr_t b,
                          const mpfr_t span, long num, long den, mpfr_t t) {
    if (num == 0)        { mpfr_set(x, a, MPFR_RNDN); return; }
    if (num == den)      { mpfr_set(x, b, MPFR_RNDN); return; }
    mpfr_mul_si(t, span, num, MPFR_RNDN);
    mpfr_div_si(t, t, den, MPFR_RNDN);
    mpfr_add(x, a, t, MPFR_RNDN);
}

static bool ncr_estimate_mpfr(DeQuadSampleMPFR f, void* ctx,
                              const mpfr_t a, const mpfr_t b, const mpfr_t span,
                              NcrRuleKind kind, int p, long n, long bits,
                              mpfr_t out_re, mpfr_t out_im) {
    mpfr_t accr, acci, x, vr, vi, t;
    mpfr_inits2((mpfr_prec_t)bits, accr, acci, x, vr, vi, t, (mpfr_ptr)0);
    mpfr_set_zero(accr, 1);
    mpfr_set_zero(acci, 1);
    bool ok = true;

    if (kind == NCR_NEWTONCOTES) {
        long blocks = n / p;
        for (long blk = 0; blk < blocks && ok; blk++) {
            for (int i = 0; i <= p; i++) {
                long k = blk * p + i;
                ncr_node_mpfr(x, a, b, span, k, n, t);
                if (!f(ctx, x, vr, vi)) { ok = false; break; }
                long w = NC_NUM[p][i];            /* exact integer weight */
                mpfr_mul_si(t, vr, w, MPFR_RNDN); mpfr_add(accr, accr, t, MPFR_RNDN);
                mpfr_mul_si(t, vi, w, MPFR_RNDN); mpfr_add(acci, acci, t, MPFR_RNDN);
            }
        }
        if (ok) {                                  /* fold the shared denominator */
            mpfr_div_si(accr, accr, NC_DEN[p], MPFR_RNDN);
            mpfr_div_si(acci, acci, NC_DEN[p], MPFR_RNDN);
        }
    } else if (kind == NCR_TRAPEZOIDAL) {
        if (f(ctx, a, vr, vi)) {
            mpfr_mul_d(t, vr, 0.5, MPFR_RNDN); mpfr_add(accr, accr, t, MPFR_RNDN);
            mpfr_mul_d(t, vi, 0.5, MPFR_RNDN); mpfr_add(acci, acci, t, MPFR_RNDN);
        } else ok = false;
        if (ok) {
            if (f(ctx, b, vr, vi)) {
                mpfr_mul_d(t, vr, 0.5, MPFR_RNDN); mpfr_add(accr, accr, t, MPFR_RNDN);
                mpfr_mul_d(t, vi, 0.5, MPFR_RNDN); mpfr_add(acci, acci, t, MPFR_RNDN);
            } else ok = false;
        }
        for (long k = 1; k < n && ok; k++) {
            ncr_node_mpfr(x, a, b, span, k, n, t);
            if (!f(ctx, x, vr, vi)) { ok = false; break; }
            mpfr_add(accr, accr, vr, MPFR_RNDN);
            mpfr_add(acci, acci, vi, MPFR_RNDN);
        }
    } else {
        for (long k = 0; k < n && ok; k++) {
            long num, den;
            if      (kind == NCR_RIEMANN_LEFT)  { num = k;         den = n;     }
            else if (kind == NCR_RIEMANN_RIGHT) { num = k + 1;     den = n;     }
            else                                { num = 2 * k + 1; den = 2 * n; }
            ncr_node_mpfr(x, a, b, span, num, den, t);
            if (!f(ctx, x, vr, vi)) { ok = false; break; }
            mpfr_add(accr, accr, vr, MPFR_RNDN);
            mpfr_add(acci, acci, vi, MPFR_RNDN);
        }
    }

    if (ok) {
        mpfr_div_si(t, span, n, MPFR_RNDN);     /* h = span / n */
        mpfr_mul(out_re, accr, t, MPFR_RNDN);
        mpfr_mul(out_im, acci, t, MPFR_RNDN);
    }
    mpfr_clears(accr, acci, x, vr, vi, t, (mpfr_ptr)0);
    return ok;
}

bool ncr_integrate_mpfr(DeQuadSampleMPFR f, void* ctx,
                        const mpfr_t a, const mpfr_t b,
                        NcrRuleKind kind, int nc_points, bool romberg, long bits,
                        double reltol, double abstol,
                        int max_levels, long max_eval,
                        mpfr_t out_re, mpfr_t out_im, double* abserr) {
    *abserr = 0.0;
    mpfr_set_zero(out_re, 1);
    mpfr_set_zero(out_im, 1);
    if (mpfr_equal_p(a, b)) return true;

    int p = nc_points - 1;
    if (kind == NCR_NEWTONCOTES) { if (p < 1) p = 1; if (p > 4) p = 4; }
    if (kind != NCR_TRAPEZOIDAL && kind != NCR_NEWTONCOTES) romberg = false;
    int q = ncr_order(kind, p);
    if (max_levels < 1) max_levels = 1;
    if (max_levels > NCR_MAX_LEVELS) max_levels = NCR_MAX_LEVELS;

    int L = max_levels;
    mpfr_t prev_re[NCR_MAX_LEVELS + 1], prev_im[NCR_MAX_LEVELS + 1];
    mpfr_t cur_re[NCR_MAX_LEVELS + 1],  cur_im[NCR_MAX_LEVELS + 1];
    mpfr_t span, best_re, best_im, bp_re, bp_im, d_re, d_im, mag;
    for (int m = 0; m <= L; m++) {
        mpfr_init2(prev_re[m], (mpfr_prec_t)bits); mpfr_init2(prev_im[m], (mpfr_prec_t)bits);
        mpfr_init2(cur_re[m],  (mpfr_prec_t)bits); mpfr_init2(cur_im[m],  (mpfr_prec_t)bits);
    }
    mpfr_inits2((mpfr_prec_t)bits, span, best_re, best_im, bp_re, bp_im,
                d_re, d_im, mag, (mpfr_ptr)0);
    mpfr_sub(span, b, a, MPFR_RNDN);

    bool have_prev = false, converged = false;
    long evals = 0;

    for (int j = 0; j <= L; j++) {
        long n = ncr_panels(kind, p, j);
        long lvl_evals = ncr_level_evals(kind, p, n);
        if (have_prev && max_eval > 0 && evals + lvl_evals > max_eval) break;

        if (!ncr_estimate_mpfr(f, ctx, a, b, span, kind, p, n, bits,
                               cur_re[0], cur_im[0])) {
            if (!have_prev) {                   /* non-numeric: signal NaN    */
                mpfr_set_nan(out_re); mpfr_set_nan(out_im);
                goto done;
            }
            break;                              /* keep best so far           */
        }
        evals += lvl_evals;

        for (int m = 1; m <= j; m++) {
            double factor = pow(2.0, (double)(q + 2 * (m - 1)));
            mpfr_sub(d_re, cur_re[m - 1], prev_re[m - 1], MPFR_RNDN);
            mpfr_div_d(d_re, d_re, factor - 1.0, MPFR_RNDN);
            mpfr_add(cur_re[m], cur_re[m - 1], d_re, MPFR_RNDN);
            mpfr_sub(d_im, cur_im[m - 1], prev_im[m - 1], MPFR_RNDN);
            mpfr_div_d(d_im, d_im, factor - 1.0, MPFR_RNDN);
            mpfr_add(cur_im[m], cur_im[m - 1], d_im, MPFR_RNDN);
        }
        int bcol = romberg ? j : 0;
        mpfr_set(best_re, cur_re[bcol], MPFR_RNDN);
        mpfr_set(best_im, cur_im[bcol], MPFR_RNDN);

        if (have_prev) {
            mpfr_sub(d_re, best_re, bp_re, MPFR_RNDN);
            mpfr_sub(d_im, best_im, bp_im, MPFR_RNDN);
            mpfr_hypot(d_re, d_re, d_im, MPFR_RNDN);
            mpfr_hypot(mag, best_re, best_im, MPFR_RNDN);
            double err = mpfr_get_d(d_re, MPFR_RNDN);
            double bmag = mpfr_get_d(mag, MPFR_RNDN);
            *abserr = err;
            if (err <= abstol || err <= reltol * bmag) { converged = true; break; }
        }
        for (int m = 0; m <= j; m++) {
            mpfr_set(prev_re[m], cur_re[m], MPFR_RNDN);
            mpfr_set(prev_im[m], cur_im[m], MPFR_RNDN);
        }
        mpfr_set(bp_re, best_re, MPFR_RNDN);
        mpfr_set(bp_im, best_im, MPFR_RNDN);
        have_prev = true;
    }

    mpfr_set(out_re, best_re, MPFR_RNDN);
    mpfr_set(out_im, best_im, MPFR_RNDN);

done:
    for (int m = 0; m <= L; m++) {
        mpfr_clear(prev_re[m]); mpfr_clear(prev_im[m]);
        mpfr_clear(cur_re[m]);  mpfr_clear(cur_im[m]);
    }
    mpfr_clears(span, best_re, best_im, bp_re, bp_im, d_re, d_im, mag, (mpfr_ptr)0);
    return converged;
}

#endif /* USE_MPFR */
