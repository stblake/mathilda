/*
 * oscint.c — oscillatory-integrand quadrature.  See oscint.h.
 *
 * The local half-period is estimated by counting sign changes of the dominant
 * (larger-amplitude) component over a window that is grown until it spans
 * several oscillations.  The range is then covered by panels about one
 * half-period wide, each integrated with the adaptive Gauss-Kronrod rule.  For
 * a finite range the panels are summed; for a half line the running partial
 * sums are extrapolated with Wynn's epsilon algorithm.
 */

#include "oscint.h"
#include "seqaccel.h"

#ifdef USE_MPFR
#  include "denint.h"   /* denint_tanhsinh_mpfr — full-precision lobe panels */
#endif

#include <math.h>
#include <stdlib.h>

/* Estimate the local half-period of f near a by counting sign changes of the
 * larger-amplitude component over a window grown until it holds >= 4 of them.
 * Also reports which component (0 = Re, 1 = Im) carries the oscillation. */
static bool osc_halfperiod_comp(GkSampleMachine f, void* ctx, double a, double b,
                                bool infinite, double* hp, int* comp_out) {
    double win = infinite ? 1.0 : fmin(b - a, fmax(1.0, fabs(a) + 1.0));
    double maxwin = infinite ? 1e6 : (b - a);
    for (int grow = 0; grow < 44 && win <= maxwin * 1.0001; grow++) {
        const int N = 65;
        int scRe = 0, scIm = 0;
        double pRe = 0, pIm = 0; bool have = false;
        double minRe = 1e300, maxRe = -1e300, minIm = 1e300, maxIm = -1e300;
        for (int i = 0; i < N; i++) {
            /* Sample at cell centres so a sample never lands exactly on an
             * endpoint: a singular endpoint (1/x at x=0) would otherwise raise a
             * spurious Power::infy and contribute a discarded non-finite value. */
            double x = a + win * ((double)i + 0.5) / (double)N;
            double _Complex v;
            if (!f(ctx, x, &v)) continue;
            double re = creal(v), im = cimag(v);
            if (re < minRe) minRe = re;
            if (re > maxRe) maxRe = re;
            if (im < minIm) minIm = im;
            if (im > maxIm) maxIm = im;
            if (have) {
                if ((pRe <= 0 && re > 0) || (pRe >= 0 && re < 0)) scRe++;
                if ((pIm <= 0 && im > 0) || (pIm >= 0 && im < 0)) scIm++;
            }
            pRe = re; pIm = im; have = true;
        }
        double ampRe = maxRe - minRe, ampIm = maxIm - minIm;
        int comp = (ampRe >= ampIm) ? 0 : 1;
        int sc = comp ? scIm : scRe;
        if (sc >= 4) { *hp = win / (double)sc; *comp_out = comp; return true; }
        win *= 2.0;
    }
    return false;
}

/* Locate the next sign change (zero) of component `comp` of f strictly after x.
 *
 * The local oscillation frequency can drift markedly from one panel to the next
 * — most sharply for an oscillatory endpoint singularity, where it grows without
 * bound.  A fixed marching step (sized from a single global half-period estimate)
 * then either crawls or, worse, leaps over whole lobes once the frequency rises.
 * Instead we probe outward from x with a step that starts at a small fraction of
 * the last zero-to-zero gap `scale` and grows geometrically until a sign change
 * is bracketed, then bisect.  Starting small guarantees a closer-than-expected
 * next zero is never skipped; the geometric growth keeps the probe count low
 * when the frequency instead falls. */
static bool osc_next_zero(GkSampleMachine f, void* ctx, int comp,
                          double x, double scale, double* z) {
    if (!(scale > 0.0)) return false;
    double _Complex v;
    /* First sample just inside the next lobe (the previous zero sits at x). */
    double h = 0.02 * scale;
    double pos = x + h;
    if (!f(ctx, pos, &v)) {
        bool ok = false;
        for (int k = 0; k < 8; k++) { h *= 1.7; pos = x + h; if (f(ctx, pos, &v)) { ok = true; break; } }
        if (!ok) return false;
    }
    double f0 = comp ? cimag(v) : creal(v);
    for (int i = 0; i < 100000; i++) {
        double nx = pos + h;
        double _Complex vn;
        if (!f(ctx, nx, &vn)) return false;
        double f1 = comp ? cimag(vn) : creal(vn);
        if ((f0 <= 0 && f1 >= 0) || (f0 >= 0 && f1 <= 0)) {
            double lo = pos, hi = nx, flo = f0;
            for (int bj = 0; bj < 60; bj++) {
                double m = 0.5 * (lo + hi);
                double _Complex vm;
                if (!f(ctx, m, &vm)) return false;
                double fm = comp ? cimag(vm) : creal(vm);
                if ((flo <= 0 && fm >= 0) || (flo >= 0 && fm <= 0)) hi = m;
                else { lo = m; flo = fm; }
            }
            *z = 0.5 * (lo + hi);
            return true;
        }
        pos = nx; f0 = f1;
        h *= 1.3;
    }
    return false;
}

/* Integrate one panel [x0,x1] with a lightly-adaptive Gauss-Kronrod rule. */
static bool osc_panel(GkSampleMachine f, void* ctx, double x0, double x1,
                      double reltol, double _Complex* out) {
    GkResult R = gk_integrate_machine(f, ctx, x0, x1, 0.0, reltol, 8, 0, false);
    if (R.status == GK_NONNUMERIC) return false;
    *out = R.value;
    return true;
}

/* Sum half-period panels of width `w` across the finite range [a,b].  Writes the
 * total and returns whether every panel was numeric and the range was covered. */
static bool osc_finite_sum(GkSampleMachine f, void* ctx, double a, double b,
                           double w, double panel_tol, int max_panels,
                           double _Complex* out) {
    double _Complex total = 0.0;
    double x = a;
    int np = 0;
    bool all_ok = true;
    while (x < b && np < max_panels) {
        double x1 = x + w;
        if (x1 > b) x1 = b;
        double _Complex p;
        if (osc_panel(f, ctx, x, x1, panel_tol, &p)) total += p;
        else all_ok = false;
        x = x1; np++;
    }
    *out = total;
    return all_ok && x >= b;
}

bool osc_integrate_machine(GkSampleMachine f, void* ctx, double a, double b,
                           bool infinite, double reltol, int max_panels,
                           double _Complex* result, double* abserr) {
    *result = 0.0; *abserr = INFINITY;
    if (reltol <= 0.0) reltol = 1e-10;
    double hp; int comp;
    if (!osc_halfperiod_comp(f, ctx, a, b, infinite, &hp, &comp) || !(hp > 0.0)) return false;
    double w = 0.9 * hp;
    double panel_tol = fmax(reltol * 0.01, 1e-13);

    if (infinite) {
        /* Panels run between successive zeros of the oscillation, so each
         * carries a single-signed lobe and the running partial sums alternate
         * cleanly — exactly the sequence Wynn's epsilon accelerates.  The local
         * lobe width (`scale`) is seeded from the half-period estimate and then
         * tracked from the last zero-to-zero gap, so the zero finder stays in
         * step even when the frequency drifts (e.g. an oscillatory singularity
         * whose lobes shrink without bound). */
        double scale = hp;
        double _Complex* S = malloc((size_t)max_panels * sizeof(double _Complex));
        if (!S) return false;
        int n = 0;
        double _Complex total = 0.0, best = 0.0;
        double err = INFINITY, maxS = 0.0;
        bool conv = false;
        double x = a;
        for (int k = 0; k < max_panels; k++) {
            double z;
            if (!osc_next_zero(f, ctx, comp, x, scale, &z)) break;
            double _Complex p;
            if (osc_panel(f, ctx, x, z, panel_tol, &p)) total += p;
            scale = z - x;
            x = z;
            S[n++] = total;
            if (cabs(total) > maxS) maxS = cabs(total);
            if (n >= 5) {
                int deg = (n - 1) / 2; if (deg > 6) deg = 6; if (deg < 1) deg = 1;
                double _Complex acc; double sstep;
                if (seqaccel_wynn_machine(S, n, deg, &acc, &sstep)
                    && isfinite(creal(acc)) && isfinite(cimag(acc)) && isfinite(sstep)
                    /* reject Wynn's singular-denominator blow-ups, and keep the
                     * best (smallest-error) estimate — late high-degree entries
                     * degrade as the tableau corner goes singular. */
                    && cabs(acc) <= 1e3 * maxS + 1.0 && sstep < err) {
                    best = acc; err = sstep;
                    if (sstep <= reltol * cabs(acc) + 1e-14) { conv = true; break; }
                }
            }
        }
        free(S);
        /* If the extrapolation never reached the tolerance, the raw partial sum
         * over many half-periods is itself a (slowly-)convergent estimate. */
        *result = conv ? best : total;
        *abserr = err;
        return conv;
    }

    /* Finite range: cover [a,b] with half-period panels and sum.  The sum is
     * computed at panel width w and again at w/2; their difference is a genuine
     * error estimate (the two agree once the panels resolve the oscillation, and
     * diverge when they do not — e.g. an oscillatory endpoint singularity, where
     * no fixed width resolves the ever-finer lobes).  Reporting that honest error
     * (rather than a blanket 0) keeps a poorly-resolved estimate from being
     * mistaken for a converged one by the caller. */
    double _Complex s1, s2;
    bool ok1 = osc_finite_sum(f, ctx, a, b, w,       panel_tol, max_panels, &s1);
    bool ok2 = osc_finite_sum(f, ctx, a, b, 0.5 * w, panel_tol, max_panels, &s2);
    *result = s2;
    *abserr = cabs(s2 - s1);
    return ok1 && ok2 && *abserr <= reltol * cabs(s2) + 1e-12;
}

#ifdef USE_MPFR
/* Magnitude of an MPFR complex value as a double (for the geometry-only gates). */
static double osc_mpfr_mag(const mpfr_t re, const mpfr_t im) {
    return hypot(mpfr_get_d(re, MPFR_RNDN), mpfr_get_d(im, MPFR_RNDN));
}

bool osc_integrate_mpfr(GkSampleMachine fgeo, void* geoctx,
                        DeQuadSampleMPFR fmp, void* mpctx,
                        double a, long bits, double reltol,
                        int max_panels, int panel_levels,
                        mpfr_t out_re, mpfr_t out_im, double* abserr) {
    mpfr_set_zero(out_re, 1); mpfr_set_zero(out_im, 1);
    *abserr = INFINITY;
    if (reltol <= 0.0) reltol = 1e-10;
    if (panel_levels < 8) panel_levels = 8;

    /* Locate the oscillation geometry with the cheap machine sampler. */
    double hp; int comp;
    if (!osc_halfperiod_comp(fgeo, geoctx, a, 0.0, true, &hp, &comp) || !(hp > 0.0))
        return false;

    /* Each lobe is integrated to (a little tighter than) the target tolerance so
     * per-panel error does not accumulate across the hundreds of cancelling
     * lobes.  Its endpoints are zeros of the oscillation, so a lobe is smooth
     * and single-signed — exactly what tanh-sinh resolves quickly. */
    double panel_tol = reltol;

    mpfr_t* Sr = malloc((size_t)max_panels * sizeof(mpfr_t));
    mpfr_t* Si = malloc((size_t)max_panels * sizeof(mpfr_t));
    if (!Sr || !Si) { free(Sr); free(Si); return false; }

    mpfr_t total_r, total_i, pr, pi, x0m, x1m, accr, acci, bestr, besti;
    mpfr_inits2((mpfr_prec_t)bits, total_r, total_i, pr, pi, x0m, x1m,
                accr, acci, bestr, besti, (mpfr_ptr)0);
    mpfr_set_zero(total_r, 1); mpfr_set_zero(total_i, 1);
    mpfr_set_zero(bestr, 1);   mpfr_set_zero(besti, 1);

    /* Wynn's epsilon extracts a bounded number of digits per column, so the
     * machine path's degree-6 cap (ample for ~15 digits) would floor the
     * extrapolation error far above a high-precision tolerance and never
     * converge — exhausting every panel.  In MPFR's wide working precision the
     * high-degree tableau corner is numerically safe, so scale the cap to the
     * requested digits (~one Wynn degree per four digits). */
    int dcap = (int)(-log10(reltol) / 4.0) + 6;
    if (dcap < 6) dcap = 6;
    if (dcap > 80) dcap = 80;

    int n = 0;
    double err = INFINITY, maxS = 0.0;
    bool conv = false;
    double scale = hp, x = a;
    for (int k = 0; k < max_panels; k++) {
        double z;
        if (!osc_next_zero(fgeo, geoctx, comp, x, scale, &z)) break;
        mpfr_set_d(x0m, x, MPFR_RNDN);
        mpfr_set_d(x1m, z, MPFR_RNDN);
        double perr;
        denint_tanhsinh_mpfr(fmp, mpctx, x0m, x1m, bits, panel_tol,
                             panel_levels, pr, pi, &perr);
        /* Accept the panel's best estimate whenever it is a finite number
         * (mirrors the machine osc_panel, which keeps a non-converged panel). */
        if (mpfr_number_p(pr) && mpfr_number_p(pi)) {
            mpfr_add(total_r, total_r, pr, MPFR_RNDN);
            mpfr_add(total_i, total_i, pi, MPFR_RNDN);
        }
        scale = z - x;
        x = z;

        mpfr_init2(Sr[n], (mpfr_prec_t)bits); mpfr_init2(Si[n], (mpfr_prec_t)bits);
        mpfr_set(Sr[n], total_r, MPFR_RNDN);  mpfr_set(Si[n], total_i, MPFR_RNDN);
        n++;
        double mag = osc_mpfr_mag(total_r, total_i);
        if (mag > maxS) maxS = mag;

        if (n >= 5) {
            int deg = (n - 1) / 2; if (deg > dcap) deg = dcap; if (deg < 1) deg = 1;
            double sstep; bool fin;
            if (seqaccel_wynn_mpfr((const mpfr_t*)Sr, (const mpfr_t*)Si, n, deg,
                                   bits, accr, acci, &sstep, &fin)
                && fin && isfinite(sstep)) {
                double amag = osc_mpfr_mag(accr, acci);
                /* Reject Wynn's singular-denominator blow-ups; keep the
                 * smallest-step estimate (late high-degree entries degrade as
                 * the tableau corner goes singular). */
                if (amag <= 1e3 * maxS + 1.0 && sstep < err) {
                    mpfr_set(bestr, accr, MPFR_RNDN);
                    mpfr_set(besti, acci, MPFR_RNDN);
                    err = sstep;
                    /* Precision-scaled gate (no fixed double-precision floor). */
                    if (sstep <= reltol * amag + 1e-300) { conv = true; break; }
                }
            }
        }
    }

    /* Fall back to the raw partial sum if extrapolation never met tolerance. */
    if (conv) { mpfr_set(out_re, bestr,   MPFR_RNDN); mpfr_set(out_im, besti,   MPFR_RNDN); }
    else      { mpfr_set(out_re, total_r, MPFR_RNDN); mpfr_set(out_im, total_i, MPFR_RNDN); }
    *abserr = err;

    for (int i = 0; i < n; i++) { mpfr_clear(Sr[i]); mpfr_clear(Si[i]); }
    free(Sr); free(Si);
    mpfr_clears(total_r, total_i, pr, pi, x0m, x1m,
                accr, acci, bestr, besti, (mpfr_ptr)0);
    return conv;
}
#endif /* USE_MPFR */
