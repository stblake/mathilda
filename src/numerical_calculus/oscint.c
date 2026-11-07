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

#include <float.h>
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
                /* at double resolution: further halvings cost a sample each
                 * and move nothing (60 a zero was a third of the cost of a
                 * machine lobe sum) */
                if (hi - lo <= 4.0 * DBL_EPSILON * fmax(fabs(lo), fabs(hi))) break;
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

/* The first lobe of a half line starts at the range's own endpoint a, which
 * may carry an integrable singularity of its own: the exponential endpoint map
 * of NIntegrate sends the OTHER end of a finite range to a, so x^a Sin[a/(2-x)]
 * on [0, 2] puts the branch point of x^0.14 at the start of the first lobe.
 * Eight bisections of the plain rule leave ~(1/256)^1.14 of that lobe's
 * singular part unresolved (4e-3 relative, silently, for a = 0.14), so the
 * first lobe gets the full adaptive rule with QAGS extrapolation, which
 * converges on an algebraic endpoint singularity and costs one panel more than
 * the plain rule on a smooth lobe. */
static bool osc_panel_first(GkSampleMachine f, void* ctx, double x0, double x1,
                            double reltol, double _Complex* out) {
    GkResult R = gk_integrate_machine(f, ctx, x0, x1, 0.0, reltol, 200, 0, true);
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
        bool conv = false, have_cand = false;
        double _Complex cand = 0.0;
        double x = a;
        for (int k = 0; k < max_panels; k++) {
            double z;
            if (!osc_next_zero(f, ctx, comp, x, scale, &z)) break;
            double _Complex p;
            if ((k == 0 ? osc_panel_first : osc_panel)(f, ctx, x, z, panel_tol, &p))
                total += p;
            scale = z - x;
            x = z;
            S[n++] = total;
            if (cabs(total) > maxS) maxS = cabs(total);
            if (n >= 5) {
                int deg = (n - 1) / 2; if (deg > 6) deg = 6; if (deg < 1) deg = 1;
                double _Complex acc; double sstep;
                bool valid = seqaccel_wynn_machine(S, n, deg, &acc, &sstep)
                    && isfinite(creal(acc)) && isfinite(cimag(acc)) && isfinite(sstep)
                    /* reject Wynn's singular-denominator blow-ups */
                    && cabs(acc) <= 1e3 * maxS + 1.0;
                /* A step at the tolerance is only a CANDIDATE: on a few lobes
                 * the tableau's last two entries can agree by coincidence (the
                 * lobe sum of x^3.14 Sin[3.14/(2-x)] jumped from a 2e-2 step to
                 * 3e-12 at seven lobes, 7e-5 away from the limit, and stopped
                 * there silently).  The next lobe's extrapolation must
                 * reproduce it to the tolerance before the sum is declared
                 * converged; a genuine geometric convergence costs one lobe. */
                if (valid && have_cand) {
                    double tol = reltol * cabs(acc) + 1e-14;
                    double dif = cabs(acc - cand);
                    if (dif <= tol) {
                        if (sstep < err) { best = acc; err = sstep; }
                        if (dif > err) err = dif;
                        conv = true; break;
                    }
                    have_cand = false;
                }
                /* keep the best (smallest-error) estimate -- late high-degree
                 * entries degrade as the tableau corner goes singular */
                if (valid && sstep < err) {
                    best = acc; err = sstep;
                    if (sstep <= reltol * cabs(acc) + 1e-14) { cand = acc; have_cand = true; }
                }
            }
        }
        free(S);
        /* The best extrapolation, converged or not (*abserr is its step); the
         * raw partial sum over the half-periods only when none was accepted. */
        *result = isfinite(err) ? best : total;
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

/* Refine the double zero z (in/out) of component `comp` of f to the working
 * precision by the secant method in MPFR.  The partial sum up to a zero that
 * is off by d differs from the true one by ~f'(z) d^2 / 2: with d at double
 * precision that is ~1e-32 f', a floor under every partial sum that Wynn's
 * epsilon then amplifies (the 50-digit lobe sum of Cos[t E^t] stalled at
 * 33 digits).  Three or four MPFR samples a zero; `half` bounds how far the
 * refined zero may move (half the lobe), and on any failure z keeps its
 * double value (the lobes still tile the range exactly). */
static void osc_refine_zero_mpfr(DeQuadSampleMPFR f, void* ctx, int comp,
                                 mpfr_t z, double half, long bits) {
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_t z0, z1, f0, f1, vr, vi, d, t, zb, fb;
    mpfr_inits2(p, z0, z1, f0, f1, vr, vi, d, t, zb, fb, (mpfr_ptr)0);
    mpfr_set(z0, z, MPFR_RNDN);
    double h = fabs(mpfr_get_d(z, MPFR_RNDN)) * 1e-13 + fabs(half) * 1e-12;
    bool ok = h > 0.0 && f(ctx, z0, vr, vi);
    if (ok) mpfr_set(f0, comp ? vi : vr, MPFR_RNDN);
    /* the best point so far, by |f|: the iteration stops at the noise floor
     * of the samples (an integrand with a 60-digit parameter carries fewer
     * bits than the working precision), where a further secant step is
     * garbage -- the best iterate is kept, not the double zero */
    bool have_best = ok;
    if (ok) { mpfr_set(zb, z0, MPFR_RNDN); mpfr_abs(fb, f0, MPFR_RNDN); }
    mpfr_add_d(z1, z0, h, MPFR_RNDN);
    ok = ok && f(ctx, z1, vr, vi);
    if (ok) mpfr_set(f1, comp ? vi : vr, MPFR_RNDN);
    double lstop = -(double)bits * 0.75;   /* d^2 below the working precision */
    for (int it = 0; ok && it < 8; it++) {
        mpfr_abs(t, f1, MPFR_RNDN);
        if (mpfr_less_p(t, fb)) { mpfr_set(zb, z1, MPFR_RNDN); mpfr_set(fb, t, MPFR_RNDN); }
        mpfr_sub(t, f1, f0, MPFR_RNDN);
        if (mpfr_zero_p(t) || mpfr_zero_p(f1)) break;
        mpfr_sub(d, z1, z0, MPFR_RNDN);
        mpfr_mul(d, d, f1, MPFR_RNDN);
        mpfr_div(d, d, t, MPFR_RNDN);              /* step = f1 (z1-z0)/(f1-f0) */
        mpfr_set(z0, z1, MPFR_RNDN); mpfr_set(f0, f1, MPFR_RNDN);
        mpfr_sub(z1, z1, d, MPFR_RNDN);
        mpfr_sub(t, z1, z, MPFR_RNDN);
        if (!mpfr_number_p(z1) || fabs(mpfr_get_d(t, MPFR_RNDN)) > fabs(half)) break;
        double ld = mpfr_zero_p(d) ? -INFINITY
                  : (double)mpfr_get_exp(d) - log2(fabs(half) + 1e-300);
        if (ld < lstop) {                          /* converged: z1 is the zero */
            mpfr_set(zb, z1, MPFR_RNDN); mpfr_set_zero(fb, 1);
            break;
        }
        if (!f(ctx, z1, vr, vi)) break;
        mpfr_set(f1, comp ? vi : vr, MPFR_RNDN);
    }
    if (have_best && mpfr_number_p(zb)) mpfr_set(z, zb, MPFR_RNDN);
    mpfr_clears(z0, z1, f0, f1, vr, vi, d, t, zb, fb, (mpfr_ptr)0);
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

    /* Each lobe is integrated to nearly the working precision, not just the
     * target tolerance: Wynn's high-degree corner amplifies the noise of the
     * partial sums by many orders (lobes good to the 1e-50 goal floored the
     * 50-digit extrapolation of Cos[t E^t] at 1e-35), and the guard bits
     * above the target are there to absorb it.  A lobe runs between zeros of
     * the oscillation, so it is analytic and single-signed: Clenshaw-Curtis
     * resolves it in 33-129 nested points, where tanh-sinh (built for
     * endpoint singularities) pays several hundred; tanh-sinh stays the rule
     * for a lobe that is not smooth (a singular first lobe at a). */
    double panel_tol = ldexp(1.0, -(int)(bits - 24));
    if (panel_tol > reltol) panel_tol = reltol;
    if (!(panel_tol > 0.0)) panel_tol = reltol;
    CcCache* cc = ccquad_cache_new(bits);

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
     * requested digits.  A degree buys about 2.5 digits on a lobe sum whose
     * terms decay like 1/n (an oscillatory endpoint singularity spread onto a
     * half line, or Sin[x]/x): one degree per two digits, so the cap is never
     * what stops a 50-digit extrapolation (one per four floored the 50-digit
     * lobe sum of Cos[t E^t] near 1e-46). */
    int dcap = (int)(-log10(reltol) / 2.0) + 6;
    if (dcap < 6) dcap = 6;
    if (dcap > 80) dcap = 80;

    int n = 0;
    double err = INFINITY, maxS = 0.0;
    bool conv = false;
    double scale = hp, x = a;
    mpfr_set_d(x0m, a, MPFR_RNDN);
    for (int k = 0; k < max_panels; k++) {
        double z;
        if (!osc_next_zero(fgeo, geoctx, comp, x, scale, &z)) break;
        /* x0m holds the previous (refined) zero, so the lobes tile exactly */
        mpfr_set_d(x1m, z, MPFR_RNDN);
        osc_refine_zero_mpfr(fmp, mpctx, comp, x1m, 0.5 * (z - x), bits);
        double perr;
        bool ccok = ccquad_mpfr_upto(fmp, mpctx, x0m, x1m, bits, panel_tol, NULL, cc, 4,
                              pr, pi, &perr);
        if (!ccok)
            denint_tanhsinh_mpfr(fmp, mpctx, x0m, x1m, bits, panel_tol,
                                 panel_levels, pr, pi, &perr);
        /* Accept the panel's best estimate whenever it is a finite number
         * (mirrors the machine osc_panel, which keeps a non-converged panel). */
        if (mpfr_number_p(pr) && mpfr_number_p(pi)) {
            mpfr_add(total_r, total_r, pr, MPFR_RNDN);
            mpfr_add(total_i, total_i, pi, MPFR_RNDN);
        }
        mpfr_set(x0m, x1m, MPFR_RNDN);
        scale = z - x;
        x = z;

        mpfr_init2(Sr[n], (mpfr_prec_t)bits); mpfr_init2(Si[n], (mpfr_prec_t)bits);
        mpfr_set(Sr[n], total_r, MPFR_RNDN);  mpfr_set(Si[n], total_i, MPFR_RNDN);
        n++;
        double mag = osc_mpfr_mag(total_r, total_i);
        if (mag > maxS) maxS = mag;

        if (n >= 5) {
            int deg = (n - 1) / 2; if (deg > dcap) deg = dcap; if (deg < 1) deg = 1;
            /* the extrapolation reads the most recent partial sums: a window
             * of 4 deg + 8 keeps each call O(deg^2) however many lobes */
            int w = n; if (w > 4 * deg + 8) w = 4 * deg + 8;
            double sstep; bool fin;
            if (seqaccel_wynn_mpfr((const mpfr_t*)Sr + (n - w), (const mpfr_t*)Si + (n - w), w, deg,
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

    /* The best extrapolation, converged or not: *abserr is its step, so the
     * value it describes is that one (the raw partial sum of a 1/n lobe
     * series is far worse: 4 digits where the extrapolation held 30).  The
     * raw sum only when no extrapolation was ever accepted. */
    if (isfinite(err)) { mpfr_set(out_re, bestr,   MPFR_RNDN); mpfr_set(out_im, besti,   MPFR_RNDN); }
    else               { mpfr_set(out_re, total_r, MPFR_RNDN); mpfr_set(out_im, total_i, MPFR_RNDN); }
    *abserr = err;
    ccquad_cache_free(cc);

    for (int i = 0; i < n; i++) { mpfr_clear(Sr[i]); mpfr_clear(Si[i]); }
    free(Sr); free(Si);
    mpfr_clears(total_r, total_i, pr, pi, x0m, x1m,
                accr, acci, bestr, besti, (mpfr_ptr)0);
    return conv;
}
#endif /* USE_MPFR */
