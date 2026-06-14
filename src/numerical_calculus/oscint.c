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
            double x = a + win * (double)i / (N - 1);
            double _Complex v;
            if (!f(ctx, x, &v)) continue;
            double re = creal(v), im = cimag(v);
            if (re < minRe) minRe = re; if (re > maxRe) maxRe = re;
            if (im < minIm) minIm = im; if (im > maxIm) maxIm = im;
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

/* Locate the next sign change (zero) of component `comp` of f strictly after x,
 * marching forward in steps of `step` (< the half-period) and bisecting. */
static bool osc_next_zero(GkSampleMachine f, void* ctx, int comp,
                          double x, double step, int max_steps, double* z) {
    double _Complex v;
    /* Start a quarter-step past x: the previous zero sits at x, where f ≈ 0, so
     * sampling there would spuriously re-detect it as the next crossing. */
    double cur = x + 0.25 * step;
    if (!f(ctx, cur, &v)) { cur += 0.01 * step; if (!f(ctx, cur, &v)) return false; }
    double f0 = comp ? cimag(v) : creal(v);
    for (int i = 0; i < max_steps; i++) {
        double nx = cur + step;
        double _Complex vn;
        if (!f(ctx, nx, &vn)) return false;
        double f1 = comp ? cimag(vn) : creal(vn);
        if ((f0 <= 0 && f1 >= 0) || (f0 >= 0 && f1 <= 0)) {
            double lo = cur, hi = nx, flo = f0;
            for (int bj = 0; bj < 50; bj++) {
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
        cur = nx; f0 = f1;
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
         * cleanly — exactly the sequence Wynn's epsilon accelerates. */
        double step = 0.4 * hp;
        double _Complex* S = malloc((size_t)max_panels * sizeof(double _Complex));
        if (!S) return false;
        int n = 0;
        double _Complex total = 0.0, best = 0.0;
        double err = INFINITY, maxS = 0.0;
        bool conv = false, have_best = false;
        double x = a;
        for (int k = 0; k < max_panels; k++) {
            double z;
            if (!osc_next_zero(f, ctx, comp, x, step, 100000, &z)) break;
            double _Complex p;
            if (osc_panel(f, ctx, x, z, panel_tol, &p)) total += p;
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
                    best = acc; err = sstep; have_best = true;
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

    /* Finite range: cover [a,b] with half-period panels and sum. */
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
    *result = total;
    *abserr = 0.0;
    return all_ok && x >= b;
}
