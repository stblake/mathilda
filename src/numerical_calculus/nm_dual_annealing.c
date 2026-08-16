/* nm_dual_annealing.c — NMinimize DualAnnealing (Generalized SA) global engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"   /* pulls in <math.h> */

/* M_PI is POSIX, not C99: glibc hides it under -std=c99 (macOS exposes it
 * anyway). Provide the standard fallback so the Tsallis visiting factors below
 * compile on Linux. See CLAUDE.md §10 / src/trig.c. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static void da_factors_init(DaFactors* F, double qv) {
    double factor2 = exp((4.0 - qv) * log(qv - 1.0));
    double factor3 = exp((2.0 - qv) * log(2.0) / (qv - 1.0));
    double factor5 = 1.0 / (qv - 1.0) - 0.5;
    double d1 = 2.0 - factor5;
    F->qv = qv;
    F->factor4p = sqrt(M_PI) * factor2 / (factor3 * (3.0 - qv));
    F->factor6  = M_PI * (1.0 - factor5) / sin(M_PI * (1.0 - factor5))
                / exp(lgamma(d1));
}

/* One raw visiting jump at temperature T (Tsallis "Visita" formula). */
static double da_visit_fn(NmRng* rng, const DaFactors* F, double T) {
    double qv = F->qv;
    double x = nm_rng_normal(rng);
    double y = nm_rng_normal(rng);
    if (fabs(y) < 1e-300) y = (y < 0.0) ? -1e-300 : 1e-300;
    double factor1 = exp(log(T) / (qv - 1.0));
    double factor4 = F->factor4p * factor1;
    x *= exp(-(qv - 1.0) * log(F->factor6 / factor4) / (3.0 - qv));
    double den = exp((qv - 1.0) * log(fabs(y)) / (3.0 - qv));
    return x / den;
}

/* Fold a jumped coordinate back into [lo, hi] periodically (scipy's fmod-based
 * wrapping), nudged off the exact lower bound. */
static double da_wrap(double v, double lo, double hi) {
    double range = hi - lo;
    if (!(range > 0.0)) return lo;         /* degenerate fixed coordinate */
    double a = v - lo;
    double b = fmod(a, range) + range;
    double w = fmod(b, range) + lo;
    if (fabs(w - lo) < NM_DA_MIN_VISIT) w += NM_DA_MIN_VISIT;
    if (w > hi) w = hi;
    if (w < lo) w = lo;
    return w;
}

/* Build x_visit from x: jump all coordinates (step < n) or only coordinate
 * step-n (step >= n), with tail clamping and box wrapping. */
static void da_visiting(NmRng* rng, const DaFactors* F, double T,
                        const double* lo, const double* hi, size_t n,
                        const double* x, size_t step, double* xv) {
    if (step < n) {
        for (size_t j = 0; j < n; j++) {
            double vis = da_visit_fn(rng, F, T);
            if (!isfinite(vis) || vis > NM_DA_TAIL_LIMIT)
                vis = NM_DA_TAIL_LIMIT * nm_rng_unif(rng);
            else if (vis < -NM_DA_TAIL_LIMIT)
                vis = -NM_DA_TAIL_LIMIT * nm_rng_unif(rng);
            xv[j] = da_wrap(x[j] + vis, lo[j], hi[j]);
        }
    } else {
        size_t idx = step - n;
        for (size_t j = 0; j < n; j++) xv[j] = x[j];
        double vis = da_visit_fn(rng, F, T);
        if (!isfinite(vis) || vis > NM_DA_TAIL_LIMIT)
            vis = NM_DA_TAIL_LIMIT * nm_rng_unif(rng);
        else if (vis < -NM_DA_TAIL_LIMIT)
            vis = -NM_DA_TAIL_LIMIT * nm_rng_unif(rng);
        xv[idx] = da_wrap(x[idx] + vis, lo[idx], hi[idx]);
    }
}

void nm_dual_annealing(NmDriver* D, const NmConfig* nc, NmRng* rng,
                              double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    *fbest = 1e300; *penbest = 1e300;
    bool have = false;

    if (n == 0) {
        double d0 = 0.0, f, p;
        nm_eval(D, &d0, &f, &p);
        *fbest = f; *penbest = p;
        return;
    }

    /* Parameters (validated at parse time; clamp qv defensively for the log/pow
     * domain in case an engine is reached with an unvalidated config). */
    double qv = nc->da_visit;
    if (!(qv > 1.0 && qv <= 3.0)) qv = NM_DA_VISIT;
    double qa = nc->da_accept;
    double T0 = nc->da_init_temp > 0.0 ? nc->da_init_temp : NM_DA_INIT_TEMP;
    double restart_ratio = (nc->da_restart_ratio > 0.0 && nc->da_restart_ratio < 1.0)
                         ? nc->da_restart_ratio : NM_DA_RESTART_RATIO;
    /* Inner local search runs unless disabled OR PostProcess -> False (which
     * asks for the raw, unpolished global-search result everywhere). */
    bool local_on = (nc->da_local_search != 0) && (nc->post_process != 0);

    /* Temperature-step budget. An explicit MaxIterations wins; otherwise use
     * Dual Annealing's own default (scipy's maxiter=1000), NOT NMinimize's
     * generic MaxIterations default of 100 — a 10x-too-small budget starves the
     * anneal on the larger boxes and is the difference between missing and
     * reaching the global (matches scipy). */
    int64_t maxiter = D->opts->max_iter_set ? D->opts->max_iter : NM_DA_MAXITER;
    int64_t K = nc->search_points > 0 ? (int64_t)nc->search_points : 1;

    DaFactors F;
    da_factors_init(&F, qv);
    double t1 = exp((qv - 1.0) * log(2.0)) - 1.0;   /* schedule numerator     */
    double T_restart = T0 * restart_ratio;

    double* xcur = (double*)malloc(sizeof(double) * n);
    double* xv   = (double*)malloc(sizeof(double) * n);
    double* xp   = (double*)malloc(sizeof(double) * n);  /* local-search scratch */

    for (int64_t chain = 0; chain < K; chain++) {
        for (size_t j = 0; j < n; j++) xcur[j] = nm_rng_range(rng, rlo[j], rhi[j]);
        double fcur, pcur;
        nm_eval(D, xcur, &fcur, &pcur);
        double Ecur = fcur + NM_PENALTY_MU * pcur;
        if (!have || nm_better(fcur, pcur, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = xcur[j];
            *fbest = fcur; *penbest = pcur; have = true;
        }

        int64_t si = 0, iteration = 0, nfev = 0, stall = 0;
        while (iteration < maxiter && nfev < NM_DA_MAXFUN) {
            double s = (double)si + 2.0;
            double t2 = exp((qv - 1.0) * log(s)) - 1.0;
            double T = T0 * t1 / t2;

            if (T < T_restart) {
                /* Reanneal: fresh current point; the temperature climbs back to
                 * T0 at si = 0. The best is preserved. Does not consume the
                 * iteration budget (matches scipy's reset-then-restart loop). */
                for (size_t j = 0; j < n; j++)
                    xcur[j] = nm_rng_range(rng, rlo[j], rhi[j]);
                nm_eval(D, xcur, &fcur, &pcur); nfev++;
                Ecur = fcur + NM_PENALTY_MU * pcur;
                if (nm_better(fcur, pcur, *fbest, *penbest)) {
                    for (size_t j = 0; j < n; j++) xbest[j] = xcur[j];
                    *fbest = fcur; *penbest = pcur;
                }
                si = 0;
                continue;
            }

            double Tstep = T / (double)(si + 1);
            /* scipy forces a local search after the hottest chain of each
             * temperature ramp (step 0) regardless of improvement, seeding an
             * early basin. */
            bool improved = (si == 0);

            for (size_t step = 0; step < 2 * n; step++) {
                da_visiting(rng, &F, T, rlo, rhi, n, xcur, step, xv);
                double fvv, pvv;
                nm_eval(D, xv, &fvv, &pvv); nfev++;
                double Ev = fvv + NM_PENALTY_MU * pvv;
                double d = Ev - Ecur;
                bool accept;
                if (d <= 0.0) {
                    accept = true;
                } else {
                    double base = 1.0 - (1.0 - qa) * d / Tstep;
                    double pqv = (base <= 0.0) ? 0.0 : exp(log(base) / (1.0 - qa));
                    accept = nm_rng_unif(rng) < pqv;
                }
                if (accept) {
                    for (size_t j = 0; j < n; j++) xcur[j] = xv[j];
                    Ecur = Ev; fcur = fvv; pcur = pvv;
                }
                if (nm_better(fvv, pvv, *fbest, *penbest)) {
                    for (size_t j = 0; j < n; j++) xbest[j] = xv[j];
                    *fbest = fvv; *penbest = pvv;
                    have = true; improved = true;
                }
                if (nfev >= NM_DA_MAXFUN) break;
            }

            if (local_on) {
                if (improved) {
                    /* Global best moved this chain (or the forced step-0 search):
                     * polish from the best, keep the refined basin as the new
                     * current (scipy's improved branch). */
                    for (size_t j = 0; j < n; j++) xp[j] = xbest[j];
                    double fp = *fbest, pp = *penbest;
                    nm_local_polish(D, xp, &fp, &pp); nfev++;
                    if (nm_better(fp, pp, *fbest, *penbest)) {
                        for (size_t j = 0; j < n; j++) { xbest[j] = xp[j]; xcur[j] = xp[j]; }
                        *fbest = fp; *penbest = pp;
                        fcur = fp; pcur = pp; Ecur = fp + NM_PENALTY_MU * pp;
                    }
                    stall = 0;
                } else if (++stall >= NM_DA_STALL_MAX) {
                    /* Long stagnation: force a polish from the current point and
                     * continue from it (scipy's not_improved_max_idx branch). */
                    for (size_t j = 0; j < n; j++) xp[j] = xcur[j];
                    double fp = fcur, pp = pcur;
                    nm_local_polish(D, xp, &fp, &pp); nfev++;
                    if (nm_better(fp, pp, *fbest, *penbest)) {
                        for (size_t j = 0; j < n; j++) xbest[j] = xp[j];
                        *fbest = fp; *penbest = pp; have = true;
                    }
                    for (size_t j = 0; j < n; j++) xcur[j] = xp[j];
                    fcur = fp; pcur = pp; Ecur = fp + NM_PENALTY_MU * pp;
                    stall = 0;
                }
            }
            si++; iteration++;
        }
    }

    /* Snap integer coordinates and clamp to the box, then re-score so the
     * reported point and value are consistent (nm_eval rounds integers during
     * scoring, but the stored coordinates are otherwise continuous). */
    if (have) {
        nm_project(D, xbest);
        nm_eval(D, xbest, fbest, penbest);
    }
    free(xcur); free(xv); free(xp);
}
