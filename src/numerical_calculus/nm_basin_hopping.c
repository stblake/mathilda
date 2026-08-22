/* nm_basin_hopping.c — NMinimize BasinHopping global engine.
 * Monte-Carlo minimization (Wales & Doye 1997), mirroring
 * scipy.optimize.basinhopping. Split from the original findmin.c; shared
 * declarations in findmin_internal.h. Do not add cross-file helpers here without
 * a prototype in that header. */
#include "findmin_internal.h"   /* pulls in <math.h> */

/* One basin-hopping run from a random start: quench the start, then hop
 * (random displacement + local minimization + Metropolis accept) for `niter`
 * steps with an adaptive step size, tracking the Deb-best in xbest/fbest/penbest
 * (which the caller has already seeded across earlier runs). Returns the run's
 * accepted current-point energy is not needed by the caller. */
static void bh_run(NmDriver* D, NmRng* rng, double T, double step0, int interval,
                   double target, double factor, int64_t niter,
                   int64_t niter_success, double* xcur, double* xtrial,
                   double* xbest, double* fbest, double* penbest, bool* have) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    double beta = (T > 0.0) ? 1.0 / T : HUGE_VAL;   /* T -> 0 => only downhill  */

    /* Initial point: uniform in the box, then local-minimize (scipy quenches the
     * user's x0 before the first hop; here the start is a fresh random point). */
    for (size_t j = 0; j < n; j++) xcur[j] = nm_rng_range(rng, rlo[j], rhi[j]);
    double fcur, pcur;
    nm_local_polish(D, xcur, &fcur, &pcur);
    double Ecur = fcur + NM_PENALTY_MU * pcur;
    if (!*have || nm_better(fcur, pcur, *fbest, *penbest)) {
        for (size_t j = 0; j < n; j++) xbest[j] = xcur[j];
        *fbest = fcur; *penbest = pcur; *have = true;
    }

    double step = step0;
    int64_t naccept = 0;         /* accepted hops so far (for the adaptation)     */
    int64_t nstep = 0;           /* hops taken so far                             */
    int64_t since_improve = 0;   /* hops since the global best last improved      */

    for (int64_t it = 0; it < niter; it++) {
        /* Adaptive step-size update, exactly scipy's AdaptiveStepsize: increment
         * the step counter first, then adjust on every `interval`-th hop using
         * the running acceptance rate. rate > target means we accept too many
         * (trapped in a basin) so take bigger steps; else take smaller ones. */
        nstep++;
        if (interval > 0 && (nstep % interval) == 0) {
            double rate = (double)naccept / (double)nstep;
            if (rate > target) step /= factor;
            else               step *= factor;
        }

        /* Random displacement of the current point, clipped to the box (a bounded
         * local minimizer would pull an out-of-box step back anyway). */
        for (size_t j = 0; j < n; j++) {
            double v = xcur[j] + nm_rng_range(rng, -step, step);
            if (v < rlo[j]) v = rlo[j];
            if (v > rhi[j]) v = rhi[j];
            xtrial[j] = v;
        }
        double ft, pt;
        nm_local_polish(D, xtrial, &ft, &pt);        /* quench the trial point   */
        double Et = ft + NM_PENALTY_MU * pt;

        /* Generalized Metropolis on the penalized, locally-minimized energies. */
        bool accept;
        double d = Et - Ecur;
        if (d <= 0.0) {
            accept = true;
        } else {
            double w = exp(-d * beta);               /* in (0, 1]; T -> 0 => 0   */
            accept = nm_rng_unif(rng) < w;
        }
        if (accept) {
            for (size_t j = 0; j < n; j++) xcur[j] = xtrial[j];
            Ecur = Et; fcur = ft; pcur = pt;
            naccept++;
        }

        /* Track the global best over every quenched point (accepted or not). */
        if (nm_better(ft, pt, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = xtrial[j];
            *fbest = ft; *penbest = pt; *have = true;
            since_improve = 0;
        } else {
            since_improve++;
        }
        if (niter_success > 0 && since_improve >= niter_success) break;
    }
}

void nm_basin_hopping(NmDriver* D, const NmConfig* nc, NmRng* rng,
                             double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    *fbest = 1e300; *penbest = 1e300;
    bool have = false;

    if (n == 0) {
        double d0 = 0.0, f, p;
        nm_eval(D, &d0, &f, &p);
        *fbest = f; *penbest = p;
        return;
    }

    /* Parameters (validated at parse time; clamp defensively in case an engine is
     * reached with an unvalidated config). */
    double T        = nc->bh_temp > 0.0 ? nc->bh_temp : NM_BH_TEMP;
    double step0    = nc->bh_step > 0.0 ? nc->bh_step : NM_BH_STEP;
    int    interval = nc->bh_interval > 0 ? nc->bh_interval : NM_BH_INTERVAL;
    double target   = (nc->bh_target_accept > 0.0 && nc->bh_target_accept < 1.0)
                    ? nc->bh_target_accept : NM_BH_TARGET_ACCEPT;
    double factor   = (nc->bh_step_factor > 0.0 && nc->bh_step_factor < 1.0)
                    ? nc->bh_step_factor : NM_BH_STEP_FACTOR;
    /* Hop count: an explicit MaxIterations wins; otherwise scipy's niter default
     * of 100, which is exactly NMinimize's generic MaxIterations default — so no
     * DA-style budget special-casing is needed here. */
    int64_t niter = D->opts->max_iter > 0 ? D->opts->max_iter : NM_BH_NITER;
    int64_t niter_success = nc->bh_niter_success > 0 ? nc->bh_niter_success : 0;
    int64_t K = nc->search_points > 0 ? (int64_t)nc->search_points : 1;

    double* xcur   = (double*)malloc(sizeof(double) * n);
    double* xtrial = (double*)malloc(sizeof(double) * n);

    for (int64_t run = 0; run < K; run++)
        bh_run(D, rng, T, step0, interval, target, factor, niter, niter_success,
               xcur, xtrial, xbest, fbest, penbest, &have);

    /* Snap integer coordinates and clamp to the box, then re-score so the
     * reported point and value are consistent (the driver re-polishes this best
     * unless "PostProcess" -> False). */
    if (have) {
        nm_project(D, xbest);
        nm_eval(D, xbest, fbest, penbest);
    }
    free(xcur); free(xtrial);
}
