/* nm_sa.c — NMinimize SimulatedAnnealing global engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* SimulatedAnnealing "BoltzmannExponent" -> f. Evaluate f[i, df, f0] to the
 * real exponent whose Exp is the Metropolis acceptance probability for an
 * uphill move: i is the (1-based) iteration, df ≥ 0 the objective increase, f0
 * the current objective. Evaluator numeric diagnostics are muted, as everywhere
 * in the trial-point loop. Returns false — so the caller falls back to the
 * built-in geometric-cooling exponent — if f does not yield a finite real. */
static bool nm_boltzmann_exponent(Expr* bf, int64_t i, double df, double f0,
                                  double* out) {
    Expr* a[3];
    a[0] = expr_new_integer(i);
    a[1] = expr_new_real(df);
    a[2] = expr_new_real(f0);
    Expr* call = expr_new_function(expr_copy(bf), a, 3);
    arith_warnings_mute_push();
    Expr* v = eval_and_free(call);
    arith_warnings_mute_pop();
    bool ok = v && fm_expr_to_double_real(v, out) && isfinite(*out);
    expr_free(v);
    return ok;
}

/* SimulatedAnnealing with geometric cooling; tracks the best point seen.
 *
 * Honors three "SimulatedAnnealing" sub-options:
 *   "SearchPoints" -> K       run K independent annealing chains from random
 *                             starts and keep the best local minimum (default
 *                             Automatic = Min[2 n, 50], following Mathematica);
 *   "PerturbationScale" -> s  multiply the trial-step size by s (default 1.0);
 *   "BoltzmannExponent" -> f  use Exp[f[i, df, f0]] as the acceptance
 *                             probability for an uphill move (default -df/T).
 *
 * A rugged, many-basin landscape (Griewank, Rastrigin, ...) is not solved by a
 * single annealing walk: which basin one walk lands in is close to luck. So the
 * default runs Min[2 n, 50] independent chains and — crucially — polishes the
 * best raw point of *each* chain into its basin minimum before ranking them,
 * exactly as RandomSearch does per restart. Ranking chains by their local
 * minima rather than by their random-walk lows is what makes the reported
 * optimum improve, not degrade, as SearchPoints grows: the lowest point a walk
 * happens to visit is often in a shallower basin than a slightly higher point
 * that sits above a deeper one. The RNG is still drawn in the original order
 * within each chain, so a seeded single-chain run (SearchPoints -> 1) anneals
 * identically to before and reaches the same polished point the driver's
 * post-process already produced. */
void nm_sa(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    double pscale = nc->perturb_scale > 0.0 ? nc->perturb_scale : 1.0;
    Expr*  bf     = nc->boltzmann_fn;

    /* "SearchPoints" -> K independent annealing chains. Default Automatic =
     * Min[Max[2 n, 12], 50]. The 2 n scaling follows Mathematica, but a rugged,
     * many-basin landscape in low dimension (Eggholder, Schwefel, Griewank, ...)
     * has far more basins than 2 n = 4 starts can cover, so which basin the run
     * reports is left to luck. A floor of 12 independent starts is what turns
     * these from "sometimes finds the global" into "reliably finds it" at a cost
     * of a few hundredths of a second — SA is opt-in (the automatic method is
     * DifferentialEvolution), so only a caller who asked for it pays. */
    int64_t K;
    if (nc->search_points > 0) {
        K = (int64_t)nc->search_points;
    } else {
        K = 2 * (int64_t)n;
        if (K < 12) K = 12;
        if (K > 50) K = 50;
    }
    if (K < 1) K = 1;

    /* Per-chain iteration budget = "LevelIterations" trial moves at each of the
     * MaxIterations temperature levels, so per_chain = MaxIterations *
     * LevelIterations — Mathematica's semantics, where LevelIterations is the
     * dwell at each level. The old hard-coded 50 was exactly that implicit
     * default. An explicit "LevelIterations" is honored verbatim (like
     * "SearchPoints"): the caller asked for that budget, so the automatic
     * runtime caps below are skipped. The automatic budget keeps runtime bounded
     * — a single chain keeps the original schedule, and many search points share
     * a bounded aggregate, with a floor so each chain still anneals. */
    int64_t level_it  = nc->level_iterations > 0 ? (int64_t)nc->level_iterations : 50;
    int64_t base_iter = D->opts->max_iter > 0 ? D->opts->max_iter : 100;
    int64_t per_chain = base_iter * level_it;
    if (nc->level_iterations <= 0) {
        if (per_chain > 20000) per_chain = 20000;
        if (K > 1 && K * per_chain > NM_SA_TOTAL_CAP) {
            per_chain = NM_SA_TOTAL_CAP / K;
            if (per_chain < 300) per_chain = 300;
        }
    }

    double* x  = (double*)malloc(sizeof(double) * n);
    double* xn = (double*)malloc(sizeof(double) * n);
    double* xc = (double*)malloc(sizeof(double) * n);  /* best raw point in chain */
    bool have = false;

    for (int64_t chain = 0; chain < K; chain++) {
        for (size_t j = 0; j < n; j++) x[j] = nm_rng_range(rng, rlo[j], rhi[j]);
        nm_project(D, x);
        double fx, px;
        nm_eval(D, x, &fx, &px);
        double phi = fx + NM_PENALTY_MU * px;
        /* Track this chain's best raw point separately, so it can be polished on
         * its own and compared against the other chains' basin minima. */
        for (size_t j = 0; j < n; j++) xc[j] = x[j];
        double fc = fx, pc = px;

        /* Adaptive temperature scale. The acceptance exponent -d/T only anneals
         * when T sits on the scale of the objective differences d. With a fixed
         * T in [1e-4, 1] and an objective ranging over hundreds (Eggholder,
         * Schwefel, ...), every uphill move is rejected and the "anneal"
         * degenerates into greedy descent from the start point — so only many
         * random restarts, never the walk itself, ever crossed a ridge into a
         * deeper basin. Probe a handful of full-temperature trial steps to
         * measure the typical |Δφ| and use it as the scale, so a move that
         * worsens φ by one characteristic step is accepted with probability
         * e^{-1/T}. The probe draws from its own RNG so the annealing walk's
         * stream is untouched — a seeded chain takes exactly the trajectory it
         * did before the scale was introduced. xn is scratch, overwritten by the
         * main loop's first proposal. */
        double scale = 0.0;
        {
            NmRng prng;
            nm_rng_seed(&prng, nc->seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(chain + 1)));
            int nb = 0;
            for (int b = 0; b < NM_SA_BURN_IN; b++) {
                for (size_t j = 0; j < n; j++) {
                    double span = rhi[j] - rlo[j];
                    xn[j] = x[j] + pscale * (span * 0.1 * 1.1 * nm_rng_normal(&prng));
                }
                nm_project(D, xn);
                double fb, pb;
                nm_eval(D, xn, &fb, &pb);
                double db = fabs((fb + NM_PENALTY_MU * pb) - phi);
                if (isfinite(db)) { scale += db; nb++; }
            }
            scale = nb > 0 ? scale / nb : 0.0;
            if (!(scale > 0.0) || !isfinite(scale))
                scale = fabs(phi) > 1.0 ? fabs(phi) : 1.0;
        }

        double T = 1.0;
        for (int64_t it = 0; it < per_chain; it++) {
            for (size_t j = 0; j < n; j++) {
                double span = rhi[j] - rlo[j];
                xn[j] = x[j]
                      + pscale * (span * 0.1 * (0.1 + T) * nm_rng_normal(rng));
            }
            nm_project(D, xn);
            double fn2, pn2;
            nm_eval(D, xn, &fn2, &pn2);
            double phin = fn2 + NM_PENALTY_MU * pn2;
            double d = phin - phi;
            bool accept;
            if (d < 0.0) {
                accept = true;
            } else {
                double expo;
                if (!bf || !nm_boltzmann_exponent(bf, it + 1, d, phi, &expo))
                    expo = -d / (T * scale + 1e-12);
                accept = nm_rng_unif(rng) < exp(expo);   /* NaN prob ⇒ reject */
            }
            if (accept) {
                for (size_t j = 0; j < n; j++) x[j] = xn[j];
                phi = phin; fx = fn2; px = pn2;
                if (nm_better(fx, px, fc, pc)) {
                    for (size_t j = 0; j < n; j++) xc[j] = x[j];
                    fc = fx; pc = px;
                }
            }
            T *= 0.995;
            if (T < 1e-4) T = 1e-4;
        }

        /* Polish this chain's best raw point into its basin minimum, then keep
         * the global best of the polished candidates. A BFGS step can overshoot
         * a bound-projected point, so fall back to the pre-polish value if the
         * polish came out worse by Deb's rules (never worsens the chain).
         *
         * Skipped when the caller disabled polishing with "PostProcess" -> False:
         * then the per-chain raw best is carried straight through and the global
         * best over chains is the global raw best — bit-for-bit what a single
         * global-raw-best pass over the same RNG sequence produced before. */
        double fp = fc, pp = pc;
        for (size_t j = 0; j < n; j++) x[j] = xc[j];   /* x reused as polish buffer */
        if (nc->post_process != 0) {
            nm_local_polish(D, x, &fp, &pp);
            if (nm_better(fc, pc, fp, pp)) {
                for (size_t j = 0; j < n; j++) x[j] = xc[j];
                fp = fc; pp = pc;
            }
        }
        if (!have || nm_better(fp, pp, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = x[j];
            *fbest = fp; *penbest = pp; have = true;
        }
    }
    free(x); free(xn); free(xc);
}
