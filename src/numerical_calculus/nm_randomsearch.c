/* nm_randomsearch.c — NMinimize RandomSearch global engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* RandomSearch: multiple random starts, each refined by the local solver, best
 * local minimum kept. This is already the "polish each start, rank by basin
 * depth" pattern; the two other multi-start engines were brought in line with it.
 * Pure multi-start local search has no global move, so on a search box far wider
 * than the optimum's basin (e.g. Griewank over [-600, 600], where the central
 * bowl is ~1e-11 of the volume in 10-D) no attainable number of random starts
 * reaches it — DifferentialEvolution / SimulatedAnnealing are the engines for
 * that shape. */
void nm_randomsearch(NmDriver* D, const NmConfig* nc, NmRng* rng,
                            double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* Explicit "SearchPoints" is honored verbatim — it was silently capped at 40,
     * so SearchPoints -> 1000 was a no-op that returned the 40-start result. The
     * automatic default keeps a bound (runtime is one local solve per start). */
    int K;
    if (nc->search_points > 0) {
        K = nc->search_points;          /* honored verbatim (was capped at 40) */
        if (K < 1) K = 1;
    } else {
        K = n > 1 ? (int)(8 * n) : 12;  /* automatic default, bounded for runtime */
        if (K < 4)  K = 4;
        if (K > 40) K = 40;
    }
    double* x = (double*)malloc(sizeof(double) * n);
    bool have = false;
    for (int k = 0; k < K; k++) {
        for (size_t j = 0; j < n; j++) x[j] = nm_rng_range(rng, rlo[j], rhi[j]);
        nm_project(D, x);
        double f, p;
        nm_local_polish(D, x, &f, &p);
        if (!have || nm_better(f, p, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = x[j];
            *fbest = f; *penbest = p; have = true;
        }
    }
    free(x);
}
