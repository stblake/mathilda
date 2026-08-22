/* nm_de.c — NMinimize DifferentialEvolution global engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Global engines                                                     *
 * ------------------------------------------------------------------ */

/* Seed the leading members of a DE population from a user "InitialPoints" list
 * and return how many were seeded (the rest are filled at random by the caller).
 * Each seed must be a length-n list of reals (evaluated + numericalized, so Pi
 * etc. are fine); it is projected into the search region and its integer
 * coordinates are rounded. Unlike nm_simplex_from_points this is per-member
 * rather than all-or-nothing: a malformed point ends the seeding and the rest of
 * the population is random, so one bad point is not fatal. Extra points beyond
 * the population size are ignored. Returns 0 (→ fully random) when the list is
 * absent or empty. */
static size_t nm_population_from_points(NmDriver* D, const Expr* pts,
                                        size_t n, size_t NP, double* pop) {
    if (!pts || !nm_is_head(pts, SYM_List) || pts->data.function.arg_count == 0)
        return 0;
    size_t npts = pts->data.function.arg_count;
    size_t use  = npts < NP ? npts : NP;
    size_t seeded = 0;
    for (size_t i = 0; i < use; i++) {
        Expr* p = pts->data.function.args[i];
        if (!nm_is_head(p, SYM_List) || p->data.function.arg_count != n) break;
        double* row = &pop[i * n];
        bool ok = true;
        for (size_t j = 0; j < n; j++) {
            Expr* c  = eval_and_free(expr_copy(p->data.function.args[j]));
            Expr* cn = c ? numericalize(c, numeric_machine_spec()) : NULL;
            double v = 0.0;
            ok = cn && fm_expr_to_double_real(cn, &v);
            expr_free(c); expr_free(cn);
            if (!ok) break;
            row[j] = v;
        }
        if (!ok) break;
        nm_project(D, row);
        for (size_t j = 0; j < n; j++)
            if (D->is_int[j]) row[j] = round(row[j]);
        seeded++;
    }
    return seeded;
}

/* DifferentialEvolution, DE/rand/1/bin, with Deb-rule selection. */
void nm_de(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* An explicit "SearchPoints" is honored verbatim (only floored at the DE
     * minimum of 4 members its DE/rand/1 mutation needs). The automatic
     * population is Storn & Price's 10n — the textbook default — clamped to
     * [15, 200]: the floor keeps a tiny problem's population workable, and the
     * ceiling 200 = 10·20 bounds the per-generation cost on high-dimensional
     * problems. This sizing is the same whether the method is an explicit
     * "DifferentialEvolution" or Method -> Automatic (there is no longer a lower
     * historical clamp for the explicit form). NelderMead restarts and
     * RandomSearch starts honor SearchPoints verbatim. */
    size_t NP;
    if (nc->search_points > 0) {
        NP = (size_t)nc->search_points;
        if (NP < 4) NP = 4;
    } else {
        NP = 10 * n;
        if (NP < 15)  NP = 15;
        if (NP > 200) NP = 200;
    }
    double F  = nc->F  > 0.0 ? nc->F  : 0.6;
    double CR = nc->CR >= 0.0 ? nc->CR : 0.9;
    int64_t maxgen = D->opts->max_iter > 0 ? D->opts->max_iter : 100;
    /* Method -> Automatic with no explicit MaxIterations scales the generation
     * budget with dimension (150n). A deceptive multimodal landscape such as
     * Michalewicz-10 (Sin[...]^20 ridges, near-flat elsewhere, so the post-polish
     * cannot rescue a weak basin) needs far more than the flat 100 generations to
     * find a good basin. This is free on easy problems: the convergence
     * early-break below stops as soon as the population collapses. An explicit
     * "DifferentialEvolution" (or a user MaxIterations) keeps the flat budget. */
    if (nc->method == NM_AUTO && !D->opts->max_iter_set)
        maxgen = 150 * (int64_t)n;

    double* pop   = (double*)malloc(sizeof(double) * NP * n);
    double* fpop  = (double*)malloc(sizeof(double) * NP);
    double* ppop  = (double*)malloc(sizeof(double) * NP);
    double* trial = (double*)malloc(sizeof(double) * n);

    /* "InitialPoints" seeds the leading population members; the remainder are
     * random. When no list is supplied seeded == 0 and every member is random —
     * the RNG stream is then identical to before, so seeded runs reproduce. */
    size_t seeded = nm_population_from_points(D, nc->init_points, n, NP, pop);
    for (size_t p = 0; p < NP; p++) {
        if (p >= seeded) {
            for (size_t j = 0; j < n; j++) {
                double v = nm_rng_range(rng, rlo[j], rhi[j]);
                if (D->is_int[j]) v = round(v);
                pop[p * n + j] = v;
            }
        }
        nm_eval(D, &pop[p * n], &fpop[p], &ppop[p]);
    }
    size_t bi = 0;
    for (size_t p = 1; p < NP; p++)
        if (nm_better(fpop[p], ppop[p], fpop[bi], ppop[bi])) bi = p;
    for (size_t j = 0; j < n; j++) xbest[j] = pop[bi * n + j];
    *fbest = fpop[bi];
    *penbest = ppop[bi];

    for (int64_t g = 0; g < maxgen; g++) {
        for (size_t p = 0; p < NP; p++) {
            size_t r1, r2, r3;
            if (NP < 4) break;
            do { r1 = nm_rng_next(rng) % NP; } while (r1 == p);
            do { r2 = nm_rng_next(rng) % NP; } while (r2 == p || r2 == r1);
            do { r3 = nm_rng_next(rng) % NP; } while (r3 == p || r3 == r1 || r3 == r2);
            size_t jr = nm_rng_next(rng) % n;
            for (size_t j = 0; j < n; j++) {
                if (nm_rng_unif(rng) < CR || j == jr) {
                    double base = pop[r1 * n + j];
                    double v = base + F * (pop[r2 * n + j] - pop[r3 * n + j]);
                    /* Bounce-back on a bound violation instead of clamping to
                     * the bound. Clamping strands the search: once several
                     * members share the exact boundary value for coordinate j,
                     * their pairwise mutation differentials for j collapse to
                     * zero, no mutant can lift the coordinate off the wall, and
                     * it freezes there for the rest of the run — the Schwefel
                     * failure, where half the coordinates pinned to ±500 while
                     * the rest found the true basin. Bounce-back (Price, Storn &
                     * Lampinen 2005) reflects the overshoot to a random point
                     * between the base vector and the violated bound: feasible,
                     * but almost never exactly on the boundary, so the
                     * population keeps its diversity. */
                    if (v < rlo[j])      v = base + nm_rng_unif(rng) * (rlo[j] - base);
                    else if (v > rhi[j]) v = base + nm_rng_unif(rng) * (rhi[j] - base);
                    if (D->is_int[j]) v = round(v);
                    trial[j] = v;
                } else {
                    trial[j] = pop[p * n + j];
                }
            }
            double ft, pt;
            nm_eval(D, trial, &ft, &pt);
            if (nm_better(ft, pt, fpop[p], ppop[p])) {
                for (size_t j = 0; j < n; j++) pop[p * n + j] = trial[j];
                fpop[p] = ft; ppop[p] = pt;
                if (nm_better(ft, pt, *fbest, *penbest)) {
                    for (size_t j = 0; j < n; j++) xbest[j] = trial[j];
                    *fbest = ft; *penbest = pt;
                }
            }
        }
        /* Converged if the best is feasible and the feasible sub-population's
         * objective spread has collapsed to the requested tolerance. */
        if (*penbest <= NM_FEAS_RANK) {   /* loose: convergence probe, not a claim */
            double fmin = 1e300, fmax = -1e300;
            size_t cnt = 0;
            for (size_t p = 0; p < NP; p++) {
                if (ppop[p] <= NM_FEAS_RANK) {
                    if (fpop[p] < fmin) fmin = fpop[p];
                    if (fpop[p] > fmax) fmax = fpop[p];
                    cnt++;
                }
            }
            /* Convergence tolerance on the feasible sub-population's objective
             * spread. An explicit "Tolerance" sub-option sets it directly;
             * otherwise it is derived from AccuracyGoal / PrecisionGoal. */
            double tol = nc->tolerance > 0.0
                       ? nc->tolerance
                       : pow(10.0, -D->opts->acc_goal_digits)
                         + (1.0 + fabs(*fbest)) * pow(10.0, -D->opts->prec_goal_digits);
            if (cnt >= NP / 2 && (fmax - fmin) <= tol) break;
        }
    }

    /* Multi-start local polish over the final population. Polishing only the
     * single global best strands DE in whichever basin that one member happened
     * to occupy, so a larger population or a shorter run — both of which leave
     * the population less converged — can report a WORSE optimum (Griewank-10:
     * SearchPoints -> 100 gave 0.197 where the default 40 gave 0.044). This is
     * the same non-monotonicity the SimulatedAnnealing per-chain polish removed.
     * Polish the best Min[2 n, 50] *distinct* members — skipping any that sits in
     * an already-polished basin — and keep the deepest local minimum, ranking
     * basins by their minima rather than by the raw fitness of the member that
     * found them. Gated by "PostProcess" -> False (then the raw global best is
     * carried through as before) and confined to continuous problems, where the
     * BFGS polish is cheap; a mixed-integer run keeps its single driver polish so
     * its integer-descent cost is unchanged. */
    if (nc->post_process != 0 && !D->any_int && NP >= 1) {
        int64_t polish_cap = 2 * (int64_t)n < 50 ? 2 * (int64_t)n : 50;
        if (polish_cap > (int64_t)NP) polish_cap = (int64_t)NP;
        if (polish_cap < 1) polish_cap = 1;

        /* Two members belong to the same basin for dedup purposes when they agree
         * to within 1e-3 of the search span on every coordinate (floored, so a
         * degenerate zero-width span still separates distinct points). */
        double* btol = (double*)malloc(sizeof(double) * n);
        for (size_t j = 0; j < n; j++) {
            double t = 1e-3 * (rhi[j] - rlo[j]);
            btol[j] = t > 1e-6 ? t : 1e-6;
        }
        unsigned char* used = (unsigned char*)calloc(NP, 1);
        double* seeds = (double*)malloc(sizeof(double) * (size_t)polish_cap * n);
        double* xp = (double*)malloc(sizeof(double) * n);
        size_t nseeds = 0;
        int64_t done = 0;

        while (used && seeds && xp && btol && done < polish_cap) {
            /* best not-yet-considered member by Deb's rules */
            size_t best = NP;
            for (size_t p = 0; p < NP; p++) {
                if (used[p]) continue;
                if (best == NP || nm_better(fpop[p], ppop[p], fpop[best], ppop[best]))
                    best = p;
            }
            if (best == NP) break;
            used[best] = 1;

            /* Skip a member already represented by a polished basin seed. */
            bool dup = false;
            for (size_t s = 0; s < nseeds && !dup; s++) {
                bool same = true;
                for (size_t j = 0; j < n; j++)
                    if (fabs(pop[best * n + j] - seeds[s * n + j]) > btol[j]) {
                        same = false; break;
                    }
                if (same) dup = true;
            }
            if (dup) continue;   /* does not count against polish_cap */
            for (size_t j = 0; j < n; j++) seeds[nseeds * n + j] = pop[best * n + j];
            nseeds++;

            for (size_t j = 0; j < n; j++) xp[j] = pop[best * n + j];
            double fp = fpop[best], pp = ppop[best];
            nm_local_polish(D, xp, &fp, &pp);
            /* TIGHT from here down: these two comparisons decide what leaves DE
             * and becomes the answer. Under the loose threshold the raw
             * population member (violating by ~1e-4, so "feasible", objective
             * 1.9998) beat its own polished counterpart (feasible, 2.0) — the
             * shipped bug, at this exact line. */
            if (nm_better_return(fpop[best], ppop[best], fp, pp)) {  /* overshoot guard */
                for (size_t j = 0; j < n; j++) xp[j] = pop[best * n + j];
                fp = fpop[best]; pp = ppop[best];
            }
            if (nm_better_return(fp, pp, *fbest, *penbest)) {
                for (size_t j = 0; j < n; j++) xbest[j] = xp[j];
                *fbest = fp; *penbest = pp;
            }
            done++;
        }
        free(btol); free(used); free(seeds); free(xp);
    }

    free(pop); free(fpop); free(ppop); free(trial);
}
