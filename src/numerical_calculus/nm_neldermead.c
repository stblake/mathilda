/* nm_neldermead.c — NMinimize NelderMead global engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* Build the (n+1)-vertex initial simplex from a user "InitialPoints" list. Each
 * element must be a length-n list of reals (evaluated + numericalized, so Pi
 * etc. are fine). If fewer than n+1 points are given, the remaining vertices are
 * seeded by perturbing the first point along successive axes; extra points are
 * ignored. Returns false (→ fall back to a random simplex) if the list is
 * malformed or any used point has the wrong dimension or a non-numeric entry. */
bool nm_simplex_from_points(NmDriver* D, const Expr* pts, size_t n, double* V) {
    if (!pts || !nm_is_head(pts, SYM_List) || pts->data.function.arg_count == 0)
        return false;
    size_t npts = pts->data.function.arg_count;
    size_t use  = npts < (n + 1) ? npts : (n + 1);
    for (size_t i = 0; i < use; i++) {
        Expr* p = pts->data.function.args[i];
        if (!nm_is_head(p, SYM_List) || p->data.function.arg_count != n) return false;
        for (size_t j = 0; j < n; j++) {
            Expr* c  = eval_and_free(expr_copy(p->data.function.args[j]));
            Expr* cn = c ? numericalize(c, numeric_machine_spec()) : NULL;
            double v;
            bool ok = cn && fm_expr_to_double_real(cn, &v);
            expr_free(c); expr_free(cn);
            if (!ok) return false;
            V[i * n + j] = v;
        }
    }
    for (size_t i = use; i <= n; i++) {
        for (size_t j = 0; j < n; j++) V[i * n + j] = V[j];
        size_t d = i - 1;
        double step = (D->reg_hi[d] - D->reg_lo[d]) * 0.1;
        if (step == 0.0) step = 1.0;
        V[i * n + d] += step;
    }
    for (size_t i = 0; i <= n; i++) nm_project(D, &V[i * n]);
    return true;
}

/* NelderMead downhill simplex on the penalized objective, with restarts. Each
 * restart's converged best vertex is polished into its basin minimum and the
 * restarts are ranked by those minima (default Min[2 n, 20] restarts), so a
 * multimodal surface improves with more restarts instead of being decided by the
 * raw simplex vertices. Gated by "PostProcess" -> False. */
void nm_neldermead(NmDriver* D, const NmConfig* nc, NmRng* rng,
                          double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* Explicit "SearchPoints" is honored verbatim (it was silently clamped at 20
     * before, contradicting the documented "honored verbatim"). The automatic
     * default runs Min[2 n, 20] random-simplex restarts rather than a flat 4, so
     * a multimodal surface is probed from several basins — the same "sample many
     * independent starts, keep the deepest" principle the other engines use. */
    int restarts;
    if (nc->search_points > 0) restarts = nc->search_points;
    else restarts = n > 1 ? (2 * (int)n < 20 ? 2 * (int)n : 20) : 2;
    if (restarts < 1) restarts = 1;
    /* Simplex coefficients: reflection (ReflectRatio, default 1), expansion
     * (ExpandRatio, default 2), contraction toward the centroid (ContractRatio,
     * default 0.5), shrink toward the best vertex (ShrinkRatio, default 0.5).
     * Tolerance is the objective-spread convergence threshold. */
    double rr  = nc->reflect_ratio  > 0.0 ? nc->reflect_ratio  : 1.0;
    double er  = nc->expand_ratio   > 0.0 ? nc->expand_ratio   : 2.0;
    double cr  = nc->contract_ratio > 0.0 ? nc->contract_ratio : 0.5;
    double sr  = nc->shrink_ratio   > 0.0 ? nc->shrink_ratio   : 0.5;
    double tol = nc->tolerance      > 0.0 ? nc->tolerance      : 1e-12;
    int64_t maxit = D->opts->max_iter > 0 ? D->opts->max_iter * (int64_t)(5 * n)
                                          : 200 * (int64_t)n;
    if (maxit < 100) maxit = 100;

    /* Domain-convergence scale: the simplex must shrink to a small fraction of
     * the search region — not just reach a flat objective — before we declare
     * convergence. Without this, a broad plateau (the flat outer region of the
     * Easom function, where f ≈ 0 everywhere away from a narrow spike) trips the
     * objective-spread test on the first iteration and the simplex never moves. */
    double rscale = 1.0;
    for (size_t j = 0; j < n; j++) { double e = rhi[j] - rlo[j]; if (e > rscale) rscale = e; }
    double xdtol = 1e-6 * rscale;

    double* V  = (double*)malloc(sizeof(double) * (n + 1) * n);
    double* fv = (double*)malloc(sizeof(double) * (n + 1));
    double* xc = (double*)malloc(sizeof(double) * n);
    double* xr = (double*)malloc(sizeof(double) * n);
    double* xe = (double*)malloc(sizeof(double) * n);
    bool have = false;

    for (int rs = 0; rs < restarts; rs++) {
        /* Restart 0 uses the user's "InitialPoints" simplex when supplied and
         * valid; every other restart (and the fallback) is a random simplex. */
        bool seeded = (rs == 0) && nm_simplex_from_points(D, nc->init_points, n, V);
        if (!seeded) {
            for (size_t j = 0; j < n; j++) V[j] = nm_rng_range(rng, rlo[j], rhi[j]);
            nm_project(D, &V[0]);
            for (size_t i = 1; i <= n; i++) {
                for (size_t j = 0; j < n; j++) V[i * n + j] = V[j];
                size_t d = i - 1;
                double step = (rhi[d] - rlo[d]) * 0.1;
                if (step == 0.0) step = 1.0;
                V[i * n + d] += step;
                nm_project(D, &V[i * n]);
            }
        }
        for (size_t i = 0; i <= n; i++) fv[i] = nm_phi(D, &V[i * n]);

        for (int64_t it = 0; it < maxit; it++) {
            size_t lo = 0, hi = 0, nh = 0;
            for (size_t i = 1; i <= n; i++) {
                if (fv[i] < fv[lo]) lo = i;
                if (fv[i] > fv[hi]) hi = i;
            }
            nh = (hi == 0) ? 1 : 0;
            for (size_t i = 0; i <= n; i++)
                if (i != hi && fv[i] > fv[nh]) nh = i;
            double xspread = 0.0;
            for (size_t i = 0; i <= n; i++)
                for (size_t j = 0; j < n; j++) {
                    double dd = fabs(V[i * n + j] - V[lo * n + j]);
                    if (dd > xspread) xspread = dd;
                }
            if (fabs(fv[hi] - fv[lo]) <= tol * (1.0 + fabs(fv[lo])) && xspread <= xdtol)
                break;

            for (size_t j = 0; j < n; j++) xc[j] = 0.0;
            for (size_t i = 0; i <= n; i++)
                if (i != hi)
                    for (size_t j = 0; j < n; j++) xc[j] += V[i * n + j];
            for (size_t j = 0; j < n; j++) xc[j] /= (double)n;

            for (size_t j = 0; j < n; j++) xr[j] = xc[j] + rr * (xc[j] - V[hi * n + j]);
            nm_project(D, xr);
            double frr = nm_phi(D, xr);
            if (frr < fv[lo]) {
                for (size_t j = 0; j < n; j++)
                    xe[j] = xc[j] + er * (xc[j] - V[hi * n + j]);
                nm_project(D, xe);
                double fe = nm_phi(D, xe);
                if (fe < frr) { for (size_t j = 0; j < n; j++) V[hi * n + j] = xe[j]; fv[hi] = fe; }
                else          { for (size_t j = 0; j < n; j++) V[hi * n + j] = xr[j]; fv[hi] = frr; }
            } else if (frr < fv[nh]) {
                for (size_t j = 0; j < n; j++) V[hi * n + j] = xr[j];
                fv[hi] = frr;
            } else {
                for (size_t j = 0; j < n; j++)
                    xe[j] = xc[j] + cr * (V[hi * n + j] - xc[j]);
                nm_project(D, xe);
                double fc = nm_phi(D, xe);
                if (fc < fv[hi]) { for (size_t j = 0; j < n; j++) V[hi * n + j] = xe[j]; fv[hi] = fc; }
                else {
                    for (size_t i = 0; i <= n; i++) {
                        if (i == lo) continue;
                        for (size_t j = 0; j < n; j++)
                            V[i * n + j] = V[lo * n + j] + sr * (V[i * n + j] - V[lo * n + j]);
                        nm_project(D, &V[i * n]);
                        fv[i] = nm_phi(D, &V[i * n]);
                    }
                }
            }
        }
        size_t lo = 0;
        for (size_t i = 1; i <= n; i++) if (fv[i] < fv[lo]) lo = i;
        double f, p;
        for (size_t j = 0; j < n; j++) xr[j] = V[lo * n + j];  /* xr: polish buffer */
        nm_eval(D, xr, &f, &p);
        /* Polish this restart's best vertex into its basin minimum before ranking
         * restarts, rather than ranking them by their raw simplex vertices and
         * polishing only the single winner afterward. The lowest converged vertex
         * across restarts need not sit in the deepest basin, so ranking by local
         * minima is what lets more restarts help rather than hurt. Skipped under
         * "PostProcess" -> False (then the raw vertex is ranked, as before); a
         * BFGS overshoot falls back to the raw vertex. */
        if (nc->post_process != 0) {
            double fr = f, pr = p;
            nm_local_polish(D, xr, &fr, &pr);
            if (nm_better(f, p, fr, pr)) {
                for (size_t j = 0; j < n; j++) xr[j] = V[lo * n + j];
            } else { f = fr; p = pr; }
        }
        if (!have || nm_better(f, p, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = xr[j];
            *fbest = f; *penbest = p; have = true;
        }
    }
    free(V); free(fv); free(xc); free(xr); free(xe);
}
