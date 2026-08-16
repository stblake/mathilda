/* findmin_neldermead.c — Nelder-Mead simplex derivative-free local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


bool fm_run_neldermead(Expr* f, Expr** vars, size_t n,
                              FmVarBind* binds, Expr** g_exprs,
                              double* x, /* in/out */
                              const FmGenCon* gens, size_t ngens, double mu,
                              const FmBox* boxes,
                              const FmOpts* opts,
                              double* fx_out) {
    (void)vars; (void)g_exprs; (void)gens; (void)ngens; (void)mu;

    /* n == 1: a 2-vertex simplex is just a crude line search; delegate to the
     * exact Brent 1-D minimiser (identical optimum, cheaper), as Powell does. */
    if (n == 1) {
        if (boxes) fm_project_box(x, 1, boxes);
        double a, b, c;
        const FmBox* box1 = boxes ? &boxes[0] : NULL;
        if (box1 && box1->has_lo && box1->has_hi) {
            a = box1->lo; c = box1->hi; b = 0.5 * (a + c);
            if (x[0] > a && x[0] < c) b = x[0];
        } else if (!fm_bracket(f, binds, opts, x[0], box1, &a, &b, &c)) {
            fm_warn(g_fm_name, "nlnum", "bracket-finding failed");
            return false;
        }
        double xm, fmv;
        bool ok1 = fm_brent_min(f, binds, opts, a, b, c, box1, &xm, &fmv);
        if (ok1) { x[0] = xm; *fx_out = fmv; }
        return ok1;
    }

    const double rho = 1.0, chi = 2.0, psi = 0.5, sigma = 0.5;
    size_t np = n + 1;
    double* sim  = (double*)malloc(np * n * sizeof(double)); /* row k = vertex k */
    double* fsim = (double*)malloc(np * sizeof(double));
    size_t* idx  = (size_t*)malloc(np * sizeof(size_t));
    double* xbar = (double*)malloc(n * sizeof(double));      /* centroid       */
    double* xr   = (double*)malloc(n * sizeof(double));      /* reflection     */
    double* xt   = (double*)malloc(n * sizeof(double));      /* expand/contract*/
    bool ok = false;
    if (!sim || !fsim || !idx || !xbar || !xr || !xt) goto cleanup;

    if (boxes) fm_project_box(x, n, boxes);

    /* Initial simplex: vertex 0 = start; vertex k+1 = start with coordinate k
     * perturbed by 5% (or 0.00025 if that coordinate is 0), scipy's rule. */
    for (size_t j = 0; j < n; j++) sim[j] = x[j];
    for (size_t k = 0; k < n; k++) {
        double* v = &sim[(k + 1) * n];
        for (size_t j = 0; j < n; j++) v[j] = x[j];
        v[k] = (x[k] != 0.0) ? (1.05 * x[k]) : 0.00025;
        if (boxes) fm_project_box(v, n, boxes);
    }
    for (size_t k = 0; k < np; k++) {
        if (!fm_eval_scalar(f, binds, &sim[k * n], n, opts, &fsim[k])) {
            fm_warn(g_fm_name, "nlnum", "objective evaluation failed building the simplex");
            goto cleanup;
        }
        idx[k] = k;
    }
    fm_nm_sort_idx(idx, fsim, np);

    double fatol = pow(10.0, -opts->acc_goal_digits);
    double xatol = pow(10.0, -opts->prec_goal_digits);

    for (int64_t it = 0; it < opts->max_iter; it++) {
        size_t best = idx[0], worst = idx[n], second = idx[n - 1];

        /* Convergence: worst-to-best spread in both f and x below tolerance.
         * (idx is sorted, so the max f-gap is fsim[worst]-fsim[best].) */
        double fspread = fabs(fsim[worst] - fsim[best]);
        double xspread = 0.0;
        {
            const double* bpt = &sim[best * n];
            for (size_t k = 1; k < np; k++) {
                const double* v = &sim[idx[k] * n];
                for (size_t j = 0; j < n; j++) {
                    double dx = fabs(v[j] - bpt[j]);
                    if (dx > xspread) xspread = dx;
                }
            }
        }
        if (fspread <= fatol && xspread <= xatol) { ok = true; break; }

        /* Centroid of the n best vertices (all but the worst). */
        for (size_t j = 0; j < n; j++) xbar[j] = 0.0;
        for (size_t k = 0; k < n; k++) {
            const double* v = &sim[idx[k] * n];
            for (size_t j = 0; j < n; j++) xbar[j] += v[j];
        }
        for (size_t j = 0; j < n; j++) xbar[j] /= (double)n;

        const double* worstv = &sim[worst * n];
        /* Reflection xr = (1+rho)*xbar - rho*worst. */
        for (size_t j = 0; j < n; j++) xr[j] = (1.0 + rho) * xbar[j] - rho * worstv[j];
        if (boxes) fm_project_box(xr, n, boxes);
        double fxr;
        bool okr = fm_eval_scalar(f, binds, xr, n, opts, &fxr);
        if (!okr) fxr = HUGE_VAL;   /* a non-finite reflection is "very bad" → contracts */

        bool doshrink = false;
        if (okr && fxr < fsim[best]) {
            /* Expansion xe = (1+rho*chi)*xbar - rho*chi*worst; keep the better. */
            for (size_t j = 0; j < n; j++) xt[j] = (1.0 + rho * chi) * xbar[j] - rho * chi * worstv[j];
            if (boxes) fm_project_box(xt, n, boxes);
            double fxe;
            bool oke = fm_eval_scalar(f, binds, xt, n, opts, &fxe);
            if (oke && fxe < fxr) {
                for (size_t j = 0; j < n; j++) sim[worst * n + j] = xt[j];
                fsim[worst] = fxe;
            } else {
                for (size_t j = 0; j < n; j++) sim[worst * n + j] = xr[j];
                fsim[worst] = fxr;
            }
        } else if (okr && fxr < fsim[second]) {
            /* Reflection accepted. */
            for (size_t j = 0; j < n; j++) sim[worst * n + j] = xr[j];
            fsim[worst] = fxr;
        } else if (okr && fxr < fsim[worst]) {
            /* Outside contraction xc = (1+psi*rho)*xbar - psi*rho*worst. */
            for (size_t j = 0; j < n; j++) xt[j] = (1.0 + psi * rho) * xbar[j] - psi * rho * worstv[j];
            if (boxes) fm_project_box(xt, n, boxes);
            double fxc;
            bool okc = fm_eval_scalar(f, binds, xt, n, opts, &fxc);
            if (okc && fxc <= fxr) {
                for (size_t j = 0; j < n; j++) sim[worst * n + j] = xt[j];
                fsim[worst] = fxc;
            } else doshrink = true;
        } else {
            /* Inside contraction xcc = (1-psi)*xbar + psi*worst. */
            for (size_t j = 0; j < n; j++) xt[j] = (1.0 - psi) * xbar[j] + psi * worstv[j];
            if (boxes) fm_project_box(xt, n, boxes);
            double fxcc;
            bool okcc = fm_eval_scalar(f, binds, xt, n, opts, &fxcc);
            if (okcc && fxcc < fsim[worst]) {
                for (size_t j = 0; j < n; j++) sim[worst * n + j] = xt[j];
                fsim[worst] = fxcc;
            } else doshrink = true;
        }

        if (doshrink) {
            /* Shrink every vertex but the best toward the best. `best` = idx[0]
             * is never rewritten, so `bpt` aliases no modified row. */
            const double* bpt = &sim[best * n];
            for (size_t k = 1; k < np; k++) {
                size_t vi = idx[k];
                double* v = &sim[vi * n];
                for (size_t j = 0; j < n; j++) v[j] = bpt[j] + sigma * (v[j] - bpt[j]);
                if (boxes) fm_project_box(v, n, boxes);
                if (!fm_eval_scalar(f, binds, v, n, opts, &fsim[vi])) fsim[vi] = HUGE_VAL;
            }
        }

        fm_nm_sort_idx(idx, fsim, np);
        fm_fire_monitor(opts->step_monitor);
    }

    {
        size_t best = idx[0];
        for (size_t j = 0; j < n; j++) x[j] = sim[best * n + j];
        *fx_out = fsim[best];
    }
    ok = true;
cleanup:
    free(sim); free(fsim); free(idx); free(xbar); free(xr); free(xt);
    return ok;
}
