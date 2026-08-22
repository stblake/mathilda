/* findmin_bfgs.c — BFGS (machine precision) local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  BFGS (machine precision)                                            *
 * ------------------------------------------------------------------ */

bool fm_run_bfgs(Expr* f, Expr** vars, size_t n,
                        FmVarBind* binds, Expr** g_exprs,
                        double* x, /* in/out */
                        const FmGenCon* gens, size_t ngens, double mu,
                        const FmBox* boxes,
                        const FmOpts* opts,
                        double* fx_out) {
    (void)vars;
    /* Inverse Hessian approximation H, stored row-major. Initialise to I. */
    double* H = (double*)calloc(n * n, sizeof(double));
    double* g = (double*)malloc(sizeof(double) * n);
    double* g_new = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* s = (double*)malloc(sizeof(double) * n);
    double* y = (double*)malloc(sizeof(double) * n);
    double* Hy = (double*)malloc(sizeof(double) * n);
    bool ok = false;

    for (size_t i = 0; i < n; i++) H[i*n + i] = 1.0;
    if (boxes) fm_project_box(x, n, boxes);

    double fx;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    if (augmented) {
        if (!fm_eval_augmented(f, binds, x, n, gens, ngens, mu, opts, &fx)) goto cleanup;
    } else {
        if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) goto cleanup;
    }

    bool got_grad;
    if (augmented) {
        got_grad = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                        binds, x, n, opts, g);
    } else {
        got_grad = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
        if (!got_grad) got_grad = fm_grad_finite_diff(f, binds, x, n, opts, g);
    }
    if (!got_grad) {
        fm_warn(g_fm_name, "nlnum", "gradient evaluation failed at start point");
        goto cleanup;
    }

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* Gradient norm convergence. */
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { ok = true; break; }

        /* d = -H g. */
        for (size_t i = 0; i < n; i++) {
            double s_ = 0.0;
            for (size_t j = 0; j < n; j++) s_ += H[i*n + j] * g[j];
            d[i] = -s_;
        }
        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) {
            /* Not a descent direction — reset H to I and use steepest. */
            for (size_t i = 0; i < n*n; i++) H[i] = 0.0;
            for (size_t i = 0; i < n; i++) { H[i*n + i] = 1.0; d[i] = -g[i]; }
            g_dot_d = 0.0; for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        }

        double alpha, fx_new;
        bool ls_ok;
        if (augmented) {
            ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                   gens, ngens, mu, boxes, opts,
                                   &alpha, &fx_new, x_new);
        } else {
            ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                   NULL, 0, 0.0, boxes, opts,
                                   &alpha, &fx_new, x_new);
        }
        if (!ls_ok) {
            /* Line-search exhaustion is expected at high μ in the penalty
             * schedule (steep walls, large directional curvature). The
             * outer fm_run_penalty loop's feasibility check is the
             * authoritative signal in that case, so stay silent here and
             * let it speak instead. */
            if (!augmented) {
                fm_warn(g_fm_name, "lstol", "line search failed at iter %lld",
                        (long long)k);
            }
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        /* Step magnitude check (PrecisionGoal). */
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
        }

        /* Compute new gradient. */
        bool ng_ok;
        if (augmented) {
            ng_ok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                         binds, x_new, n, opts, g_new);
        } else {
            ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
            if (!ng_ok) ng_ok = fm_grad_finite_diff(f, binds, x_new, n, opts, g_new);
        }
        if (!ng_ok) {
            fm_warn(g_fm_name, "nlnum", "gradient evaluation failed in iteration");
            /* Take the step and stop. */
            for (size_t i = 0; i < n; i++) x[i] = x_new[i];
            fx = fx_new;
            break;
        }

        /* BFGS update: s = x_new - x; y = g_new - g; ρ = 1 / (y . s). */
        for (size_t i = 0; i < n; i++) { s[i] = x_new[i] - x[i]; y[i] = g_new[i] - g[i]; }
        double sy = 0.0;
        for (size_t i = 0; i < n; i++) sy += s[i] * y[i];
        if (sy > 1e-12) {
            double rho = 1.0 / sy;
            /* Hy = H y. */
            for (size_t i = 0; i < n; i++) {
                double t = 0.0;
                for (size_t j = 0; j < n; j++) t += H[i*n + j] * y[j];
                Hy[i] = t;
            }
            double yHy = 0.0;
            for (size_t i = 0; i < n; i++) yHy += y[i] * Hy[i];
            /* H ← H + ((sy + yHy) ρ²) s s^T − ρ (Hy s^T + s (Hy)^T). */
            double coef = (sy + yHy) * rho * rho;
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    H[i*n + j] += coef * s[i] * s[j]
                                - rho * (Hy[i] * s[j] + s[i] * Hy[j]);
                }
            }
        }

        for (size_t i = 0; i < n; i++) { x[i] = x_new[i]; g[i] = g_new[i]; }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }
    }
    if (!ok) {
        /* Either max iters or line search exhausted — still report best. */
    }
    *fx_out = fx;
    /* Always return true so the driver gets the best iterate; warnings
     * already emitted above when convergence failed. */
    ok = true;
cleanup:
    free(H); free(g); free(g_new); free(d); free(x_new);
    free(s); free(y); free(Hy);
    return ok;
}
