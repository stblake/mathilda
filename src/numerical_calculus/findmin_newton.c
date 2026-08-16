/* findmin_newton.c — Newton local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Newton (machine precision, modified Cholesky)                       *
 * ------------------------------------------------------------------ */

bool fm_run_newton(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs, Expr*** H_exprs,
                          double* x,
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars;
    double* g = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* H = (double*)malloc(sizeof(double) * n * n);
    double* Hcopy = (double*)malloc(sizeof(double) * n * n);
    double* neg_g = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    if (boxes) fm_project_box(x, n, boxes);

    double fx;
    if (augmented) {
        if (!fm_eval_augmented(f, binds, x, n, gens, ngens, mu, opts, &fx)) goto cleanup;
    } else {
        if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) goto cleanup;
    }

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        bool gok;
        if (augmented) {
            gok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                       binds, x, n, opts, g);
        } else {
            gok = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
            if (!gok) gok = fm_grad_finite_diff(f, binds, x, n, opts, g);
        }
        if (!gok) {
            fm_warn(g_fm_name, "nlnum", "gradient failed during Newton");
            goto cleanup;
        }
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { ok = true; break; }

        bool Hok = H_exprs && fm_eval_hessian(H_exprs, binds, x, n, opts, H);
        if (!Hok) {
            /* Fall back to BFGS-style steepest. */
            for (size_t i = 0; i < n; i++) d[i] = -g[i];
        } else {
            /* Try Cholesky with increasing τ. */
            double tau = 0.0;
            bool factored = false;
            for (int t = 0; t < 30 && !factored; t++) {
                for (size_t i = 0; i < n*n; i++) Hcopy[i] = H[i];
                factored = fm_chol_factor(Hcopy, n, tau);
                if (!factored) tau = (tau == 0.0) ? 1e-3 : tau * 2.0;
            }
            if (!factored) {
                fm_warn(g_fm_name, "dsing", "Hessian not positive definite");
                for (size_t i = 0; i < n; i++) d[i] = -g[i];
            } else {
                for (size_t i = 0; i < n; i++) neg_g[i] = -g[i];
                fm_chol_solve(Hcopy, n, neg_g, d);
            }
        }
        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) { for (size_t i = 0; i < n; i++) d[i] = -g[i];
                              g_dot_d = 0.0;
                              for (size_t i = 0; i < n; i++) g_dot_d += g[i]*d[i]; }
        double alpha, fx_new;
        bool ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                    augmented ? gens : NULL,
                                    augmented ? ngens : 0,
                                    augmented ? mu : 0.0,
                                    boxes, opts, &alpha, &fx_new, x_new);
        if (!ls_ok) {
            if (!augmented) fm_warn(g_fm_name, "lstol", "Newton line search failed");
            break;
        }
        fm_fire_monitor(opts->step_monitor);
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
            x[i] = x_new[i];
        }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }
    }
    *fx_out = fx;
    ok = true;
cleanup:
    free(g); free(d); free(x_new); free(H); free(Hcopy); free(neg_g);
    return ok;
}
