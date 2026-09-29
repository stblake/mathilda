/* findmin_cg.c — conjugate-gradient local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Conjugate gradient (Polak-Ribière+ with restart)                    *
 * ------------------------------------------------------------------ *
 * Nonlinear CG only works as well as its line search: the PR+ direction is
 * built on the assumption that the previous step ended (nearly) at a line
 * minimum, i.e. that g_new . d_old ~ 0. The shared Armijo backtracking search
 * (fm_line_search) never lengthens a step past ||alpha d|| = 1 and stops at the
 * first sufficient decrease, so on a curved valley such as Rosenbrock's CG
 * degenerated into short, poorly-conjugate steps and ran out its 500
 * iterations at f ~ 6e-4. The unconstrained, unbounded path therefore uses a
 * STRONG-WOLFE line search (Nocedal & Wright, Alg. 3.5/3.6, with c2 = 0.1 as
 * recommended for CG), seeded by the N&W (3.60) initial step
 * alpha_prev * (g_prev.d_prev)/(g.d), and restarts on Powell's orthogonality
 * test |g_new.g| >= 0.2 |g_new|^2. The penalty-wrapped (augmented) and
 * box-bounded solves keep the original projected Armijo search and n-periodic
 * restart unchanged, since the Wolfe curvature test is not meaningful on a
 * projected path. */

/* f and its gradient at x: exact symbolic (compiled) gradient first, central
 * differences as the fallback -- the same pattern as the main loop. */
static bool fm_cg_eval_fg(Expr* f, FmVarBind* binds, Expr** g_exprs,
                          const double* x, size_t n, const FmOpts* opts,
                          double* fv, double* g) {
    if (!fm_eval_scalar(f, binds, x, n, opts, fv) || !isfinite(*fv)) return false;
    bool ok = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
    if (!ok) ok = fm_grad_finite_diff(f, binds, x, n, opts, g);
    return ok;
}

/* Safeguarded minimiser of the cubic interpolating (a, fa, da) and
 * (b, fb, db) on the interval between a and b (N&W eq. 3.59); falls back to
 * bisection when the cubic is degenerate or its minimiser lies too close to
 * an endpoint. */
static double fm_cg_interp(double a, double fa, double da,
                           double b, double fb, double db) {
    double lo = (a < b) ? a : b, hi = (a < b) ? b : a;
    double d1 = da + db - 3.0 * (fa - fb) / (a - b);
    double disc = d1 * d1 - da * db;
    double t = 0.5 * (a + b);
    if (disc >= 0.0 && isfinite(disc)) {
        double d2 = sqrt(disc);
        if (b < a) d2 = -d2;
        double den = db - da + 2.0 * d2;
        if (den != 0.0) {
            double c = b - (b - a) * (db + d2 - d1) / den;
            if (isfinite(c)) t = c;
        }
    }
    double margin = 0.1 * (hi - lo);
    if (t < lo + margin || t > hi - margin) t = 0.5 * (a + b);
    return t;
}

/* Strong-Wolfe line search along d from x (f0, dphi0 = g.d < 0). On success
 * writes the accepted point to xt, its value to *ft and gradient to gt, the
 * step to *alpha_out, and returns true. When no point meets the strong Wolfe
 * conditions within the budget but some trial decreased f (Armijo), that best
 * trial is returned instead -- still a descent step. Returns false only when
 * no trial decreased f at all. */
static bool fm_cg_wolfe(Expr* f, FmVarBind* binds, Expr** g_exprs, size_t n,
                        const double* x, const double* d, double f0, double dphi0,
                        double alpha_init, const FmOpts* opts,
                        double* xt, double* gt, double* ft, double* alpha_out,
                        double* xbest, double* gbest) {
    const double c1 = 1e-4, c2 = 0.1;
    double a_prev = 0.0, f_prev = f0, dp_prev = dphi0;
    double alpha = alpha_init;
    double best_a = 0.0, best_f = f0;

#define FM_CG_TRY(a_, fa_, da_, okv_)                                          \
    do {                                                                       \
        for (size_t i_ = 0; i_ < n; i_++) xt[i_] = x[i_] + (a_) * d[i_];      \
        okv_ = fm_cg_eval_fg(f, binds, g_exprs, xt, n, opts, &(fa_), gt);      \
        if (okv_) {                                                            \
            (da_) = 0.0;                                                       \
            for (size_t i_ = 0; i_ < n; i_++) (da_) += gt[i_] * d[i_];         \
            if ((fa_) < best_f && (fa_) <= f0 + c1 * (a_) * dphi0) {           \
                best_f = (fa_); best_a = (a_);                                 \
                for (size_t i_ = 0; i_ < n; i_++) { xbest[i_] = xt[i_]; gbest[i_] = gt[i_]; } \
            }                                                                  \
        }                                                                      \
    } while (0)

    double lo = 0.0, flo = f0, dlo = dphi0, hi = 0.0, fhi = 0.0, dhi = 0.0;
    bool zoom = false;
    for (int it = 0; it < 25 && !zoom; it++) {
        double fa = 0.0, da = 0.0; bool okv;
        FM_CG_TRY(alpha, fa, da, okv);
        if (!okv) {                       /* outside the domain: shorten */
            alpha = a_prev + 0.5 * (alpha - a_prev);
            if (alpha - a_prev < 1e-16 * (1.0 + a_prev)) break;
            continue;
        }
        if (fa > f0 + c1 * alpha * dphi0 || (it > 0 && fa >= f_prev)) {
            lo = a_prev; flo = f_prev; dlo = dp_prev;
            hi = alpha;  fhi = fa;     dhi = da;
            zoom = true; break;
        }
        if (fabs(da) <= -c2 * dphi0) {
            *alpha_out = alpha; *ft = fa;           /* xt, gt already hold it */
            return true;
        }
        if (da >= 0.0) {
            lo = alpha;  flo = fa;     dlo = da;
            hi = a_prev; fhi = f_prev; dhi = dp_prev;
            zoom = true; break;
        }
        a_prev = alpha; f_prev = fa; dp_prev = da;
        alpha *= 2.0;
    }
    if (zoom) {
        for (int it = 0; it < 30; it++) {
            double aj = fm_cg_interp(lo, flo, dlo, hi, fhi, dhi);
            if (fabs(hi - lo) < 1e-16 * (1.0 + fabs(lo))) break;
            double fj = 0.0, dj = 0.0; bool okv;
            FM_CG_TRY(aj, fj, dj, okv);
            if (!okv) { hi = aj; fhi = HUGE_VAL; dhi = 0.0; continue; }
            if (fj > f0 + c1 * aj * dphi0 || fj >= flo) {
                hi = aj; fhi = fj; dhi = dj;
            } else {
                if (fabs(dj) <= -c2 * dphi0) {
                    *alpha_out = aj; *ft = fj;
                    return true;
                }
                if (dj * (hi - lo) >= 0.0) { hi = lo; fhi = flo; dhi = dlo; }
                lo = aj; flo = fj; dlo = dj;
            }
        }
    }
#undef FM_CG_TRY
    if (best_a > 0.0) {
        for (size_t i = 0; i < n; i++) { xt[i] = xbest[i]; gt[i] = gbest[i]; }
        *alpha_out = best_a; *ft = best_f;
        return true;
    }
    return false;
}

bool fm_run_cg(Expr* f, Expr** vars, size_t n,
                      FmVarBind* binds, Expr** g_exprs,
                      double* x,
                      const FmGenCon* gens, size_t ngens, double mu,
                      const FmBox* boxes,
                      const FmOpts* opts,
                      double* fx_out) {
    (void)vars;
    double* g = (double*)malloc(sizeof(double) * n);
    double* g_new = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* xb = (double*)malloc(sizeof(double) * n);   /* Wolfe scratch */
    double* gb = (double*)malloc(sizeof(double) * n);
    bool ok = false;

    if (boxes) fm_project_box(x, n, boxes);
    double fx;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    bool any_box = false;
    if (boxes)
        for (size_t i = 0; i < n; i++)
            if (boxes[i].has_lo || boxes[i].has_hi) { any_box = true; break; }
    bool use_wolfe = !augmented && !any_box;
    bool* act = any_box ? (bool*)calloc(n ? n : 1, sizeof(bool)) : NULL;
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
        fm_warn(g_fm_name, "nlnum", "gradient failed at start point");
        goto cleanup;
    }
    for (size_t i = 0; i < n; i++) d[i] = -g[i];

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    bool converged = false;
    bool stopped = false;          /* line-search / gradient failure: warned */
    double alpha_last = 0.0, gd_last = 0.0;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* Binding box bounds are held fixed (d_A = 0) and left out of the
         * gradient test, as in fm_run_bfgs: a projected step that points into
         * an active bound otherwise clips to nothing and stalls the solve. */
        size_t nact = act ? fm_box_binding_mask(x, g, n, boxes, act) : 0;
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) if (!nact || !act[i]) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { converged = true; break; }
        if (nact) for (size_t i = 0; i < n; i++) if (act[i]) d[i] = 0.0;

        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) {
            /* Restart with steepest descent. */
            for (size_t i = 0; i < n; i++) d[i] = (nact && act[i]) ? 0.0 : -g[i];
            g_dot_d = 0.0; for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        }
        double alpha, fx_new;
        bool ls_ok;
        bool have_g_new = false;
        if (use_wolfe) {
            double dnorm = sqrt(fm_dot(d, d, n));
            double a0 = (k == 0 || alpha_last <= 0.0)
                ? ((dnorm > 1.0) ? 1.0 / dnorm : 1.0)
                : alpha_last * gd_last / g_dot_d;           /* N&W (3.60) */
            if (!(a0 > 0.0) || !isfinite(a0)) a0 = (dnorm > 1.0) ? 1.0 / dnorm : 1.0;
            ls_ok = fm_cg_wolfe(f, binds, g_exprs, n, x, d, fx, g_dot_d, a0, opts,
                                x_new, g_new, &fx_new, &alpha, xb, gb);
            have_g_new = ls_ok;
            if (ls_ok) { alpha_last = alpha; gd_last = g_dot_d; }
        } else {
            ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                   augmented ? gens : NULL,
                                   augmented ? ngens : 0,
                                   augmented ? mu : 0.0,
                                   boxes, opts, &alpha, &fx_new, x_new);
        }
        if (!ls_ok) {
            if (!augmented) fm_warn(g_fm_name, "lstol", "line search failed");
            stopped = true;
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        bool ng_ok = have_g_new;
        if (!ng_ok) {
            if (augmented) {
                ng_ok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                             binds, x_new, n, opts, g_new);
            } else {
                ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
                if (!ng_ok) ng_ok = fm_grad_finite_diff(f, binds, x_new, n, opts, g_new);
            }
        }
        if (!ng_ok) {
            for (size_t i = 0; i < n; i++) x[i] = x_new[i];
            fx = fx_new;
            stopped = true;
            break;
        }
        /* Polak-Ribière+. */
        double num = 0.0, den = 0.0, gg = 0.0, gn2 = 0.0;
        for (size_t i = 0; i < n; i++) {
            num += g_new[i] * (g_new[i] - g[i]);
            den += g[i] * g[i];
            gg  += g_new[i] * g[i];
            gn2 += g_new[i] * g_new[i];
        }
        double beta = (den > 0.0) ? num / den : 0.0;
        if (beta < 0.0) beta = 0.0;
        if (use_wolfe) {
            /* Powell's restart: successive gradients far from orthogonal. */
            if (fabs(gg) >= 0.2 * gn2) beta = 0.0;
        } else if ((k + 1) % n == 0) {
            beta = 0.0;
        }
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
        }
        for (size_t i = 0; i < n; i++) {
            d[i] = -g_new[i] + beta * d[i];
            x[i] = x_new[i];
            g[i] = g_new[i];
        }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { converged = true; break; }
    }
    if (!converged && !stopped && !augmented) fm_warn_maxit(opts);
    *fx_out = fx;
    ok = true;
cleanup:
    free(g); free(g_new); free(d); free(x_new); free(xb); free(gb); free(act);
    return ok;
}
