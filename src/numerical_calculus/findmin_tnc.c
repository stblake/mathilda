/* findmin_tnc.c — truncated-Newton (TNC) local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Truncated Newton with bound Constraints (TNC, Hessian-free)         *
 * ------------------------------------------------------------------ *
 * Nash's truncated (inexact) Newton method: each outer iteration
 * approximately solves the Newton system H*p = -g by an INNER
 * conjugate-gradient loop that accesses the Hessian only through
 * Hessian-vector products Hv, computed by finite-differencing the
 * (exact, compiled) gradient -- so the Hessian is never formed. Box
 * bounds use the same active-set projection as fm_run_lbfgsb (masked
 * gradient, projected-gradient KKT test, free-subspace direction,
 * alpha_max face cap, strong-Wolfe line search). The niche vs the other
 * methods: true Newton curvature at O(n) memory, so it scales to large n
 * where "Newton" would form an O(n^2) symbolic Hessian, while converging
 * in fewer outer iterations than L-BFGS's low-rank model on
 * ill-conditioned problems -- Mathilda's exact compiled gradient makes
 * the finite-difference Hv both cheap and accurate. Machine precision
 * only (WorkingPrecision > MachinePrecision falls back to QuasiNewton).
 * Exposed as Method -> "TNC" (alias "TruncatedNewton"). Refs: Nash 1984;
 * Dembo & Steihaug 1983 (truncation); Nocedal & Wright 2nd ed. ch. 7. */

/* Gradient (only) at point p: plain or augmented, exact-then-FD. The
 * gradient half of fm_lbfgs_fg, for the two evaluations of a
 * Hessian-vector product (the objective value is not needed there). */
bool fm_tnc_grad(const FmLbfgsCtx* c, const double* p, double* g) {
    if (c->augmented)
        return fm_eval_aug_gradient(c->f, c->g_exprs, c->gens, c->ngens, c->mu,
                                    c->binds, p, c->n, c->opts, g);
    bool gok = c->g_exprs && fm_eval_gradient(c->g_exprs, c->binds, p, c->n, c->opts, g);
    if (!gok) gok = fm_grad_finite_diff(c->f, c->binds, p, c->n, c->opts, g);
    return gok;
}

/* Hessian-vector product Hv ~= (grad f(x + h v) - grad f(x)) / h, restricted
 * to the free subspace (v is zero on active coords; the result is masked to
 * zero there so the inner CG stays in the free subspace, realising P H P).
 * g_base = grad f(x) is cached by the caller. h fixes the perturbation NORM
 * at sqrt(eps)*max(1,||x||) regardless of ||v||: sqrt(eps) is the optimal
 * forward-difference step for a gradient whose noise is ~eps (the compiled
 * gradient is exact). Returns false on a non-finite perturbed gradient. */
bool fm_tnc_hessvec(const FmLbfgsCtx* c, const double* x, const double* g_base,
                           const double* v, const bool* active,
                           double* xpert, double* gpert, double* Hv) {
    size_t n = c->n;
    double vn = 0.0, xn = 0.0;
    for (size_t i = 0; i < n; i++) { vn += v[i] * v[i]; xn += x[i] * x[i]; }
    vn = sqrt(vn); xn = sqrt(xn);
    if (vn < 1e-300) { for (size_t i = 0; i < n; i++) Hv[i] = 0.0; return true; }
    const double sqrt_eps = 1.4901161193847656e-08;   /* sqrt(2^-52) */
    double h = sqrt_eps * (xn > 1.0 ? xn : 1.0) / vn;
    for (size_t i = 0; i < n; i++) xpert[i] = x[i] + h * v[i];
    if (!fm_tnc_grad(c, xpert, gpert)) return false;
    for (size_t i = 0; i < n; i++) {
        Hv[i] = active[i] ? 0.0 : (gpert[i] - g_base[i]) / h;
        if (!isfinite(Hv[i])) return false;
    }
    return true;
}

/* Inner truncated conjugate-gradient: approximately solve H p = -gm over the
 * free subspace (gm is the masked gradient; p, r, d stay zero on active coords
 * throughout). Truncates on (1) non-positive curvature, (2) the inexact-Newton
 * forcing sequence ||r|| <= eta ||g|| with eta = min(0.5, sqrt(||g||)), or
 * (3) maxcg iterations. The returned p is always a descent direction (g.p < 0):
 * CG from p=0 on the SPD Newton system with rhs -gm gives -gm.p_k > 0. scratch
 * r,d,Hd,xpert,gpert are caller-owned (length n). */
void fm_tnc_cg(const FmLbfgsCtx* c, const double* x, const double* g_base,
                      const double* gm, const bool* active,
                      double* p, double* r, double* d, double* Hd,
                      double* xpert, double* gpert) {
    size_t n = c->n;
    const double curv_eps = 2.220446049250313e-16;
    for (size_t i = 0; i < n; i++) { p[i] = 0.0; r[i] = -gm[i]; d[i] = r[i]; }
    double rr = 0.0;
    for (size_t i = 0; i < n; i++) rr += r[i] * r[i];
    double gm_norm = sqrt(rr);
    if (gm_norm == 0.0) return;                        /* p = 0 (KKT / all active) */
    double eta = sqrt(gm_norm);
    if (eta > 0.5) eta = 0.5;                          /* min(0.5, sqrt(||g||)) */
    /* Cap the inner iterations at min(50, n) rather than scipy's min(50, n/2):
     * scipy's n/2 relies on its diagonal preconditioner to get a near-Newton
     * direction in few inner steps, whereas the identity-preconditioned CG here
     * must solve the Newton system more fully. For small n this reaches the exact
     * Newton step (linear CG converges in <= n steps), so a 2-D curved valley
     * (Rosenbrock) gets a true Newton direction instead of a lone
     * steepest-descent step that would zigzag. A diagonal preconditioner is the
     * v2 refinement that would let this drop back toward n/2. */
    size_t maxcg = n; if (maxcg > 50) maxcg = 50; if (maxcg < 1) maxcg = 1;

    for (size_t j = 0; j < maxcg; j++) {
        if (!fm_tnc_hessvec(c, x, g_base, d, active, xpert, gpert, Hd)) {
            if (j == 0) for (size_t i = 0; i < n; i++) p[i] = -gm[i];  /* steepest */
            return;                                    /* else keep accumulated p */
        }
        double dHd = 0.0, dd = 0.0;
        for (size_t i = 0; i < n; i++) { dHd += d[i] * Hd[i]; dd += d[i] * d[i]; }
        if (dHd <= curv_eps * dd) {                    /* non-positive curvature */
            if (j == 0) for (size_t i = 0; i < n; i++) p[i] = -gm[i];
            return;
        }
        double alpha = rr / dHd;
        for (size_t i = 0; i < n; i++) { p[i] += alpha * d[i]; r[i] -= alpha * Hd[i]; }
        double rr_new = 0.0;
        for (size_t i = 0; i < n; i++) rr_new += r[i] * r[i];
        if (sqrt(rr_new) <= eta * gm_norm) return;     /* forcing sequence met */
        double beta = rr_new / rr;                     /* Fletcher-Reeves */
        for (size_t i = 0; i < n; i++) d[i] = r[i] + beta * d[i];
        rr = rr_new;
    }
    /* maxcg reached: every curvature used was positive, so p is a descent step. */
}

bool fm_run_tnc(Expr* f, Expr** vars, size_t n,
                       FmVarBind* binds, Expr** g_exprs,
                       double* x, /* in/out */
                       const FmGenCon* gens, size_t ngens, double mu,
                       const FmBox* boxes,
                       const FmOpts* opts,
                       double* fx_out) {
    (void)vars;
    /* gm and active use calloc, not malloc: the per-iteration loop below fills all n elements
       before passing them (as const*) to fm_tnc_cg, but GCC cannot relate the trip count n to the
       reads inside fm_tnc_cg and flags a phantom uninitialized read on the degenerate n==0 path.
       Zero-init settles it at negligible once-per-solve cost. The other buffers are written inside
       the CG inner loop before use, so they do not warn. */
    double* g      = (double*)malloc(sizeof(double) * n);  /* base grad, cached  */
    double* gm     = (double*)calloc(n, sizeof(double));   /* masked gradient    */
    double* g_new  = (double*)malloc(sizeof(double) * n);
    double* p      = (double*)malloc(sizeof(double) * n);  /* CG step / direction*/
    double* r      = (double*)malloc(sizeof(double) * n);  /* CG residual        */
    double* d      = (double*)malloc(sizeof(double) * n);  /* CG direction       */
    double* Hd     = (double*)malloc(sizeof(double) * n);  /* Hessian-vector prod*/
    double* xpert  = (double*)malloc(sizeof(double) * n);  /* Hv scratch         */
    double* gpert  = (double*)malloc(sizeof(double) * n);  /* Hv scratch         */
    double* x_new  = (double*)malloc(sizeof(double) * n);
    bool*   active = (bool*)calloc(n, sizeof(bool));
    bool ok = false;
    double fx = 0.0;
    if (!g || !gm || !g_new || !p || !r || !d || !Hd || !xpert || !gpert
        || !x_new || !active)
        goto cleanup;

    FmLbfgsCtx ctx;
    ctx.f = f; ctx.g_exprs = g_exprs; ctx.binds = binds; ctx.n = n;
    ctx.gens = gens; ctx.ngens = ngens; ctx.mu = mu; ctx.boxes = boxes;
    ctx.opts = opts; ctx.augmented = (mu > 0.0 && gens && ngens > 0);

    if (boxes) fm_project_box(x, n, boxes);
    if (!fm_lbfgs_fg(&ctx, x, &fx, g)) {
        fm_warn(g_fm_name, "nlnum", "objective/gradient evaluation failed at start point");
        goto cleanup;
    }

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* Active set + masked gradient + projected-gradient KKT norm
         * (identical to fm_run_lbfgsb). */
        double pgnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xi = x[i] - g[i];
            active[i] = false;
            if (boxes) {
                double lt = (boxes[i].has_lo) ? 1e-12 * (1.0 + fabs(boxes[i].lo)) : 0.0;
                double ut = (boxes[i].has_hi) ? 1e-12 * (1.0 + fabs(boxes[i].hi)) : 0.0;
                if (boxes[i].has_lo && xi < boxes[i].lo) xi = boxes[i].lo;
                if (boxes[i].has_hi && xi > boxes[i].hi) xi = boxes[i].hi;
                if ((boxes[i].has_lo && x[i] <= boxes[i].lo + lt && g[i] > 0.0) ||
                    (boxes[i].has_hi && x[i] >= boxes[i].hi - ut && g[i] < 0.0))
                    active[i] = true;
            }
            double pg = fabs(x[i] - xi);
            if (pg > pgnorm) pgnorm = pg;
            gm[i] = active[i] ? 0.0 : g[i];
        }
        if (pgnorm < tol_acc) { ok = true; break; }

        /* Inner truncated-CG for the Hessian-free Newton direction p. */
        fm_tnc_cg(&ctx, x, g, gm, active, p, r, d, Hd, xpert, gpert);

        /* Keep p in the feasible cone (identical clamp to fm_run_lbfgsb). */
        if (boxes) {
            for (size_t i = 0; i < n; i++) {
                if (active[i]) { p[i] = 0.0; continue; }
                double lt = (boxes[i].has_lo) ? 1e-12 * (1.0 + fabs(boxes[i].lo)) : 0.0;
                double ut = (boxes[i].has_hi) ? 1e-12 * (1.0 + fabs(boxes[i].hi)) : 0.0;
                if (boxes[i].has_hi && x[i] >= boxes[i].hi - ut && p[i] > 0.0) p[i] = 0.0;
                if (boxes[i].has_lo && x[i] <= boxes[i].lo + lt && p[i] < 0.0) p[i] = 0.0;
            }
        }
        double dphi0 = 0.0;
        for (size_t i = 0; i < n; i++) dphi0 += g[i] * p[i];
        if (dphi0 >= 0.0) {
            /* FD noise let an indefinite direction through: masked steepest. */
            for (size_t i = 0; i < n; i++) p[i] = -gm[i];
            dphi0 = 0.0; for (size_t i = 0; i < n; i++) dphi0 += g[i] * p[i];
            if (dphi0 >= 0.0) { ok = true; break; }    /* no feasible descent → KKT */
        }

        /* alpha_max: nearest box face along p (identical to fm_run_lbfgsb). */
        double alpha_max = HUGE_VAL;
        if (boxes) {
            for (size_t i = 0; i < n; i++) {
                if (p[i] > 0.0 && boxes[i].has_hi) {
                    double amx = (boxes[i].hi - x[i]) / p[i];
                    if (amx < alpha_max) alpha_max = amx;
                } else if (p[i] < 0.0 && boxes[i].has_lo) {
                    double amx = (boxes[i].lo - x[i]) / p[i];
                    if (amx < alpha_max) alpha_max = amx;
                }
            }
            if (alpha_max < 0.0) alpha_max = 0.0;
        }
        if (alpha_max <= 0.0) { ok = true; break; }    /* pinned at a corner */

        double a, fx_new;
        if (!fm_lbfgs_linesearch(&ctx, x, p, fx, dphi0, alpha_max,
                                 x_new, g_new, &fx_new, &a)) {
            if (!ctx.augmented)
                fm_warn(g_fm_name, "lstol", "line search failed at iter %lld",
                        (long long)k);
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        /* Convergence: the projected-gradient KKT test (pgnorm, above) is the
         * primary criterion; here a relative function-value stall detector
         * guards against spinning at the machine-noise floor. A step-SIZE test
         * is deliberately NOT used -- on a narrow curved valley (Rosenbrock) a
         * legitimate small step trips it far from the optimum, the same trap
         * fm_run_lbfgsb documents. The line search left g_new = grad(x_new);
         * reuse it as the next base gradient (no extra eval). TNC keeps no
         * curvature memory between iterations. */
        double fdenom = fabs(fx);
        if (fabs(fx_new) > fdenom) fdenom = fabs(fx_new);
        if (fdenom < 1.0) fdenom = 1.0;
        double f_rel = (fx - fx_new) / fdenom;
        for (size_t i = 0; i < n; i++) { x[i] = x_new[i]; g[i] = g_new[i]; }
        fx = fx_new;
        if (f_rel >= 0.0 && f_rel < 1e-14) { ok = true; break; }
    }
    *fx_out = fx;
    ok = true;
cleanup:
    free(g); free(gm); free(g_new); free(p); free(r); free(d); free(Hd);
    free(xpert); free(gpert); free(x_new); free(active);
    return ok;
}
