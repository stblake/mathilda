/* findmin_lbfgsb.c — L-BFGS-B bound-constrained local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  L-BFGS-B — limited-memory BFGS with bound constraints              *
 * ------------------------------------------------------------------ *
 *
 * A Mathilda extension to FindMinimum's Method set (Mathematica exposes no
 * such method name), selected by Method -> "LBFGSB" (aliases "LBFGS",
 * "LimitedMemoryBFGS"). Unlike the full-memory QuasiNewton solver
 * (fm_run_bfgs), which stores a dense n×n inverse-Hessian approximation and is
 * therefore O(n^2) in memory and per-iteration cost, this keeps only the m most
 * recent correction pairs and forms the search direction by the Nocedal
 * two-loop recursion at O(m·n) — so it scales to n in the thousands where the
 * dense method does not. Reference: Byrd, Lu, Nocedal & Zhu 1995, with the
 * Morales–Nocedal 2011 correction (the variant scipy's L-BFGS-B ships).
 *
 * This is milestone M1 (limited-memory core): the direction is the
 * unconstrained two-loop step and bound constraints are honoured by the
 * projecting line search (fm_line_search projects each trial onto the box),
 * exactly as fm_run_bfgs handles boxes today. When no bound is active this is
 * L-BFGS-B; the true generalized-Cauchy-point / subspace-minimization path
 * (M2) refines the bound-active case. The augmented-objective branch (mu > 0)
 * makes L-BFGS-B usable as the inner solver of the augmented-Lagrangian outer
 * loop for general constraints (M3, via fm_run_penalty). */

#ifndef FM_LBFGS_DEFAULT_M
#define FM_LBFGS_DEFAULT_M 10   /* history depth (scipy L-BFGS-B maxcor) */
#endif

/* d = -H_k · g by the Nocedal two-loop recursion. `alpha` (length m) and `q`
 * (length n) are caller-provided scratch. With an empty memory (cnt == 0) this
 * degrades to d = -gamma·g (scaled steepest descent), so it is always valid. */
static void fm_lbfgs_direction(const FmLbfgsMem* mem, const double* g,
                               double* alpha, double* q, double* d) {
    size_t n = mem->n, m = mem->m;
    for (size_t i = 0; i < n; i++) q[i] = g[i];
    /* First loop, newest → oldest. */
    for (size_t t = 0; t < mem->cnt; t++) {
        size_t j = (mem->head + m - 1 - t) % m;   /* t-th newest pair */
        const double* s = mem->S + j * n;
        const double* y = mem->Y + j * n;
        double sq = 0.0;
        for (size_t i = 0; i < n; i++) sq += s[i] * q[i];
        double a = mem->rho[j] * sq;
        alpha[j] = a;
        for (size_t i = 0; i < n; i++) q[i] -= a * y[i];
    }
    /* r = H0 q, H0 = gamma·I (gamma == 1 before any pair is stored). */
    double g0 = (mem->cnt > 0) ? mem->gamma : 1.0;
    for (size_t i = 0; i < n; i++) d[i] = g0 * q[i];
    /* Second loop, oldest → newest. */
    for (size_t t = 0; t < mem->cnt; t++) {
        size_t j = (mem->head + m - mem->cnt + t) % m;   /* t-th oldest pair */
        const double* s = mem->S + j * n;
        const double* y = mem->Y + j * n;
        double yr = 0.0;
        for (size_t i = 0; i < n; i++) yr += y[i] * d[i];
        double coef = alpha[j] - mem->rho[j] * yr;
        for (size_t i = 0; i < n; i++) d[i] += coef * s[i];
    }
    for (size_t i = 0; i < n; i++) d[i] = -d[i];
}

/* Insert a correction pair (s = x_new - x, y = g_new - g), skipping it when the
 * curvature s·y is not sufficiently positive (keeps H0 and the stored pairs
 * positive-definite — the standard L-BFGS damping-free skip). */
static void fm_lbfgs_push(FmLbfgsMem* mem, const double* s, const double* y) {
    size_t n = mem->n;
    double sy = 0.0, yy = 0.0;
    for (size_t i = 0; i < n; i++) { sy += s[i] * y[i]; yy += y[i] * y[i]; }
    const double curv_eps = 2.220446049250313e-16;   /* machine epsilon */
    if (yy <= 0.0 || sy <= curv_eps * yy) return;    /* skip, keep memory */
    size_t j = mem->head;
    double* Sj = mem->S + j * n;
    double* Yj = mem->Y + j * n;
    for (size_t i = 0; i < n; i++) { Sj[i] = s[i]; Yj[i] = y[i]; }
    mem->rho[j] = 1.0 / sy;
    mem->gamma  = sy / yy;
    mem->head = (mem->head + 1) % mem->m;
    if (mem->cnt < mem->m) mem->cnt++;
}

/* Objective value AND gradient at point p (length n): plain or augmented. */
bool fm_lbfgs_fg(const FmLbfgsCtx* c, const double* p, double* fval, double* g) {
    bool ok = c->augmented
        ? fm_eval_augmented(c->f, c->binds, p, c->n, c->gens, c->ngens, c->mu, c->opts, fval)
        : fm_eval_scalar(c->f, c->binds, p, c->n, c->opts, fval);
    if (!ok) return false;
    if (c->augmented)
        return fm_eval_aug_gradient(c->f, c->g_exprs, c->gens, c->ngens, c->mu,
                                    c->binds, p, c->n, c->opts, g);
    bool gok = c->g_exprs && fm_eval_gradient(c->g_exprs, c->binds, p, c->n, c->opts, g);
    if (!gok) gok = fm_grad_finite_diff(c->f, c->binds, p, c->n, c->opts, g);
    return gok;
}

/* phi(alpha) = f(x + alpha d) and its slope phi'(alpha) = grad(x+alpha d)·d.
 * xa/ga receive the trial point and its gradient. The box projection is a
 * no-op while alpha <= alpha_max (the search is bounded there), so the slope
 * is exact. */
static bool fm_lbfgs_phi(const FmLbfgsCtx* c, const double* x, const double* d,
                         double alpha, double* xa, double* ga,
                         double* phi, double* dphi) {
    size_t n = c->n;
    for (size_t i = 0; i < n; i++) xa[i] = x[i] + alpha * d[i];
    if (c->boxes) fm_project_box(xa, n, c->boxes);
    if (!fm_lbfgs_fg(c, xa, phi, ga)) return false;
    double sdot = 0.0;
    for (size_t i = 0; i < n; i++) sdot += ga[i] * d[i];
    *dphi = sdot;
    return true;
}

/* "zoom" between a bracketing pair (Nocedal & Wright, Alg. 3.6). a_lo has the
 * lower function value and satisfies Armijo; the strong-Wolfe minimiser lies
 * between a_lo and a_hi. Always returns a usable point (the bracket's low end
 * if the curvature condition is not met within the budget), so the search never
 * aborts a solve. On return xa/ga hold the accepted point and its gradient. */
static bool fm_lbfgs_zoom(const FmLbfgsCtx* c, const double* x, const double* d,
                          double f0, double dphi0, double c1, double c2,
                          double a_lo, double f_lo, double a_hi,
                          double* xa, double* ga, double* f_out, double* a_out) {
    for (int it = 0; it < 30; it++) {
        double a_j = 0.5 * (a_lo + a_hi);           /* bisection: always safe */
        double phi, dphi;
        if (!fm_lbfgs_phi(c, x, d, a_j, xa, ga, &phi, &dphi)) { a_hi = a_j; continue; }
        if (phi > f0 + c1 * a_j * dphi0 || phi >= f_lo) {
            a_hi = a_j;
        } else {
            if (fabs(dphi) <= -c2 * dphi0) { *f_out = phi; *a_out = a_j; return true; }
            if (dphi * (a_hi - a_lo) >= 0.0) a_hi = a_lo;
            a_lo = a_j; f_lo = phi;
        }
        if (fabs(a_hi - a_lo) < 1e-18) break;
    }
    /* Curvature not met in budget: settle on the low bracket point, which is a
     * genuine sufficient-decrease step (so the outer loop keeps progressing). */
    double phi, dphi;
    if (a_lo <= 0.0 || !fm_lbfgs_phi(c, x, d, a_lo, xa, ga, &phi, &dphi)) return false;
    *f_out = phi; *a_out = a_lo;
    return true;
}

/* Line search for the L-BFGS direction: the quasi-Newton UNIT step (alpha = 1)
 * is tried FIRST and the step is EXPANDED (alpha doubled) while the objective
 * keeps decreasing steeply — essential for following a curved valley (extended
 * Rosenbrock), and the reason the shared fm_line_search (which only backtracks
 * from a 1/||d||-capped step) stalls here. Implements the strong-Wolfe
 * bracketing of Nocedal & Wright Alg. 3.5, but is made robust: it falls back to
 * the best sufficient-decrease point rather than ever aborting the solve.
 * Restricted to [0, alpha_max]; a step reaching alpha_max lands on a box face
 * (that coordinate joins the active set next iteration). On success xa/ga hold
 * the accepted point and its gradient. */
bool fm_lbfgs_linesearch(const FmLbfgsCtx* c, const double* x, const double* d,
                                double f0, double dphi0, double alpha_max,
                                double* xa, double* ga, double* f_out, double* a_out) {
    if (dphi0 >= 0.0) return false;             /* not a descent direction */
    const double c1 = 1e-4, c2 = 0.9;
    double a_prev = 0.0, f_prev = f0;
    double a_cur = (alpha_max < 1.0) ? alpha_max : 1.0;
    if (a_cur <= 0.0) return false;
    bool have_prev_ok = false;                  /* a_prev is a valid Armijo point */
    for (int it = 0; it < 40; it++) {
        double phi, dphi;
        if (!fm_lbfgs_phi(c, x, d, a_cur, xa, ga, &phi, &dphi)) {
            a_cur = 0.5 * (a_prev + a_cur);     /* domain/non-finite: shrink */
            if (a_cur - a_prev < 1e-18) break;
            continue;
        }
        if (phi > f0 + c1 * a_cur * dphi0 || (it > 0 && phi >= f_prev))
            return fm_lbfgs_zoom(c, x, d, f0, dphi0, c1, c2,
                                 a_prev, f_prev, a_cur, xa, ga, f_out, a_out);
        if (fabs(dphi) <= -c2 * dphi0) { *f_out = phi; *a_out = a_cur; return true; }
        if (dphi >= 0.0)
            return fm_lbfgs_zoom(c, x, d, f0, dphi0, c1, c2,
                                 a_cur, phi, a_prev, xa, ga, f_out, a_out);
        /* Armijo holds and the slope is still steeply negative: expand. */
        if (a_cur >= alpha_max) { *f_out = phi; *a_out = a_cur; return true; }  /* on a face */
        a_prev = a_cur; f_prev = phi; have_prev_ok = true;
        a_cur = (2.0 * a_cur < alpha_max) ? 2.0 * a_cur : alpha_max;
    }
    /* Budget exhausted: fall back to the last accepted Armijo point. */
    if (have_prev_ok) {
        double phi, dphi;
        if (fm_lbfgs_phi(c, x, d, a_prev, xa, ga, &phi, &dphi)) {
            *f_out = phi; *a_out = a_prev; return true;
        }
    }
    return false;
}

bool fm_run_lbfgsb(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars;
    size_t m = FM_LBFGS_DEFAULT_M;
    FmLbfgsMem mem;
    mem.n = n; mem.m = m; mem.cnt = 0; mem.head = 0; mem.gamma = 1.0;
    mem.S   = (double*)malloc(sizeof(double) * m * n);
    mem.Y   = (double*)malloc(sizeof(double) * m * n);
    mem.rho = (double*)malloc(sizeof(double) * m);
    double* alpha = (double*)malloc(sizeof(double) * m);
    double* g     = (double*)malloc(sizeof(double) * n);
    double* gm    = (double*)malloc(sizeof(double) * n);   /* masked gradient */
    double* g_new = (double*)malloc(sizeof(double) * n);
    double* d     = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* s     = (double*)malloc(sizeof(double) * n);
    double* y     = (double*)malloc(sizeof(double) * n);
    double* q     = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    double fx = 0.0;
    if (!mem.S || !mem.Y || !mem.rho || !alpha || !g || !gm || !g_new || !d
        || !x_new || !s || !y || !q)
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

    double tol_acc = pow(10.0, -opts->acc_goal_digits);   /* projected-grad tol */

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* Projected-gradient infinity-norm convergence, and the active-set /
         * masked gradient. A coordinate at a box face with the gradient
         * pushing outward is "active": it satisfies its KKT condition, so it
         * contributes zero to the projected-gradient norm and is masked to
         * zero for the direction computation. */
        double pgnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xi = x[i] - g[i];
            bool active = false;
            if (boxes) {
                double lt = (boxes[i].has_lo) ? 1e-12 * (1.0 + fabs(boxes[i].lo)) : 0.0;
                double ut = (boxes[i].has_hi) ? 1e-12 * (1.0 + fabs(boxes[i].hi)) : 0.0;
                if (boxes[i].has_lo && xi < boxes[i].lo) xi = boxes[i].lo;
                if (boxes[i].has_hi && xi > boxes[i].hi) xi = boxes[i].hi;
                if ((boxes[i].has_lo && x[i] <= boxes[i].lo + lt && g[i] > 0.0) ||
                    (boxes[i].has_hi && x[i] >= boxes[i].hi - ut && g[i] < 0.0))
                    active = true;
            }
            double pg = fabs(x[i] - xi);
            if (pg > pgnorm) pgnorm = pg;
            gm[i] = active ? 0.0 : g[i];
        }
        if (pgnorm < tol_acc) { ok = true; break; }

        /* Search direction d = -H_k g on the free subspace. */
        fm_lbfgs_direction(&mem, gm, alpha, q, d);
        /* Keep d in the feasible cone: no motion on a masked (active) face, and
         * no component that would immediately drive a coordinate through a
         * bound it already sits on. */
        if (boxes) {
            for (size_t i = 0; i < n; i++) {
                if (gm[i] == 0.0 && g[i] != 0.0) { d[i] = 0.0; continue; }
                double lt = (boxes[i].has_lo) ? 1e-12 * (1.0 + fabs(boxes[i].lo)) : 0.0;
                double ut = (boxes[i].has_hi) ? 1e-12 * (1.0 + fabs(boxes[i].hi)) : 0.0;
                if (boxes[i].has_hi && x[i] >= boxes[i].hi - ut && d[i] > 0.0) d[i] = 0.0;
                if (boxes[i].has_lo && x[i] <= boxes[i].lo + lt && d[i] < 0.0) d[i] = 0.0;
            }
        }
        double dphi0 = 0.0;
        for (size_t i = 0; i < n; i++) dphi0 += g[i] * d[i];
        if (dphi0 >= 0.0) {
            /* Non-descent (stale memory): reset and take a projected steepest
             * step on the free coordinates. */
            mem.cnt = 0; mem.head = 0; mem.gamma = 1.0;
            for (size_t i = 0; i < n; i++) d[i] = -gm[i];
            dphi0 = 0.0; for (size_t i = 0; i < n; i++) dphi0 += g[i] * d[i];
            if (dphi0 >= 0.0) { ok = true; break; }   /* no feasible descent → KKT */
        }

        /* Cap the step at the nearest box face reached along d. */
        double alpha_max = HUGE_VAL;
        if (boxes) {
            for (size_t i = 0; i < n; i++) {
                if (d[i] > 0.0 && boxes[i].has_hi) {
                    double amx = (boxes[i].hi - x[i]) / d[i];
                    if (amx < alpha_max) alpha_max = amx;
                } else if (d[i] < 0.0 && boxes[i].has_lo) {
                    double amx = (boxes[i].lo - x[i]) / d[i];
                    if (amx < alpha_max) alpha_max = amx;
                }
            }
            if (alpha_max < 0.0) alpha_max = 0.0;
        }
        if (alpha_max <= 0.0) { ok = true; break; }   /* pinned at a corner */

        double a, fx_new;
        if (!fm_lbfgs_linesearch(&ctx, x, d, fx, dphi0, alpha_max,
                                 x_new, g_new, &fx_new, &a)) {
            /* Silent at high mu in the penalty schedule — fm_run_penalty's
             * feasibility check is the authoritative signal there. */
            if (!ctx.augmented)
                fm_warn(g_fm_name, "lstol", "line search failed at iter %lld",
                        (long long)k);
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        /* Relative function-value change (cf. scipy L-BFGS-B's "factr" test).
         * A conservative threshold is used purely as a stall detector — the
         * projected-gradient test above is the primary convergence criterion.
         * A single-step DISPLACEMENT test was tried here first and rejected: on
         * a narrow curved valley (extended Rosenbrock at large n) a legitimate
         * small step would trip it and declare convergence far from the
         * optimum. */
        double fdenom = fabs(fx);
        if (fabs(fx_new) > fdenom) fdenom = fabs(fx_new);
        if (fdenom < 1.0) fdenom = 1.0;
        double f_rel = (fx - fx_new) / fdenom;

        /* fm_lbfgs_linesearch already left g_new = grad(x_new). Store the
         * correction pair (curvature-guarded). */
        for (size_t i = 0; i < n; i++) { s[i] = x_new[i] - x[i]; y[i] = g_new[i] - g[i]; }
        fm_lbfgs_push(&mem, s, y);

        for (size_t i = 0; i < n; i++) { x[i] = x_new[i]; g[i] = g_new[i]; }
        fx = fx_new;
        /* Genuine numerical stall (relative change at the machine-noise floor).
         * Deliberately far tighter than scipy's factr so it cannot pre-empt the
         * projected-gradient test on a slowly-improving valley. */
        if (f_rel >= 0.0 && f_rel < 1e-14) { ok = true; break; }
    }
    *fx_out = fx;
    ok = true;
cleanup:
    free(mem.S); free(mem.Y); free(mem.rho); free(alpha);
    free(g); free(gm); free(g_new); free(d); free(x_new); free(s); free(y); free(q);
    return ok;
}
