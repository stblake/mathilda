/* findmin_slsqp.c — SLSQP sequential quadratic programming local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  SLSQP -- Sequential Least-Squares Quadratic Programming             *
 * ------------------------------------------------------------------ *
 * Han-Powell SQP for smooth CONSTRAINED problems (equality + inequality
 * + bounds), the faithful analogue of scipy's minimize(method="SLSQP").
 * Each outer iteration replaces f and the constraints by a quadratic /
 * linear model at the current iterate and solves the QP
 *     min_d  1/2 dᵀB d + ∇fᵀd
 *     s.t.   ∇h_jᵀd + h_j = 0            (equalities)
 *            ∇g_iᵀd + g_i ≤ 0            (inequalities, feasible ≡ g ≤ 0)
 *            lo ≤ x + d ≤ hi             (bounds)
 * where B is a BFGS approximation to the Hessian of the LAGRANGIAN kept
 * SPD by Powell's (1978) damping.  The step is accepted by an Armijo
 * backtracking line search on the L1 exact-penalty merit function, whose
 * penalties are driven above the QP multipliers by Powell's rule so d is
 * always a descent direction.  Unlike the augmented-Lagrangian penalty
 * wrapper (fm_run_penalty) used by the other gradient methods, SLSQP
 * treats the constraints DIRECTLY through the QP, giving super-linear
 * local convergence and accurate constraint satisfaction.  With no
 * constraints and no bounds it degrades gracefully to a damped-BFGS
 * line-search solve.  Machine precision only (WorkingPrecision >
 * MachinePrecision falls back to QuasiNewton in the driver).  Exposed as
 * Method -> "SLSQP" (alias "SequentialQuadraticProgramming").  Refs:
 * Kraft 1988; Powell 1978; Nocedal & Wright 2nd ed. ch. 18; Goldfarb &
 * Idnani 1983 (the dual active-set QP). */

/* Objective gradient at x: exact symbolic (compiled) when available, else
 * finite differences.  A black-box objective (a ?NumericQ helper, a
 * LinearSolve inside a Module) has no symbolic gradient, and its differences
 * ARE the cost of the method: n + 1 samples an iteration forward, 2n central.
 * Far from the solution the forward difference (step sqrt(eps)|x|, error
 * ~1e-8 relative) is as good a search direction as the central one at half
 * the price -- the choice scipy's SLSQP makes throughout.  Once the QP steps
 * shrink (`central`, set by the caller and never cleared) the final
 * iterations switch to central differences, whose ~1e-11 accuracy is what the
 * stationarity test then needs.  A forward step that would leave a box bound
 * is taken backwards instead.  fx is f(x), unscaled. */
static bool fm_slsqp_objgrad(Expr* f, Expr** g_exprs, FmVarBind* binds,
                             const double* x, size_t n, const FmOpts* opts,
                             double* g, double fx, const FmBox* boxes,
                             bool central, bool* had_exact) {
    *had_exact = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
    if (*had_exact) return true;
    if (central) return fm_grad_finite_diff(f, binds, x, n, opts, g);
    double* xp = (double*)malloc(sizeof(double) * (n ? n : 1));
    if (!xp) return false;
    for (size_t i = 0; i < n; i++) xp[i] = x[i];
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        double xi = x[i];
        double h = 1.4901161193847656e-08 * (fabs(xi) > 1.0 ? fabs(xi) : 1.0);
        if (boxes && boxes[i].has_hi && xi + h > boxes[i].hi) h = -h;
        xp[i] = xi + h;
        h = xp[i] - xi;                      /* the step actually represented */
        double f1;
        ok = fm_eval_scalar(f, binds, xp, n, opts, &f1) && isfinite(f1);
        if (ok) g[i] = (f1 - fx) / h;
        xp[i] = xi;
    }
    free(xp);
    return ok;
}

static double fm_slsqp_infnorm(const double* v, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; i++) if (fabs(v[i]) > m) m = fabs(v[i]);
    return m;
}

/* One constraint's Jacobian row ∇c_k at x: exact-then-FD, mirroring
 * fm_eval_aug_gradient's per-constraint handling. */
static bool fm_slsqp_congrad(const FmGenCon* gk, FmVarBind* binds,
                             const double* x, size_t n, const FmOpts* opts,
                             double* row) {
    bool ok = gk->grad_exprs
              && fm_eval_gradient(gk->grad_exprs, binds, x, n, opts, row);
    if (!ok) ok = fm_grad_finite_diff(gk->expr, binds, x, n, opts, row);
    return ok;
}

/* Solve the strictly-convex QP
 *     min_d  1/2 dᵀB d + gᵀd     s.t.  n_kᵀd  =  b_k   (is_eq[k])
 *                                       n_kᵀd  ≥  b_k   (inequality/bound)
 * by a dual active-set sweep on B's Cholesky factor L (B = L Lᵀ).  A dual
 * method needs NO Phase-1 / feasible start (SQP iterates are routinely
 * infeasible w.r.t. the linearisation): it begins from the equality-only
 * active set and, each round, solves the equality-constrained KKT system
 * over the active set via the range-space identity
 *     d = Σ_{k∈A} μ_k B⁻¹n_k − B⁻¹g,
 *     K μ = rhs,  K_{jl}=n_jᵀB⁻¹n_l,  rhs_j = b_j + n_jᵀB⁻¹g,
 * then either DROPS the active inequality with the most-negative μ (dual
 * step) or ADDS the most-violated inactive inequality (primal step) until
 * KKT holds.  Returns 0 on a KKT point, 1 when a violated inequality can
 * only be met by a linearly-dependent normal (inconsistent linearisation
 * -> caller relaxes), 2 on the iteration cap (d is still usable).  Fills
 * d (nv) and mult (m, 0 on inactive; μ_k ≥ 0 for active inequalities). */
int fm_slsqp_activeset(const double* L, size_t nv, const double* g,
                              const double* Nrm, const double* b,
                              const int* is_eq, size_t m,
                              double* d, double* mult) {
    int status = 2;
    size_t m1 = m ? m : 1;
    double* w_g = (double*)malloc(sizeof(double) * nv);
    double* W   = (double*)malloc(sizeof(double) * m1 * nv);
    bool*   act = (bool*)malloc(sizeof(bool) * m1);
    bool*   blk = (bool*)calloc(m1, sizeof(bool));   /* blocked (dependent) */
    size_t* A   = (size_t*)malloc(sizeof(size_t) * m1);
    double* K   = (double*)malloc(sizeof(double) * m1 * m1);
    double* Kc  = (double*)malloc(sizeof(double) * m1 * m1);
    double* rhs = (double*)malloc(sizeof(double) * m1);
    double* muA = (double*)malloc(sizeof(double) * m1);
    if (!w_g || !W || !act || !blk || !A || !K || !Kc || !rhs || !muA) {
        status = 1; goto done;
    }
    fm_chol_solve(L, nv, g, w_g);                    /* w_g = B⁻¹ g */
    for (size_t k = 0; k < m; k++) act[k] = (is_eq[k] != 0);

    const double drop_tol = -1e-9;
    const double feas_tol =  1e-9;
    size_t maxit = 20 * (m + 1) + 50;
    for (size_t it = 0; it < maxit; it++) {
        size_t q = 0;
        for (size_t k = 0; k < m; k++) { mult[k] = 0.0; if (act[k]) A[q++] = k; }
        for (size_t j = 0; j < q; j++)
            fm_chol_solve(L, nv, &Nrm[A[j] * nv], &W[j * nv]);
        for (size_t j = 0; j < q; j++) {
            const double* nj = &Nrm[A[j] * nv];
            double sg = b[A[j]];
            for (size_t t = 0; t < nv; t++) sg += nj[t] * w_g[t];
            rhs[j] = sg;
            for (size_t l = 0; l < q; l++) {
                const double* wl = &W[l * nv];
                double kk = 0.0;
                for (size_t t = 0; t < nv; t++) kk += nj[t] * wl[t];
                K[j * q + l] = kk;
            }
        }
        bool kok = (q == 0);
        if (q > 0) {
            double tau = 0.0;
            for (int attempt = 0; attempt < 6 && !kok; attempt++) {
                for (size_t t = 0; t < q * q; t++) Kc[t] = K[t];
                if (fm_chol_factor(Kc, q, tau)) {
                    fm_chol_solve(Kc, q, rhs, muA); kok = true;
                } else {
                    tau = (tau == 0.0) ? 1e-12 : tau * 100.0;
                }
            }
        }
        if (!kok) {
            /* Dependent active normals: drop the last active inequality and
             * block it from re-entry.  If only equalities remain, they are
             * themselves dependent -> signal the caller to relax. */
            size_t drop = m;
            for (size_t j = q; j-- > 0; ) if (!is_eq[A[j]]) { drop = A[j]; break; }
            if (drop == m) { status = 1; goto done; }
            act[drop] = false; blk[drop] = true;
            continue;
        }
        for (size_t t = 0; t < nv; t++) d[t] = -w_g[t];
        for (size_t j = 0; j < q; j++) {
            mult[A[j]] = muA[j];
            const double* wj = &W[j * nv];
            for (size_t t = 0; t < nv; t++) d[t] += muA[j] * wj[t];
        }
        /* Dual step: drop the most-negative active inequality multiplier. */
        size_t worst = m; double worstv = drop_tol;
        for (size_t j = 0; j < q; j++) {
            size_t k = A[j];
            if (!is_eq[k] && muA[j] < worstv) { worstv = muA[j]; worst = k; }
        }
        if (worst != m) { act[worst] = false; continue; }
        /* Primal step: add the most-violated inactive (non-blocked) inequality. */
        size_t addk = m; double maxviol = feas_tol;
        for (size_t k = 0; k < m; k++) {
            if (act[k] || is_eq[k] || blk[k]) continue;
            const double* nk = &Nrm[k * nv];
            double nd = 0.0;
            for (size_t t = 0; t < nv; t++) nd += nk[t] * d[t];
            double viol = b[k] - nd;                 /* > 0 ⇒ violated */
            if (viol > maxviol) { maxviol = viol; addk = k; }
        }
        if (addk == m) {
            /* KKT for the non-blocked set.  If a BLOCKED inequality is still
             * violated, the linearisation is genuinely inconsistent. */
            status = 0;
            for (size_t k = 0; k < m; k++) {
                if (!blk[k]) continue;
                const double* nk = &Nrm[k * nv];
                double nd = 0.0;
                for (size_t t = 0; t < nv; t++) nd += nk[t] * d[t];
                if (b[k] - nd > 1e-7) { status = 1; break; }
            }
            goto done;
        }
        act[addk] = true;
    }
done:
    free(w_g); free(W); free(act); free(blk); free(A); free(K); free(Kc);
    free(rhs); free(muA);
    return status;
}

/* Assemble and solve the SQP step QP at the current iterate.  L is the
 * Cholesky factor of the Lagrangian-Hessian approximation B (n×n); g=∇f;
 * c[k]=gens[k].expr(x) and Jac (ngens×n, row k = ∇c_k) are the general
 * constraint values / Jacobian; boxes + x give the bound rows.  Fills the
 * step d (n) and the SIGNED general-constraint multipliers lam[k] for the
 * Lagrangian L_ag = f + Σ lam_k c_k.  If the plain QP is inconsistent it
 * re-solves Kraft's relaxed QP (one slack ξ∈[0,1] making d=0 feasible at
 * ξ=1) so a step always exists.  Returns false only on allocation
 * failure. */
bool fm_slsqp_qp(const double* L, size_t n, const double* g,
                        const FmGenCon* gens, const double* c,
                        const double* Jac, size_t ngens,
                        const FmBox* boxes, const double* x,
                        double* d, double* lam) {
    size_t nb = 0;
    for (size_t i = 0; i < n; i++) {
        if (boxes && boxes[i].has_lo) nb++;
        if (boxes && boxes[i].has_hi) nb++;
    }
    size_t m = ngens + nb;
    size_t m1 = m ? m : 1;
    double* Nrm  = (double*)calloc(m1 * n, sizeof(double));
    double* b    = (double*)malloc(sizeof(double) * m1);
    int*    iseq = (int*)malloc(sizeof(int) * m1);
    double* mult = (double*)malloc(sizeof(double) * m1);
    if (!Nrm || !b || !iseq || !mult) {
        free(Nrm); free(b); free(iseq); free(mult); return false;
    }
    /* General constraints: equality a·d = -c (n=a, b=-c); inequality c≤0
     * linearised a·d ≤ -c ⟺ (-a)·d ≥ c (n=-a, b=c). */
    for (size_t k = 0; k < ngens; k++) {
        const double* ak = &Jac[k * n];
        double* row = &Nrm[k * n];
        if (gens[k].equality) {
            for (size_t t = 0; t < n; t++) row[t] = ak[t];
            b[k] = -c[k]; iseq[k] = 1;
        } else {
            for (size_t t = 0; t < n; t++) row[t] = -ak[t];
            b[k] = c[k]; iseq[k] = 0;
        }
    }
    /* Bounds: lo ≤ x+d ⟺ e_i·d ≥ lo-x_i;  x+d ≤ hi ⟺ (-e_i)·d ≥ x_i-hi. */
    size_t idx = ngens;
    for (size_t i = 0; i < n; i++) {
        if (boxes && boxes[i].has_lo) {
            Nrm[idx * n + i] = 1.0;  b[idx] = boxes[i].lo - x[i]; iseq[idx] = 0; idx++;
        }
        if (boxes && boxes[i].has_hi) {
            Nrm[idx * n + i] = -1.0; b[idx] = x[i] - boxes[i].hi; iseq[idx] = 0; idx++;
        }
    }

    int st = fm_slsqp_activeset(L, n, g, Nrm, b, iseq, m, d, mult);

    if (st == 1) {
        /* Inconsistent linearisation -> Kraft slack relaxation on n+1 vars.
         * B' = diag(B, eps_xi) so L' = diag(L, sqrt(eps_xi)); g' = [g; rho];
         * each GENERAL row gets the ξ-column b_k (so at ξ=1 it becomes
         * n·d = 0, met by d=0); bounds are NOT relaxed (d=0 already meets
         * them); ξ∈[0,1] as two extra bound rows. */
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        double rho = 100.0 * (1.0 + gnorm);
        double eps_xi = 1.0;
        size_t nv2 = n + 1, m2 = m + 2;
        double* L2   = (double*)calloc(nv2 * nv2, sizeof(double));
        double* g2   = (double*)malloc(sizeof(double) * nv2);
        double* Nrm2 = (double*)calloc(m2 * nv2, sizeof(double));
        double* b2   = (double*)malloc(sizeof(double) * m2);
        int*    iseq2= (int*)malloc(sizeof(int) * m2);
        double* mul2 = (double*)malloc(sizeof(double) * m2);
        double* d2   = (double*)malloc(sizeof(double) * nv2);
        if (L2 && g2 && Nrm2 && b2 && iseq2 && mul2 && d2) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) L2[i * nv2 + j] = L[i * n + j];
                g2[i] = g[i];
            }
            L2[n * nv2 + n] = sqrt(eps_xi);
            g2[n] = rho;
            for (size_t k = 0; k < m; k++) {
                for (size_t t = 0; t < n; t++) Nrm2[k * nv2 + t] = Nrm[k * n + t];
                Nrm2[k * nv2 + n] = (k < ngens) ? b[k] : 0.0;   /* relax generals only */
                b2[k] = b[k]; iseq2[k] = iseq[k];
            }
            Nrm2[m * nv2 + n] = 1.0;  b2[m] = 0.0;       iseq2[m] = 0;   /* ξ ≥ 0 */
            Nrm2[(m + 1) * nv2 + n] = -1.0; b2[m + 1] = -1.0; iseq2[m + 1] = 0; /* ξ ≤ 1 */
            (void)fm_slsqp_activeset(L2, nv2, g2, Nrm2, b2, iseq2, m2, d2, mul2);
            for (size_t i = 0; i < n; i++) d[i] = d2[i];
            for (size_t k = 0; k < m; k++) mult[k] = mul2[k];
        }
        free(L2); free(g2); free(Nrm2); free(b2); free(iseq2); free(mul2); free(d2);
    }

    /* Map QP multipliers back to the Lagrangian sign convention. */
    for (size_t k = 0; k < ngens; k++)
        lam[k] = gens[k].equality ? -mult[k] : mult[k];

    free(Nrm); free(b); free(iseq); free(mult);
    return true;
}

bool g_fm_slsqp_hit_maxit = false;

bool fm_gens_all_linear(const FmGenCon* gens, size_t ngens, Expr** vars, size_t n) {
    for (size_t k = 0; k < ngens; k++) {
        if (!gens[k].grad_exprs) return false;
        for (size_t i = 0; i < n; i++) {
            const Expr* gi = gens[k].grad_exprs[i];
            if (!gi) return false;
            for (size_t j = 0; j < n; j++)
                if (nm_expr_contains_symbol(gi, vars[j]->data.symbol.name)) return false;
        }
    }
    return true;
}

/* Consecutive non-improving infeasible iterations before SLSQP reports
 * FindMinimum::infeas and stops (see the stall test in the main loop). */
#define FM_SLSQP_STALL 30

bool fm_run_slsqp(Expr* f, Expr** vars, size_t n,
                         FmVarBind* binds, Expr** g_exprs,
                         double* x, /* in/out */
                         const FmGenCon* gens, size_t ngens, double mu,
                         const FmBox* boxes,
                         const FmOpts* opts,
                         double* fx_out) {
    (void)vars; (void)mu;
    size_t ng1 = ngens ? ngens : 1;
    double* B       = (double*)calloc(n * n, sizeof(double));
    double* L       = (double*)malloc(sizeof(double) * n * n);
    double* g       = (double*)malloc(sizeof(double) * n);
    double* g_new   = (double*)malloc(sizeof(double) * n);
    double* d       = (double*)malloc(sizeof(double) * n);
    double* x_new   = (double*)malloc(sizeof(double) * n);
    double* svec    = (double*)malloc(sizeof(double) * n);
    double* yvec    = (double*)malloc(sizeof(double) * n);
    double* Bs      = (double*)malloc(sizeof(double) * n);
    double* xbest   = (double*)malloc(sizeof(double) * n);
    double* lam     = (double*)calloc(ng1, sizeof(double));
    double* pen     = (double*)calloc(ng1, sizeof(double));
    double* c       = (double*)malloc(sizeof(double) * ng1);
    double* c_new   = (double*)malloc(sizeof(double) * ng1);
    double* Jac     = (double*)malloc(sizeof(double) * ng1 * n);
    double* Jac_new = (double*)malloc(sizeof(double) * ng1 * n);
    bool ok = false;
    double fx = 0.0;
    bool have_best = false; double best_f = 0.0;
    g_fm_slsqp_hit_maxit = false;
    if (!B || !L || !g || !g_new || !d || !x_new || !svec || !yvec || !Bs
        || !xbest || !lam || !pen || !c || !c_new || !Jac || !Jac_new)
        goto cleanup;

    for (size_t i = 0; i < n; i++) B[i * n + i] = 1.0;
    if (boxes) fm_project_box(x, n, boxes);

    if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) {
        fm_warn(g_fm_name, "nlnum", "objective evaluation failed at start point");
        goto cleanup;
    }
    bool central = false, exact_grad = false;
    if (!fm_slsqp_objgrad(f, g_exprs, binds, x, n, opts, g, fx, boxes,
                          central, &exact_grad)) {
        fm_warn(g_fm_name, "nlnum", "gradient evaluation failed at start point");
        goto cleanup;
    }
    for (size_t k = 0; k < ngens; k++) {
        if (!fm_eval_scalar(gens[k].expr, binds, x, n, opts, &c[k])
            || !fm_slsqp_congrad(&gens[k], binds, x, n, opts, &Jac[k * n])) {
            fm_warn(g_fm_name, "nlnum", "constraint evaluation failed at start point");
            goto cleanup;
        }
    }
    /* Objective scaling.  B starts at the identity and the stationarity test
     * is absolute, both of which assume an objective whose gradient is O(1).
     * A small-scale objective (a compliance of 6e-4 over 64 cell values of
     * size 50: |grad f| ~ 1e-6) made every QP step a 1e-6 nudge, the merit
     * descent test fail by rounding, and 500 iterations end where they began.
     * SLSQP therefore minimises fs f with fs = 1/|grad f(x0)|_inf (capped at
     * 1e8), the scaling fm_run_penalty applies, and reports f unscaled. */
    double fs = 1.0;
    {
        double gmax = 0.0;
        for (size_t i = 0; i < n; i++) if (fabs(g[i]) > gmax) gmax = fabs(g[i]);
        if (gmax > 0.0 && gmax < 1e-2 && isfinite(gmax))
            fs = (1.0 / gmax < 1e8) ? 1.0 / gmax : 1e8;
    }
    fx *= fs;
    for (size_t i = 0; i < n; i++) g[i] *= fs;

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    int    infeas_streak = 0;
    double best_viol = HUGE_VAL;       /* smallest max-violation seen so far */
    int    zero_streak = 0;

    int64_t it;
    for (it = 0; it < opts->max_iter; it++) {
        /* Factor a COPY of B (tau-retry); reset B=I if it is not SPD. */
        bool fac = false; double tau = 0.0;
        for (int attempt = 0; attempt < 6 && !fac; attempt++) {
            for (size_t t = 0; t < n * n; t++) L[t] = B[t];
            if (fm_chol_factor(L, n, tau)) fac = true;
            else tau = (tau == 0.0) ? 1e-10 : tau * 100.0;
        }
        if (!fac) {
            for (size_t t = 0; t < n * n; t++) B[t] = 0.0;
            for (size_t i = 0; i < n; i++) B[i * n + i] = 1.0;
            for (size_t t = 0; t < n * n; t++) L[t] = B[t];
            fm_chol_factor(L, n, 0.0);
        }

        if (!fm_slsqp_qp(L, n, g, gens, c, Jac, ngens, boxes, x, d, lam))
            goto cleanup;

        double dnorm = 0.0;
        for (size_t i = 0; i < n; i++) if (fabs(d[i]) > dnorm) dnorm = fabs(d[i]);

        /* Powell penalty update: pen_k = max(|lam_k|, (pen_k+|lam_k|)/2), so
         * pen_k ≥ |multiplier| ⇒ d is a merit descent direction. */
        for (size_t k = 0; k < ngens; k++) {
            double a = fabs(lam[k]);
            double avg = 0.5 * (pen[k] + a);
            double np = (a > avg) ? a : avg;
            if (np < a) np = a;
            if (np > 1e12) np = 1e12;
            pen[k] = np;
        }

        /* L1 merit at x and its directional derivative along d. */
        double phi0 = fx, dphi = 0.0;
        for (size_t k = 0; k < ngens; k++) {
            double viol = gens[k].equality ? fabs(c[k]) : (c[k] > 0.0 ? c[k] : 0.0);
            phi0 += pen[k] * viol;
            dphi -= pen[k] * viol;
        }
        for (size_t i = 0; i < n; i++) dphi += g[i] * d[i];
        if (dphi >= 0.0) {
            /* Merit ascent (penalties still too small / tiny step): fall back
             * to the guaranteed-descent estimate -dᵀBd. */
            double dBd = 0.0;
            for (size_t i = 0; i < n; i++) {
                double t = 0.0;
                for (size_t j = 0; j < n; j++) t += B[i * n + j] * d[j];
                dBd += d[i] * t;
            }
            dphi = -dBd;
        }

        /* Armijo backtracking on the merit; track the best trial step. */
        double alpha = 1.0, best_alpha = 0.0, best_phi = phi0, best_ft = fx;
        bool ls_ok = false;
        for (int bt = 0; bt < 40; bt++) {
            double ft;
            for (size_t i = 0; i < n; i++) x_new[i] = x[i] + alpha * d[i];
            bool eok = fm_eval_scalar(f, binds, x_new, n, opts, &ft);
            ft *= fs;
            bool cok = eok && isfinite(ft);
            double phit = ft;
            for (size_t k = 0; cok && k < ngens; k++) {
                double cv;
                if (!fm_eval_scalar(gens[k].expr, binds, x_new, n, opts, &cv)) { cok = false; break; }
                phit += pen[k] * (gens[k].equality ? fabs(cv) : (cv > 0.0 ? cv : 0.0));
            }
            if (cok) {
                if (phit <= phi0 + 1e-4 * alpha * dphi) {
                    best_alpha = alpha; best_phi = phit; best_ft = ft; ls_ok = true; break;
                }
                if (phit < best_phi) { best_phi = phit; best_alpha = alpha; best_ft = ft; }
            }
            alpha *= 0.5;
            if (alpha < 1e-12) break;
        }

        if (!ls_ok && best_alpha == 0.0) {
            /* No trial improved the merit -> B is a poor model: reset and
             * retry once; two in a row means we are stationary. */
            zero_streak++;
            for (size_t t = 0; t < n * n; t++) B[t] = 0.0;
            for (size_t i = 0; i < n; i++) B[i * n + i] = 1.0;
            if (zero_streak >= 2) { ok = true; break; }
            continue;
        }
        zero_streak = 0;
        alpha = best_alpha;
        for (size_t i = 0; i < n; i++) x_new[i] = x[i] + alpha * d[i];
        double fx_new = best_ft;     /* the line search's (scaled) value there */
        fm_fire_monitor(opts->step_monitor);

        /* Gradient + Jacobian at x_new (needed for the Lagrangian BFGS pair). */
        /* switch the differences to central once the steps are small */
        if (!central && !exact_grad
            && alpha * dnorm < 1e-4 * (fm_slsqp_infnorm(x, n) + 1.0))
            central = true;
        bool grads_ok = fm_slsqp_objgrad(f, g_exprs, binds, x_new, n, opts, g_new,
                                         fx_new / fs, boxes, central, &exact_grad);
        for (size_t i = 0; grads_ok && i < n; i++) g_new[i] *= fs;
        for (size_t k = 0; grads_ok && k < ngens; k++) {
            if (!fm_eval_scalar(gens[k].expr, binds, x_new, n, opts, &c_new[k])
                || !fm_slsqp_congrad(&gens[k], binds, x_new, n, opts, &Jac_new[k * n]))
                grads_ok = false;
        }
        if (!grads_ok) {
            /* Take the step and stop -- cannot form the next model. */
            for (size_t i = 0; i < n; i++) x[i] = x_new[i];
            fx = fx_new;
            ok = true; break;
        }

        /* Powell-damped BFGS on B using y = ∇L(x_new) − ∇L(x) at lam. */
        for (size_t i = 0; i < n; i++) {
            double gl_new = g_new[i], gl_old = g[i];
            for (size_t k = 0; k < ngens; k++) {
                gl_new += lam[k] * Jac_new[k * n + i];
                gl_old += lam[k] * Jac[k * n + i];
            }
            yvec[i] = gl_new - gl_old;
            svec[i] = x_new[i] - x[i];
        }
        for (size_t i = 0; i < n; i++) {
            double t = 0.0;
            for (size_t j = 0; j < n; j++) t += B[i * n + j] * svec[j];
            Bs[i] = t;
        }
        double sBs = 0.0, sy = 0.0;
        for (size_t i = 0; i < n; i++) { sBs += svec[i] * Bs[i]; sy += svec[i] * yvec[i]; }
        if (sBs > 1e-12) {
            double theta = (sy >= 0.2 * sBs) ? 1.0 : (0.8 * sBs) / (sBs - sy);
            double sr = 0.0;
            for (size_t i = 0; i < n; i++) {
                yvec[i] = theta * yvec[i] + (1.0 - theta) * Bs[i];   /* r */
                sr += svec[i] * yvec[i];
            }
            if (sr > 1e-12) {
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < n; j++)
                        B[i * n + j] += yvec[i] * yvec[j] / sr - Bs[i] * Bs[j] / sBs;
            }
        }

        /* Commit the step. */
        for (size_t i = 0; i < n; i++) { x[i] = x_new[i]; g[i] = g_new[i]; }
        fx = fx_new;
        for (size_t k = 0; k < ngens; k++) {
            c[k] = c_new[k];
            for (size_t t = 0; t < n; t++) Jac[k * n + t] = Jac_new[k * n + t];
        }

        /* Feasibility, best-feasible tracking, and convergence. */
        double viol = 0.0;
        for (size_t k = 0; k < ngens; k++) {
            double v = gens[k].equality ? fabs(c[k]) : (c[k] > 0.0 ? c[k] : 0.0);
            if (v > viol) viol = v;
        }
        if (viol <= 1e-8 && (!have_best || fx < best_f)) {
            for (size_t i = 0; i < n; i++) xbest[i] = x[i];
            best_f = fx; have_best = true;
        }
        double gLnorm = 0.0, xnorm = 0.0;
        for (size_t i = 0; i < n; i++) {
            double gl = g[i];
            for (size_t k = 0; k < ngens; k++) gl += lam[k] * Jac[k * n + i];
            if (fabs(gl) > gLnorm) gLnorm = fabs(gl);
            if (fabs(x[i]) > xnorm) xnorm = fabs(x[i]);
        }
        bool feas = (viol <= tol_acc);
        bool small_step = (alpha * dnorm < tol_prec * (xnorm + 1.0));
        if (feas && (small_step || gLnorm < tol_acc)) { ok = true; break; }

        /* Give up on feasibility only when the violation has STALLED: no 1%
         * improvement on the best violation seen for FM_SLSQP_STALL straight
         * iterations. The old rule quit after any 8 consecutive infeasible
         * iterations, even while the violation was still falling -- so an
         * infeasible start on a problem with many coupled equalities (a
         * 10-link hanging chain) was abandoned mid-approach at a point
         * violating its constraints by 2e-2. */
        if (viol > 1e-6) {
            if (viol < 0.99 * best_viol) { best_viol = viol; infeas_streak = 0; }
            else if (++infeas_streak >= FM_SLSQP_STALL) {
                fm_warn(g_fm_name, "infeas",
                        "could not satisfy constraints to tolerance");
                break;
            }
        } else { infeas_streak = 0; if (viol < best_viol) best_viol = viol; }
    }
    /* Loop exhausted without meeting the goals: say so (cvmit). */
    g_fm_slsqp_hit_maxit = (it >= opts->max_iter);
    if (it >= opts->max_iter) fm_warn_maxit(opts);

    /* If we ended infeasible but saw a feasible iterate, return the best. */
    {
        double final_viol = 0.0;
        for (size_t k = 0; k < ngens; k++) {
            double v = gens[k].equality ? fabs(c[k]) : (c[k] > 0.0 ? c[k] : 0.0);
            if (v > final_viol) final_viol = v;
        }
        if (final_viol > 1e-8 && have_best) {
            for (size_t i = 0; i < n; i++) x[i] = xbest[i];
            fx = best_f;
        }
    }
    *fx_out = fx / fs;
    ok = true;
cleanup:
    free(B); free(L); free(g); free(g_new); free(d); free(x_new);
    free(svec); free(yvec); free(Bs); free(xbest); free(lam); free(pen);
    free(c); free(c_new); free(Jac); free(Jac_new);
    return ok;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Linear constraints at WorkingPrecision                              *
 * ------------------------------------------------------------------ *
 * The machine SQP above locates the solution and, with it, the ACTIVE SET:
 * which box bounds hold with equality and which linear inequalities are
 * tight.  On that set the problem is an equality-constrained one, and linear
 * equalities can be eliminated exactly: each independent active row solves
 * (by Solve, over its exact coefficients) for one pivot variable, the bound
 * variables are pinned to their bounds, and what is left is an UNCONSTRAINED
 * problem in the remaining variables -- the objective with the pivots
 * substituted -- which the MPFR quasi-Newton (fm_run_bfgs_mpfr) solves to the
 * working precision.  The pivots are then recovered from their exact
 * expressions.  The refined point is accepted only if it keeps every pivot
 * inside its box and every inactive inequality satisfied; a wrong active set
 * fails that test and the caller falls back to the machine solution with the
 * FindMinimum::nimpl warning.  A pinned bound is the double the box parser
 * stored, taken exactly (an integer, else that double in MPFR): exact for the
 * integer and dyadic bounds such problems carry. */

static bool fm_lin_num(const Expr* e, double* out) {
    Expr* ev = eval_and_free(expr_copy((Expr*)e));
    Expr* nv = ev ? numericalize(ev, numeric_machine_spec()) : NULL;
    bool ok = nv && fm_expr_to_double_real(nv, out);
    expr_free(ev); expr_free(nv);
    return ok;
}

/* lhs -> rhs (lhs copied, rhs adopted) */
static Expr* fm_lin_rule(Expr* lhs, Expr* rhs) {
    Expr* a[2] = { expr_copy(lhs), rhs };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* e /. rules, evaluated (both borrowed) */
static Expr* fm_lin_subst(Expr* e, Expr* rules) {
    Expr* a[2] = { expr_copy(e), expr_copy(rules) };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll), a, 2));
}

bool fm_linear_refine_mpfr(Expr* f, Expr** vars, size_t n, FmVarBind* binds,
                           const FmGenCon* gens, size_t ngens,
                           const FmBox* boxes, const FmOpts* opts,
                           const double* xm, mpfr_t* x_out, mpfr_t fx_out) {
    long bits = opts->wp_bits;
    bool ok = false;
    for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    size_t n1 = n ? n : 1, g1 = ngens ? ngens : 1;
    int*    role = (int*)calloc(n1, sizeof(int));     /* 0 free, 1 lo, 2 hi, 3 pivot */
    double* A    = (double*)calloc(g1 * n1, sizeof(double));
    int*    act  = (int*)calloc(g1, sizeof(int));
    Expr**  sub  = (Expr**)malloc(sizeof(Expr*) * n1);
    FmVarBind* sb = (FmVarBind*)malloc(sizeof(FmVarBind) * n1);
    FmBox*  sbox = (FmBox*)malloc(sizeof(FmBox) * n1);
    size_t* sidx = (size_t*)malloc(sizeof(size_t) * n1);
    Expr *fixed = NULL, *sol = NULL, *allr = NULL, *fred = NULL;
    Expr** gm = NULL;
    mpfr_t* xs = NULL;
    size_t m = 0, np = 0;
    if (!role || !A || !act || !sub || !sb || !sbox || !sidx) goto done;

    /* 1. bounds active at the machine solution */
    for (size_t i = 0; i < n; i++) {
        if (boxes[i].has_lo && fabs(xm[i] - boxes[i].lo) <= 1e-7 * (1.0 + fabs(boxes[i].lo)))
            role[i] = 1;
        else if (boxes[i].has_hi && fabs(xm[i] - boxes[i].hi) <= 1e-7 * (1.0 + fabs(boxes[i].hi)))
            role[i] = 2;
    }
    /* 2. active rows (every equality, the tight inequalities) and their
     *    machine coefficients */
    for (size_t k = 0; k < ngens; k++) {
        double scale = 0.0;
        for (size_t i = 0; i < n; i++) {
            if (!fm_lin_num(gens[k].grad_exprs[i], &A[k * n + i])) goto done;
            scale += fabs(A[k * n + i] * xm[i]);
        }
        double c;
        if (!fm_eval_scalar(gens[k].expr, binds, xm, n, opts, &c)) goto done;
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
        if (gens[k].equality || fabs(c) <= 1e-7 * (1.0 + scale)) act[k] = 1;
    }
    /* 3. one pivot per independent active row, by elimination with the
     *    largest free coefficient; a dependent row (pivot ~ 0) is dropped */
    {
        double* R = (double*)malloc(sizeof(double) * g1 * n1);
        if (!R) goto done;
        for (size_t t = 0; t < ngens * n; t++) R[t] = A[t];
        double amax = 0.0;
        for (size_t t = 0; t < ngens * n; t++) if (fabs(R[t]) > amax) amax = fabs(R[t]);
        for (size_t k = 0; k < ngens; k++) {
            if (!act[k]) continue;
            size_t pc = n; double pv = 0.0;
            for (size_t i = 0; i < n; i++)
                if (role[i] == 0 && fabs(R[k * n + i]) > pv) { pv = fabs(R[k * n + i]); pc = i; }
            if (pc == n || pv <= 1e-12 * amax) { act[k] = 0; continue; }
            role[pc] = 3;
            for (size_t l = k + 1; l < ngens; l++) {
                if (!act[l] || R[l * n + pc] == 0.0) continue;
                double fct = R[l * n + pc] / R[k * n + pc];
                for (size_t i = 0; i < n; i++) R[l * n + i] -= fct * R[k * n + i];
            }
        }
        free(R);
    }
    /* 4. exact rules: the pinned bounds, then the pivots by Solve */
    {
        size_t nf = 0;
        Expr** fr = (Expr**)malloc(sizeof(Expr*) * n1);
        for (size_t i = 0; i < n; i++)
            if (role[i] == 1 || role[i] == 2) {
                double d = role[i] == 1 ? boxes[i].lo : boxes[i].hi;
                Expr* bx;
                if (d == floor(d) && fabs(d) < 9.0e15) {
                    bx = expr_new_integer((int64_t)d);
                } else {                     /* the double, exactly, in MPFR */
                    mpfr_t b; mpfr_init2(b, bits > 64 ? bits : 64);
                    mpfr_set_d(b, d, MPFR_RNDN);
                    bx = expr_new_mpfr_copy(b);
                    mpfr_clear(b);
                }
                fr[nf++] = fm_lin_rule(vars[i], bx);
            }
        fixed = expr_new_function(expr_new_symbol(SYM_List), fr, nf);
        free(fr);
    }
    {
        size_t ne = 0;
        Expr** eqs = (Expr**)malloc(sizeof(Expr*) * g1);
        Expr** pv  = (Expr**)malloc(sizeof(Expr*) * n1);
        for (size_t k = 0; k < ngens; k++) if (act[k]) {
            Expr* ea[2] = { fm_lin_subst(gens[k].expr, fixed), expr_new_integer(0) };
            eqs[ne++] = expr_new_function(expr_new_symbol(SYM_Equal), ea, 2);
        }
        for (size_t i = 0; i < n; i++) if (role[i] == 3) pv[np++] = expr_copy(vars[i]);
        if (ne > 0) {
            Expr* sa[2] = { expr_new_function(expr_new_symbol(SYM_List), eqs, ne),
                            expr_new_function(expr_new_symbol(SYM_List), pv, np) };
            sol = eval_and_free(expr_new_function(expr_new_symbol(SYM_Solve), sa, 2));
        } else {
            for (size_t i = 0; i < np; i++) expr_free(pv[i]);
        }
        free(eqs); free(pv);
        /* {{p1 -> e1, ...}}: one solution, one rule per pivot */
        if (ne > 0) {
            if (!sol || sol->type != EXPR_FUNCTION || sol->data.function.arg_count != 1)
                goto done;
            Expr* s0 = sol->data.function.args[0];
            if (s0->type != EXPR_FUNCTION || s0->data.function.arg_count != np) goto done;
        }
        size_t nfx = fixed->data.function.arg_count;
        size_t nr = nfx + (sol ? np : 0);
        Expr** ar = (Expr**)malloc(sizeof(Expr*) * (nr ? nr : 1));
        size_t q = 0;
        for (size_t i = 0; i < nfx; i++) ar[q++] = expr_copy(fixed->data.function.args[i]);
        if (sol) for (size_t i = 0; i < np; i++)
            ar[q++] = expr_copy(sol->data.function.args[0]->data.function.args[i]);
        allr = expr_new_function(expr_new_symbol(SYM_List), ar, q);
        free(ar);
    }
    /* 5. the reduced, unconstrained problem in the remaining variables */
    fred = fm_lin_subst(f, allr);
    for (size_t i = 0; i < n; i++) if (role[i] == 0) {
        sub[m] = vars[i]; sb[m] = binds[i]; sidx[m] = i;
        sbox[m].has_lo = sbox[m].has_hi = false;   /* interior at the solution */
        m++;
    }
    xs = fm_mpfr_array(m ? m : 1, bits);
    for (size_t j = 0; j < m; j++) mpfr_set_d(xs[j], xm[sidx[j]], MPFR_RNDN);
    if (m > 0) {
        gm = fm_compute_gradient(fred, sub, m);
        if (!fm_run_bfgs_mpfr(fred, sub, m, sb, gm, xs, sbox, opts, fx_out)) goto done;
        for (size_t j = 0; j < m; j++) fm_bind_clear_temp(&sb[j]);
    }
    /* 6. every coordinate at the working precision, the pivots in their box */
    {
        Expr** vr = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t j = 0; j < m; j++) vr[j] = fm_lin_rule(sub[j], expr_new_mpfr_copy(xs[j]));
        Expr* vals = expr_new_function(expr_new_symbol(SYM_List), vr, m);
        free(vr);
        NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = bits;
        bool good = true;
        for (size_t i = 0; i < n && good; i++) {
            Expr* xi = fm_lin_subst(vars[i], allr);
            Expr* xv = xi ? fm_lin_subst(xi, vals) : NULL;
            Expr* xn = xv ? numericalize(xv, sp) : NULL;
            mpfr_t im; mpfr_init2(im, bits);
            good = xn && get_approx_mpfr(xn, x_out[i], im, NULL) && mpfr_zero_p(im);
            mpfr_clear(im);
            expr_free(xi); expr_free(xv); expr_free(xn);
            if (good && role[i] == 3) {
                double d = mpfr_get_d(x_out[i], MPFR_RNDN);
                double tol = 1e-12 * (1.0 + fabs(d));
                if ((boxes[i].has_lo && d < boxes[i].lo - tol)
                    || (boxes[i].has_hi && d > boxes[i].hi + tol)) good = false;
            }
        }
        if (m == 0 && good) {
            Expr* fv = fm_lin_subst(fred, vals);
            Expr* fn = fv ? numericalize(fv, sp) : NULL;
            mpfr_t im; mpfr_init2(im, bits);
            good = fn && get_approx_mpfr(fn, fx_out, im, NULL);
            mpfr_clear(im);
            expr_free(fv); expr_free(fn);
        }
        expr_free(vals);
        if (!good) goto done;
    }
    /* 7. the inactive inequalities still hold */
    {
        double* xd = (double*)malloc(sizeof(double) * n1);
        bool good = xd != NULL;
        for (size_t i = 0; good && i < n; i++) xd[i] = mpfr_get_d(x_out[i], MPFR_RNDN);
        for (size_t k = 0; k < ngens && good; k++) {
            if (gens[k].equality) continue;
            double c;
            if (!fm_eval_scalar(gens[k].expr, binds, xd, n, opts, &c) || c > 1e-9) good = false;
        }
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
        free(xd);
        if (!good) goto done;
    }
    ok = true;
done:
    for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    if (gm) { for (size_t j = 0; j < m; j++) expr_free(gm[j]); free(gm); }
    if (xs) fm_mpfr_array_free(xs, m ? m : 1);
    expr_free(fixed); expr_free(sol); expr_free(allr); expr_free(fred);
    free(role); free(A); free(act); free(sub); free(sb); free(sbox); free(sidx);
    return ok;
}
#endif /* USE_MPFR */
