/* findmin_trust.c — trust-region family (dogleg / trust-ncg / trust-exact / trust-krylov / Newton-CG).
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* Solve ||z + τ d||² = Δ² for the two real roots τlo ≤ τhi (d ≠ 0). The
 * discriminant is non-negative whenever ||z|| ≤ Δ; clamped to 0 for safety. */
static void fm_tr_boundary(const double* z, const double* d, size_t n,
                           double Delta, double* tlo, double* thi) {
    double a = fm_dot(d, d, n);
    double b = 2.0 * fm_dot(z, d, n);
    double c = fm_dot(z, z, n) - Delta * Delta;
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) disc = 0.0;
    double s = sqrt(disc);
    if (a <= 0.0) { *tlo = 0.0; *thi = 0.0; return; }
    double t1 = (-b - s) / (2.0 * a);
    double t2 = (-b + s) / (2.0 * a);
    *tlo = t1; *thi = t2;
}

/* Build the dense Hessian B at x: exact symbolic (fm_eval_hessian) when the
 * n×n derivative expressions exist, else column by column via n Hessian-vector
 * products on the unit vectors (symmetrised). g_base = ∇f(x) is cached. */
static bool fm_tr_build_hessian(const FmLbfgsCtx* c, Expr*** H_exprs,
                                const double* x, const double* g_base, double* B,
                                double* ej, double* col,
                                double* xpert, double* gpert, const bool* active0) {
    size_t n = c->n;
    if (H_exprs && fm_eval_hessian(H_exprs, c->binds, x, n, c->opts, B))
        return true;
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) ej[i] = (i == j) ? 1.0 : 0.0;
        if (!fm_tnc_hessvec(c, x, g_base, ej, active0, xpert, gpert, col))
            return false;
        for (size_t i = 0; i < n; i++) B[i * n + j] = col[i];
    }
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++) {
            double a = 0.5 * (B[i * n + j] + B[j * n + i]);
            B[i * n + j] = a; B[j * n + i] = a;
        }
    return true;
}

/* ---- dogleg (needs dense SPD B) ---------------------------------------- */
bool fm_tr_dogleg(const FmQuad* q, const double* g, double gnorm,
                         double Delta, double* p, bool* hits) {
    size_t n = q->n;
    const double* B = q->B;
    double* L    = (double*)malloc(n * n * sizeof(double));
    double* pB   = (double*)malloc(n * sizeof(double));
    double* pU   = (double*)malloc(n * sizeof(double));
    double* tmp  = (double*)malloc(n * sizeof(double));
    double* negg = (double*)malloc(n * sizeof(double));
    bool ret = false;
    if (!L || !pB || !pU || !tmp || !negg) goto done;

    for (size_t i = 0; i < n; i++) negg[i] = -g[i];

    /* Newton point pB = -B⁻¹g via Cholesky (τ = 0: keep the true B). */
    bool have_pB = false;
    memcpy(L, B, n * n * sizeof(double));
    if (fm_chol_factor(L, n, 0.0)) {
        fm_chol_solve(L, n, negg, pB);
        have_pB = true;
        if (sqrt(fm_dot(pB, pB, n)) <= Delta) {         /* full Newton, interior */
            for (size_t i = 0; i < n; i++) p[i] = pB[i];
            *hits = false; ret = true; goto done;
        }
    }

    /* Cauchy point pU = -(gᵀg / gᵀBg) g. */
    fm_matvec(B, n, g, tmp);
    double gBg = fm_dot(g, tmp, n);
    if (!(gBg > 0.0) || gnorm == 0.0) {                 /* no positive curvature */
        double sc = (gnorm > 0.0) ? (Delta / gnorm) : 0.0;
        for (size_t i = 0; i < n; i++) p[i] = -sc * g[i];
        *hits = (gnorm > 0.0); ret = true; goto done;
    }
    double tauU = (gnorm * gnorm) / gBg;
    for (size_t i = 0; i < n; i++) pU[i] = -tauU * g[i];
    double pUn = sqrt(fm_dot(pU, pU, n));
    if (pUn >= Delta || !have_pB) {                     /* Cauchy past boundary  */
        double sc = (pUn > 0.0) ? (Delta / pUn) : 0.0;
        for (size_t i = 0; i < n; i++) p[i] = sc * pU[i];
        *hits = true; ret = true; goto done;
    }

    /* Dogleg leg pU → pB: p = pU + τ(pB - pU), ||p|| = Δ, τ ∈ [0, 1]. */
    for (size_t i = 0; i < n; i++) tmp[i] = pB[i] - pU[i];
    double tlo, thi; fm_tr_boundary(pU, tmp, n, Delta, &tlo, &thi);
    double tau = thi;
    if (tau < 0.0) tau = 0.0; else if (tau > 1.0) tau = 1.0;
    for (size_t i = 0; i < n; i++) p[i] = pU[i] + tau * tmp[i];
    *hits = true; ret = true;
done:
    free(L); free(pB); free(pU); free(tmp); free(negg);
    return ret;
}

/* ---- trust-ncg: Steihaug-Toint truncated CG (Hessian-vector only) ------ */
bool fm_tr_steihaug(const FmQuad* q, const double* g, double gnorm,
                           double Delta, double* p, bool* hits) {
    size_t n = q->n;
    double* z  = (double*)calloc(n, sizeof(double));
    double* r  = (double*)malloc(n * sizeof(double));
    double* d  = (double*)malloc(n * sizeof(double));
    double* Bd = (double*)malloc(n * sizeof(double));
    bool ret = false;
    if (!z || !r || !d || !Bd) goto done;

    *hits = false;
    if (gnorm == 0.0) { for (size_t i = 0; i < n; i++) p[i] = 0.0; ret = true; goto done; }
    for (size_t i = 0; i < n; i++) { r[i] = g[i]; d[i] = -g[i]; }
    double eta = sqrt(gnorm); if (eta > 0.5) eta = 0.5;
    double tol = eta * gnorm;
    double rr = fm_dot(r, r, n);
    size_t maxit = 2 * n + 1; if (maxit > 1000) maxit = 1000;

    for (size_t j = 0; j < maxit; j++) {
        if (!fm_quad_matvec(q, d, Bd)) {                /* Hv failed             */
            if (j == 0) {
                double sc = Delta / gnorm;
                for (size_t i = 0; i < n; i++) p[i] = -sc * g[i];
                *hits = true;
            } else {
                for (size_t i = 0; i < n; i++) p[i] = z[i];
            }
            ret = true; goto done;
        }
        double dBd = fm_dot(d, Bd, n), dd = fm_dot(d, d, n);
        if (dBd <= 1e-16 * dd) {                        /* negative curvature    */
            double tlo, thi; fm_tr_boundary(z, d, n, Delta, &tlo, &thi);
            /* pick the boundary root with the lower model, using Bz = r - g. */
            double gz = fm_dot(g, z, n), gd = fm_dot(g, d, n);
            double zBz = 0.0; for (size_t i = 0; i < n; i++) zBz += z[i] * (r[i] - g[i]);
            double zBd = fm_dot(z, Bd, n);
            double mlo = gz + tlo * gd + 0.5 * (zBz + 2.0 * tlo * zBd + tlo * tlo * dBd);
            double mhi = gz + thi * gd + 0.5 * (zBz + 2.0 * thi * zBd + thi * thi * dBd);
            double tau = (mlo < mhi) ? tlo : thi;
            for (size_t i = 0; i < n; i++) p[i] = z[i] + tau * d[i];
            *hits = true; ret = true; goto done;
        }
        double alpha = rr / dBd;
        double zn2 = 0.0;
        for (size_t i = 0; i < n; i++) { double zi = z[i] + alpha * d[i]; zn2 += zi * zi; }
        if (zn2 >= Delta * Delta) {                     /* step leaves the ball  */
            double tlo, thi; fm_tr_boundary(z, d, n, Delta, &tlo, &thi);
            for (size_t i = 0; i < n; i++) p[i] = z[i] + thi * d[i];
            *hits = true; ret = true; goto done;
        }
        for (size_t i = 0; i < n; i++) { z[i] += alpha * d[i]; r[i] += alpha * Bd[i]; }
        double rr_new = fm_dot(r, r, n);
        if (sqrt(rr_new) < tol) {                       /* interior solution     */
            for (size_t i = 0; i < n; i++) p[i] = z[i];
            ret = true; goto done;
        }
        double beta = rr_new / rr;
        for (size_t i = 0; i < n; i++) d[i] = -r[i] + beta * d[i];
        rr = rr_new;
    }
    for (size_t i = 0; i < n; i++) p[i] = z[i];          /* maxit: best interior  */
    ret = true;
done:
    free(z); free(r); free(d); free(Bd);
    return ret;
}

/* ---- trust-exact: Moré-Sorensen near-exact subproblem (dense B) --------- *
 * Find λ ≥ 0 and p solving (B + λI)p = -g with B + λI ⪰ 0 and ||p|| ≈ Δ.
 * φ(λ) = 1/||p(λ)|| - 1/Δ is smooth and near-linear, so a Newton iteration
 * λ ← λ + (||p||/||w||)²·(||p||-Δ)/Δ  (L w = p, forward solve) converges in a
 * handful of steps, safeguarded by a bracket [lo, hi] and PD enforcement via
 * the Cholesky retry. The hard case (g ⟂ the least-eigenvector, ||p|| < Δ at
 * the PD threshold) is handled by adding the least-eigenvector direction
 * (inverse iteration) out to the boundary. */
bool fm_tr_moresorensen(const FmQuad* q, const double* g, double gnorm,
                               double Delta, double* p, bool* hits) {
    (void)gnorm;                 /* bracket is discovered dynamically, not from ||g|| */
    size_t n = q->n;
    const double* B = q->B;
    double* L    = (double*)malloc(n * n * sizeof(double));
    double* negg = (double*)malloc(n * sizeof(double));
    double* wv   = (double*)malloc(n * sizeof(double));
    double* zev  = (double*)malloc(n * sizeof(double));
    double* Bp   = (double*)malloc(n * sizeof(double));
    bool ret = false;
    if (!L || !negg || !wv || !zev || !Bp) goto done;
    for (size_t i = 0; i < n; i++) negg[i] = -g[i];

    /* Gershgorin lower bound and inf-norm for the initial λ bracket. */
    double gersh_lb = HUGE_VAL, hinf = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diag = B[i * n + i], rad = 0.0, rowsum = 0.0;
        for (size_t j = 0; j < n; j++) { double a = fabs(B[i * n + j]); rowsum += a; if (j != i) rad += a; }
        if (diag - rad < gersh_lb) gersh_lb = diag - rad;
        if (rowsum > hinf) hinf = rowsum;
    }
    double lam_lo = 0.0, lam_hi = HUGE_VAL;
    double lam = (gersh_lb < 0.0) ? -gersh_lb : 0.0;    /* PD-safe starting shift */
    const double k_easy = 0.1;
    bool have_p = false; double pnorm = 0.0;

    for (int it = 0; it < 60; it++) {
        memcpy(L, B, n * n * sizeof(double));
        if (!fm_chol_factor(L, n, lam)) {               /* not PD: raise lower bd */
            if (lam > lam_lo) lam_lo = lam;
            double nl = (lam_hi < HUGE_VAL) ? 0.5 * (lam_lo + lam_hi)
                                            : (lam > 0.0 ? lam * 2.0 : 1.0);
            if (nl <= lam_lo) nl = lam_lo + (lam_lo > 0.0 ? 0.1 * lam_lo : 0.1);
            lam = nl; continue;
        }
        fm_chol_solve(L, n, negg, p);
        pnorm = sqrt(fm_dot(p, p, n)); have_p = true;
        if (lam == 0.0 && pnorm <= Delta) {             /* B PD, Newton interior */
            *hits = false; ret = true; goto done;
        }
        if (fabs(pnorm - Delta) <= k_easy * Delta) {    /* converged (easy case) */
            *hits = true; ret = true; goto done;
        }
        /* forward solve L w = p, then Newton update on 1/||p||. */
        for (size_t i = 0; i < n; i++) {
            double s = p[i];
            for (size_t kk = 0; kk < i; kk++) s -= L[i * n + kk] * wv[kk];
            wv[i] = s / L[i * n + i];
        }
        double wn2 = fm_dot(wv, wv, n);
        if (pnorm < Delta) { if (lam < lam_hi) lam_hi = lam; }
        else               { if (lam > lam_lo) lam_lo = lam; }
        /* hard case: bracket collapsed with ||p|| < Δ. */
        if (lam_hi < HUGE_VAL && (lam_hi - lam_lo) <= 1e-10 * (1.0 + lam_hi) && pnorm < Delta) {
            for (size_t i = 0; i < n; i++) zev[i] = 1.0 / sqrt((double)n);
            for (int ii = 0; ii < 3; ii++) {            /* inverse iteration     */
                fm_chol_solve(L, n, zev, wv);
                double wn = sqrt(fm_dot(wv, wv, n));
                if (wn <= 0.0) break;
                for (size_t i = 0; i < n; i++) zev[i] = wv[i] / wn;
            }
            double tlo, thi; fm_tr_boundary(p, zev, n, Delta, &tlo, &thi);
            double taus[2] = { tlo, thi }, bestm = HUGE_VAL, bestt = 0.0;
            for (int cc = 0; cc < 2; cc++) {
                for (size_t i = 0; i < n; i++) wv[i] = p[i] + taus[cc] * zev[i];
                fm_matvec(B, n, wv, Bp);
                double m = fm_dot(g, wv, n) + 0.5 * fm_dot(wv, Bp, n);
                if (m < bestm) { bestm = m; bestt = taus[cc]; }
            }
            for (size_t i = 0; i < n; i++) p[i] += bestt * zev[i];
            *hits = true; ret = true; goto done;
        }
        double lam_new = lam + (pnorm * pnorm / wn2) * ((pnorm - Delta) / Delta);
        if (!(lam_new > lam_lo && (lam_hi == HUGE_VAL || lam_new < lam_hi))) {
            lam_new = (lam_hi < HUGE_VAL) ? 0.5 * (lam_lo + lam_hi)
                                          : (lam > 0.0 ? lam * 2.0 : 1.0);
        }
        lam = lam_new;
    }
    if (have_p) { *hits = true; ret = true; }           /* cap: best-effort step */
done:
    free(L); free(negg); free(wv); free(zev); free(Bp);
    return ret;
}

/* ---- trust-krylov: GLTR (Hessian-vector only) -------------------------- *
 * Lanczos-tridiagonalize B in the Krylov space span{g, Bg, B²g, ...} with full
 * re-orthogonalization, and at each expansion solve the small tridiagonal
 * trust-region subproblem exactly via fm_tr_moresorensen. Reconstruct p = Q·y.
 * Unlike Steihaug-CG, the exact tridiagonal solve recovers the least-eigenvector
 * component, so the indefinite/hard case is handled as the Krylov space grows. */
bool fm_tr_gltr(const FmQuad* q, const double* g, double gnorm,
                       double Delta, double* p, bool* hits) {
    size_t n = q->n;
    *hits = false;
    if (gnorm == 0.0) { for (size_t i = 0; i < n; i++) p[i] = 0.0; return true; }
    size_t KMAX = n; if (KMAX > 100) KMAX = 100; if (KMAX < 1) KMAX = 1;

    double* Q     = (double*)malloc(n * KMAX * sizeof(double));
    double* alpha = (double*)malloc(KMAX * sizeof(double));
    double* betav = (double*)malloc(KMAX * sizeof(double));
    double* w     = (double*)malloc(n * sizeof(double));
    double* Hv    = (double*)malloc(n * sizeof(double));
    double* y     = (double*)malloc(KMAX * sizeof(double));
    double* T     = (double*)malloc(KMAX * KMAX * sizeof(double));
    double* gsub  = (double*)malloc(KMAX * sizeof(double));
    bool ret = false;
    if (!Q || !alpha || !betav || !w || !Hv || !y || !T || !gsub) goto done;

    for (size_t i = 0; i < n; i++) Q[i] = g[i] / gnorm;  /* q_1 = g/||g||        */
    size_t k = 0; double beta_prev = 0.0;

    for (size_t j = 0; j < KMAX; j++) {
        const double* qj = Q + j * n;
        if (!fm_quad_matvec(q, qj, Hv)) break;
        double aj = fm_dot(qj, Hv, n);
        alpha[j] = aj;
        for (size_t i = 0; i < n; i++) w[i] = Hv[i] - aj * qj[i];
        if (j > 0) { const double* qm = Q + (j - 1) * n;
                     for (size_t i = 0; i < n; i++) w[i] -= beta_prev * qm[i]; }
        for (size_t t = 0; t <= j; t++) {               /* full re-orthogonalize */
            const double* qt = Q + t * n;
            double c = fm_dot(qt, w, n);
            for (size_t i = 0; i < n; i++) w[i] -= c * qt[i];
        }
        double bj = sqrt(fm_dot(w, w, n));
        k = j + 1;

        /* Tridiagonal subproblem: min ||g|| e_1·y + ½ yᵀT_k y, ||y|| ≤ Δ. */
        for (size_t a = 0; a < k * k; a++) T[a] = 0.0;
        for (size_t a = 0; a < k; a++) T[a * k + a] = alpha[a];
        for (size_t a = 0; a + 1 < k; a++) { T[a * k + (a + 1)] = betav[a];
                                             T[(a + 1) * k + a] = betav[a]; }
        for (size_t a = 0; a < k; a++) gsub[a] = 0.0;
        gsub[0] = gnorm;
        FmQuad qT; qT.n = k; qT.B = T; qT.c = NULL; qT.xbase = NULL; qT.gbase = NULL;
        qT.xpert = NULL; qT.gpert = NULL; qT.active0 = NULL;
        bool hitsT = false;
        if (!fm_tr_moresorensen(&qT, gsub, gnorm, Delta, y, &hitsT)) break;
        *hits = hitsT;

        betav[j] = bj; beta_prev = bj;
        double resid = bj * fabs(y[k - 1]);
        if (bj <= 1e-12 * (1.0 + fabs(aj)) || resid < 1e-9 * gnorm
            || k >= n || j == KMAX - 1) { ret = true; break; }
        double* qn = Q + (j + 1) * n;
        for (size_t i = 0; i < n; i++) qn[i] = w[i] / bj;
    }
    /* ret is true iff the loop broke after a successful tridiagonal solve, so y
     * is valid exactly then; on a matvec/subproblem failure ret stays false and
     * the driver falls back to the best iterate rather than reading a stale y. */
    if (ret && k > 0) {
        for (size_t i = 0; i < n; i++) p[i] = 0.0;       /* p = Q_k y            */
        for (size_t a = 0; a < k; a++) { const double* qa = Q + a * n; double ya = y[a];
            for (size_t i = 0; i < n; i++) p[i] += ya * qa[i]; }
    } else {
        ret = false;
    }
done:
    free(Q); free(alpha); free(betav); free(w); free(Hv); free(y); free(T); free(gsub);
    return ret;
}

/* ---- shared trust-region driver --------------------------------------- */
bool fm_run_trust_region(Expr* f, Expr** vars, size_t n,
                                FmVarBind* binds, Expr** g_exprs, Expr*** H_exprs,
                                double* x, /* in/out */ const FmBox* boxes,
                                const FmOpts* opts, double* fx_out,
                                FmSubSolver solve_sub, bool needs_dense_B) {
    (void)vars;
    FmLbfgsCtx ctx;
    ctx.f = f; ctx.g_exprs = g_exprs; ctx.binds = binds; ctx.n = n;
    ctx.gens = NULL; ctx.ngens = 0; ctx.mu = 0.0; ctx.boxes = NULL;
    ctx.opts = opts; ctx.augmented = false;

    double* g     = (double*)malloc(n * sizeof(double));
    double* p     = (double*)malloc(n * sizeof(double));
    double* x_new = (double*)malloc(n * sizeof(double));
    double* Bp    = (double*)malloc(n * sizeof(double));
    double* xpert = (double*)malloc(n * sizeof(double));
    double* gpert = (double*)malloc(n * sizeof(double));
    bool*   act0  = (bool*)malloc(n * sizeof(bool));
    double* B     = needs_dense_B ? (double*)malloc(n * n * sizeof(double)) : NULL;
    double* ej    = needs_dense_B ? (double*)malloc(n * sizeof(double)) : NULL;
    double* hcol  = needs_dense_B ? (double*)malloc(n * sizeof(double)) : NULL;
    bool ok = false; double fx = 0.0;
    if (!g || !p || !x_new || !Bp || !xpert || !gpert || !act0
        || (needs_dense_B && (!B || !ej || !hcol))) goto cleanup;
    for (size_t i = 0; i < n; i++) act0[i] = false;

    if (boxes) fm_project_box(x, n, boxes);
    if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) {
        fm_warn(g_fm_name, "nlnum", "objective evaluation failed at start point"); goto cleanup;
    }
    if (!fm_tnc_grad(&ctx, x, g)) {
        fm_warn(g_fm_name, "nlnum", "gradient evaluation failed at start point"); goto cleanup;
    }
    ok = true;                                           /* have a usable iterate */

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    double Delta = FM_TR_DELTA0;

    FmQuad q; q.n = n; q.B = B; q.c = &ctx; q.xbase = x; q.gbase = g;
    q.xpert = xpert; q.gpert = gpert; q.active0 = act0;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        double gnorm = sqrt(fm_dot(g, g, n));
        if (gnorm < tol_acc) break;
        if (needs_dense_B) {
            if (!fm_tr_build_hessian(&ctx, H_exprs, x, g, B, ej, hcol, xpert, gpert, act0)) {
                fm_warn(g_fm_name, "nlnum", "Hessian evaluation failed"); break;
            }
        }
        q.xbase = x; q.gbase = g;
        bool hits = false;
        if (!solve_sub(&q, g, gnorm, Delta, p, &hits)) break;

        double gp = fm_dot(g, p, n);
        double pBp = 0.0;
        if (fm_quad_matvec(&q, p, Bp)) pBp = fm_dot(p, Bp, n);
        double pred = -(gp + 0.5 * pBp);

        for (size_t i = 0; i < n; i++) x_new[i] = x[i] + p[i];
        double fx_new;
        bool feval = fm_eval_scalar(f, binds, x_new, n, opts, &fx_new);
        if (!feval || !isfinite(fx_new)) {
            Delta *= 0.25;                               /* infeasible trial      */
        } else {
            double act = fx - fx_new;
            double tiny = 1e-16 * (1.0 + fabs(fx));
            double rho = (pred <= tiny) ? (act > 0.0 ? 1.0 : -1.0) : act / pred;
            if (rho < 0.25) Delta *= 0.25;
            else if (rho > 0.75 && hits) { Delta *= 2.0; if (Delta > FM_TR_DELTAMAX) Delta = FM_TR_DELTAMAX; }
            if (rho > FM_TR_ETA) {                       /* accept                */
                double max_step = 0.0, max_x = 0.0;
                for (size_t i = 0; i < n; i++) {
                    double ds = fabs(x_new[i] - x[i]); if (ds > max_step) max_step = ds;
                    double ax = fabs(x_new[i]);         if (ax > max_x)    max_x    = ax;
                    x[i] = x_new[i];
                }
                fx = fx_new;
                if (!fm_tnc_grad(&ctx, x, g)) {
                    fm_warn(g_fm_name, "nlnum", "gradient evaluation failed in iteration"); break;
                }
                fm_fire_monitor(opts->step_monitor);
                if (max_step < tol_prec * (max_x + 1e-300)) break;
                continue;
            }
        }
        double xinf = 0.0;
        for (size_t i = 0; i < n; i++) { double ax = fabs(x[i]); if (ax > xinf) xinf = ax; }
        if (Delta < tol_prec * (xinf + 1e-300)) break;   /* radius underflow      */
    }
    *fx_out = fx;
cleanup:
    free(g); free(p); free(x_new); free(Bp); free(xpert); free(gpert); free(act0);
    if (B)    free(B);
    if (ej)   free(ej);
    if (hcol) free(hcol);
    return ok;
}

/* ---- Newton-CG (line search, not trust region) ------------------------- *
 * Inexact Newton: fm_tnc_cg approximately solves B p = -g with the
 * Eisenstat-Walker forcing sequence and negative-curvature truncation, then a
 * unit-step-first Wolfe line search (fm_lbfgs_linesearch, NOT fm_line_search
 * whose 1/||d|| cap throttles a well-scaled Newton step). Hessian-free. */
bool fm_run_newton_cg(Expr* f, Expr** vars, size_t n, FmVarBind* binds,
                             Expr** g_exprs, double* x, /* in/out */
                             const FmGenCon* gens, size_t ngens, double mu,
                             const FmBox* boxes, const FmOpts* opts, double* fx_out) {
    (void)vars; (void)gens; (void)ngens; (void)mu;
    FmLbfgsCtx ctx;
    ctx.f = f; ctx.g_exprs = g_exprs; ctx.binds = binds; ctx.n = n;
    ctx.gens = NULL; ctx.ngens = 0; ctx.mu = 0.0; ctx.boxes = NULL;
    ctx.opts = opts; ctx.augmented = false;

    double* g     = (double*)malloc(n * sizeof(double));
    double* p     = (double*)malloc(n * sizeof(double));
    double* r     = (double*)malloc(n * sizeof(double));
    double* d     = (double*)malloc(n * sizeof(double));
    double* Hd    = (double*)malloc(n * sizeof(double));
    double* xpert = (double*)malloc(n * sizeof(double));
    double* gpert = (double*)malloc(n * sizeof(double));
    double* x_new = (double*)malloc(n * sizeof(double));
    double* g_new = (double*)malloc(n * sizeof(double));
    bool*   act0  = (bool*)malloc(n * sizeof(bool));
    bool ok = false; double fx = 0.0;
    if (!g || !p || !r || !d || !Hd || !xpert || !gpert || !x_new || !g_new || !act0)
        goto cleanup;
    for (size_t i = 0; i < n; i++) act0[i] = false;

    if (boxes) fm_project_box(x, n, boxes);
    if (!fm_lbfgs_fg(&ctx, x, &fx, g)) {
        fm_warn(g_fm_name, "nlnum", "objective/gradient evaluation failed at start point");
        goto cleanup;
    }
    ok = true;

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        double gnorm = sqrt(fm_dot(g, g, n));
        if (gnorm < tol_acc) break;
        fm_tnc_cg(&ctx, x, g, g, act0, p, r, d, Hd, xpert, gpert);   /* gm = g */
        double dphi0 = fm_dot(g, p, n);
        if (dphi0 >= 0.0) {                              /* FD noise: steepest    */
            for (size_t i = 0; i < n; i++) p[i] = -g[i];
            dphi0 = fm_dot(g, p, n);
            if (dphi0 >= 0.0) break;
        }
        double a, fx_new;
        if (!fm_lbfgs_linesearch(&ctx, x, p, fx, dphi0, HUGE_VAL,
                                 x_new, g_new, &fx_new, &a)) {
            fm_warn(g_fm_name, "lstol", "line search failed at iter %lld", (long long)k);
            break;
        }
        fm_fire_monitor(opts->step_monitor);
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]); if (ds > max_step) max_step = ds;
            double ax = fabs(x_new[i]);        if (ax > max_x)    max_x    = ax;
            x[i] = x_new[i]; g[i] = g_new[i];
        }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) break;
    }
    *fx_out = fx;
cleanup:
    free(g); free(p); free(r); free(d); free(Hd);
    free(xpert); free(gpert); free(x_new); free(g_new); free(act0);
    return ok;
}
