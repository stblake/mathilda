/* findmin_cobyqa.c — COBYQA constrained derivative-free local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  COBYQA -- Constrained Optimization BY Quadratic Approximations       *
 * ------------------------------------------------------------------ *
 * A derivative-free trust-region SQP with QUADRATIC interpolation models, the
 * analogue of scipy's minimize(method="COBYQA") (Ragonneau & Zhang 2023). Where
 * COBYLA models f and the constraints by linear approximations, COBYQA models
 * each by a full quadratic -- so it captures curvature (converging tighter on
 * smooth problems, and navigating curved valleys that COBYLA's linear models
 * cannot -- e.g. Rosenbrock) and handles equality + inequality + bound
 * constraints natively.
 *
 * Each iteration builds the quadratic models by finite differences on a
 * structured stencil around the base of side Delta: the coordinate cross
 * base +/- Delta e_i gives the gradient and the diagonal Hessian in closed form,
 * and the corner points base + Delta(e_i + e_j) give the off-diagonal Hessian
 * entries -- a fully-determined quadratic per function, no interpolation-matrix
 * bookkeeping (a minimum-Frobenius-norm 2n+1-point model is the eval-efficiency
 * refinement). The Byrd-Omojokun trust-region SQP step is realised as a single
 * inf-norm box trust-region QP handed to the SLSQP dual active-set solver
 * (fm_slsqp_activeset) with the Lagrangian-Hessian model B = H_f + sum lambda_k
 * H_{c_k}; steps are accepted by an L1 exact-penalty merit with Powell's penalty
 * update, and the trust radius Delta doubles as the stencil scale (grown on a
 * good step, halved on a poor one, terminating at Delta = 10^-PrecisionGoal).
 * Machine precision only (WorkingPrecision > MachinePrecision falls back to
 * QuasiNewton). Exposed as Method -> "COBYQA". Reference: Ragonneau & Zhang 2023;
 * Conn, Gould & Toint 2000 (derivative-free trust regions). */

/* Max constraint violation of the general constraints at a point whose raw gens
 * values are gvals[ngens] (equalities counted as |value|, inequalities as
 * max(0,value)); bounds are handled by projection, not here. */
static double fm_cobyqa_viol(const double* gvals, const FmGenCon* gens, size_t ngens) {
    double mv = 0.0;
    for (size_t k = 0; k < ngens; k++) {
        double v = gens[k].equality ? fabs(gvals[k]) : (gvals[k] > 0.0 ? gvals[k] : 0.0);
        if (v > mv) mv = v;
    }
    return mv;
}

/* Trust-region SQP step: min 1/2 dᵀ HL d + gf·d s.t. the linearised constraints
 * and ||d||_inf <= Delta, solved by the SLSQP dual active-set QP. HL is the
 * Lagrangian-Hessian model (regularised to SPD by a tau-retry, reset to I if it
 * stays indefinite). cvals[k]/Jac[k*n..] are each constraint's model value and
 * gradient at the base. Writes d[n] and the signed multipliers lam[ngens]. */
static bool fm_cobyqa_subproblem(size_t n, const double* HL, const double* gf,
                                 const FmGenCon* gens, const double* cvals,
                                 const double* Jac, size_t ngens,
                                 const FmBox* boxes, double Delta,
                                 double* d, double* lam) {
    size_t nb = 0;
    for (size_t i = 0; i < n; i++) {
        if (boxes && boxes[i].has_lo) nb++;
        if (boxes && boxes[i].has_hi) nb++;
    }
    size_t m = ngens + nb + 2 * n;
    double* L    = (double*)calloc(n * n, sizeof(double));
    double* Nrm  = (double*)calloc(m * n, sizeof(double));
    double* b    = (double*)malloc(sizeof(double) * m);
    int*    ie   = (int*)calloc(m, sizeof(int));
    double* mult = (double*)malloc(sizeof(double) * m);
    bool ok = false;
    if (!L || !Nrm || !b || !ie || !mult) goto done;

    /* B = HL + tau I -> Cholesky factor L (tau-retry; identity fallback). */
    {
        double tr = 0.0;
        for (size_t i = 0; i < n; i++) tr += fabs(HL[i * n + i]);
        double tau = 0.0, tau0 = 1e-10 * (1.0 + tr / (double)(n ? n : 1));
        bool fac = false;
        for (int attempt = 0; attempt < 6 && !fac; attempt++) {
            for (size_t t = 0; t < n * n; t++) L[t] = HL[t];
            if (fm_chol_factor(L, n, tau)) fac = true;
            else tau = (tau == 0.0) ? tau0 : tau * 100.0;
        }
        if (!fac) {
            for (size_t t = 0; t < n * n; t++) L[t] = 0.0;
            for (size_t i = 0; i < n; i++) L[i * n + i] = 1.0;
        }
    }

    /* General constraints: equality Jac_k·d = -cval (is_eq); inequality c<=0
     * linearised (-Jac_k)·d >= cval. (Same mapping as fm_slsqp_qp.) */
    for (size_t k = 0; k < ngens; k++) {
        double* row = &Nrm[k * n];
        if (gens[k].equality) {
            for (size_t i = 0; i < n; i++) row[i] = Jac[k * n + i];
            b[k] = -cvals[k]; ie[k] = 1;
        } else {
            for (size_t i = 0; i < n; i++) row[i] = -Jac[k * n + i];
            b[k] = cvals[k]; ie[k] = 0;
        }
    }
    /* Bounds. The caller passes a base-relative box (lo-x0, hi-x0), so these are
     * already the step-frame residuals: e_i·d >= lo-x0_i and -e_i·d >= x0_i-hi. */
    {
        size_t r = ngens;
        for (size_t i = 0; i < n; i++) {
            if (boxes && boxes[i].has_lo) { Nrm[r * n + i] = 1.0;  b[r] = boxes[i].lo; r++; }
            if (boxes && boxes[i].has_hi) { Nrm[r * n + i] = -1.0; b[r] = -boxes[i].hi; r++; }
        }
    }
    /* Trust box: -Delta <= d_i <= Delta. */
    {
        size_t r = ngens + nb;
        for (size_t i = 0; i < n; i++) {
            Nrm[(r) * n + i]     = 1.0;  b[r]     = -Delta;
            Nrm[(r + 1) * n + i] = -1.0; b[r + 1] = -Delta;
            r += 2;
        }
    }

    (void)fm_slsqp_activeset(L, n, gf, Nrm, b, ie, m, d, mult);
    for (size_t k = 0; k < ngens; k++)
        lam[k] = gens[k].equality ? -mult[k] : mult[k];
    ok = true;
done:
    free(L); free(Nrm); free(b); free(ie); free(mult);
    return ok;
}

bool fm_run_cobyqa(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars; (void)g_exprs; (void)mu;

    if (n == 1 && ngens == 0) {
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

    size_t ngv = ngens ? ngens : 1;
    double* base   = (double*)malloc(sizeof(double) * n);
    double* gbase  = (double*)malloc(sizeof(double) * ngv);
    double* gf     = (double*)malloc(sizeof(double) * n);
    double* Hf     = (double*)malloc(sizeof(double) * n * n);
    double* gck    = (double*)malloc(sizeof(double) * ngv * n);
    double* Hck    = (double*)malloc(sizeof(double) * ngv * n * n);
    double* HL     = (double*)malloc(sizeof(double) * n * n);
    double* lam    = (double*)calloc(ngv, sizeof(double));
    double* pen    = (double*)calloc(ngv, sizeof(double));
    double* xs     = (double*)malloc(sizeof(double) * n);
    double* gp     = (double*)malloc(sizeof(double) * ngv);
    double* gm     = (double*)malloc(sizeof(double) * ngv);
    double* gcn    = (double*)malloc(sizeof(double) * ngv);
    double* d      = (double*)malloc(sizeof(double) * n);
    double* xtrial = (double*)malloc(sizeof(double) * n);
    double* gtrial = (double*)malloc(sizeof(double) * ngv);
    double* xbest  = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    double fbase = 0.0;
    bool have_best = false; double best_f = 0.0;
    if (!base || !gbase || !gf || !Hf || !gck || !Hck || !HL || !lam || !pen
        || !xs || !gp || !gm || !gcn || !d || !xtrial || !gtrial || !xbest)
        goto cleanup;

    if (boxes) fm_project_box(x, n, boxes);
    for (size_t i = 0; i < n; i++) base[i] = x[i];
    if (!fm_cobyla_eval(f, gens, ngens, binds, opts, base, n, &fbase, gbase)) {
        fm_warn(g_fm_name, "nlnum", "objective/constraint evaluation failed at start point");
        goto cleanup;
    }

    double dend = pow(10.0, -opts->prec_goal_digits);
    double dbeg = 1.0;
    for (size_t i = 0; i < n; i++)
        if (boxes && boxes[i].has_lo && boxes[i].has_hi) {
            double half = 0.5 * (boxes[i].hi - boxes[i].lo);
            if (half > 0.0 && half < dbeg) dbeg = half;
        }
    if (dbeg < 4.0 * dend) dbeg = 4.0 * dend;
    double Delta = dbeg;
    double Delta_max = 10.0 * dbeg;

    for (int64_t it = 0; it < opts->max_iter; it++) {
        /* Build quadratic models by finite differences on the stencil of side
         * Delta.  Samples are not forced inside the box (algebraic objectives are
         * defined everywhere); a non-numeric sample shrinks Delta and retries. */
        bool model_ok = true;
        for (size_t i = 0; i < n && model_ok; i++) {
            double fp, fmv;
            for (size_t t = 0; t < n; t++) xs[t] = base[t];
            xs[i] = base[i] + Delta;
            if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xs, n, &fp, gp)) { model_ok = false; break; }
            xs[i] = base[i] - Delta;
            if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xs, n, &fmv, gm)) { model_ok = false; break; }
            gf[i] = (fp - fmv) / (2.0 * Delta);
            Hf[i * n + i] = (fp + fmv - 2.0 * fbase) / (Delta * Delta);
            for (size_t k = 0; k < ngens; k++) {
                gck[k * n + i] = (gp[k] - gm[k]) / (2.0 * Delta);
                Hck[k * n * n + i * n + i] = (gp[k] + gm[k] - 2.0 * gbase[k]) / (Delta * Delta);
            }
        }
        for (size_t i = 0; i < n && model_ok; i++) {
            for (size_t j = i + 1; j < n; j++) {
                double fc;
                for (size_t t = 0; t < n; t++) xs[t] = base[t];
                xs[i] = base[i] + Delta; xs[j] = base[j] + Delta;
                if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xs, n, &fc, gcn)) { model_ok = false; break; }
                double hij = (fc - fbase - Delta * (gf[i] + gf[j])
                              - 0.5 * Delta * Delta * (Hf[i * n + i] + Hf[j * n + j]))
                             / (Delta * Delta);
                Hf[i * n + j] = Hf[j * n + i] = hij;
                for (size_t k = 0; k < ngens; k++) {
                    double h = (gcn[k] - gbase[k] - Delta * (gck[k * n + i] + gck[k * n + j])
                                - 0.5 * Delta * Delta * (Hck[k * n * n + i * n + i] + Hck[k * n * n + j * n + j]))
                               / (Delta * Delta);
                    Hck[k * n * n + i * n + j] = Hck[k * n * n + j * n + i] = h;
                }
            }
        }
        if (!model_ok) {
            if (Delta <= dend) { ok = true; break; }
            Delta *= 0.5; if (Delta <= 1.5 * dend) Delta = dend;
            continue;
        }

        /* Lagrangian-Hessian model H_L = H_f + sum lambda_k H_{c_k}. */
        for (size_t t = 0; t < n * n; t++) {
            double v = Hf[t];
            for (size_t k = 0; k < ngens; k++) v += lam[k] * Hck[k * n * n + t];
            HL[t] = v;
        }

        /* Bound rows in the subproblem need residuals relative to the base; pass
         * a shifted box so `lo <= x0+d` becomes `d_i >= lo - x0_i`. */
        FmBox* rbox = NULL;
        if (boxes) {
            rbox = (FmBox*)malloc(sizeof(FmBox) * n);
            if (!rbox) goto cleanup;
            for (size_t i = 0; i < n; i++) {
                rbox[i].has_lo = boxes[i].has_lo; rbox[i].has_hi = boxes[i].has_hi;
                rbox[i].lo = boxes[i].lo - base[i];
                rbox[i].hi = boxes[i].hi - base[i];
            }
        }
        bool sp_ok = fm_cobyqa_subproblem(n, HL, gf, gens, gbase, gck, ngens,
                                          rbox, Delta, d, lam);
        free(rbox);
        if (!sp_ok) goto cleanup;

        double dnorm = 0.0;
        for (size_t i = 0; i < n; i++) if (fabs(d[i]) > dnorm) dnorm = fabs(d[i]);
        if (dnorm < 1e-3 * Delta) {
            if (Delta <= dend) { ok = true; break; }
            Delta *= 0.5; if (Delta <= 1.5 * dend) Delta = dend;
            continue;
        }

        /* Penalty update (Powell), then merit ratio. */
        for (size_t k = 0; k < ngens; k++) {
            double a = fabs(lam[k]);
            double avg = 0.5 * (pen[k] + a);
            double np = (a > avg) ? a : avg;
            if (np < a) np = a;
            if (np > 1e12) np = 1e12;
            pen[k] = np;
        }
        double base_viol = fm_cobyqa_viol(gbase, gens, ngens);
        /* predicted objective reduction -(gf·d + 1/2 dᵀHf d) */
        double gfd = 0.0, dHd = 0.0;
        for (size_t i = 0; i < n; i++) {
            gfd += gf[i] * d[i];
            double t = 0.0;
            for (size_t j = 0; j < n; j++) t += Hf[i * n + j] * d[j];
            dHd += d[i] * t;
        }
        double obj_red = -(gfd + 0.5 * dHd);
        double lin_viol = 0.0;
        for (size_t k = 0; k < ngens; k++) {
            double c = gbase[k];
            for (size_t i = 0; i < n; i++) c += gck[k * n + i] * d[i];
            double v = gens[k].equality ? fabs(c) : (c > 0.0 ? c : 0.0);
            if (v > lin_viol) lin_viol = v;
        }
        /* Merit uses a single penalty scale (max pen) on the L-infinity violation. */
        double penmax = 0.0;
        for (size_t k = 0; k < ngens; k++) if (pen[k] > penmax) penmax = pen[k];
        double pred = obj_red + penmax * (base_viol - lin_viol);

        for (size_t i = 0; i < n; i++) xtrial[i] = base[i] + d[i];
        if (boxes) fm_project_box(xtrial, n, boxes);
        double ftrial;
        if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xtrial, n, &ftrial, gtrial))
            goto cleanup;
        fm_fire_monitor(opts->step_monitor);

        double trial_viol = fm_cobyqa_viol(gtrial, gens, ngens);
        double ared = (fbase + penmax * base_viol) - (ftrial + penmax * trial_viol);
        double r = (pred > 1e-300) ? (ared / pred) : (ared > 0.0 ? 1.0 : -1.0);

        if (trial_viol <= 1e-8 && (!have_best || ftrial < best_f)) {
            for (size_t i = 0; i < n; i++) xbest[i] = xtrial[i];
            best_f = ftrial; have_best = true;
        }

        if (r > 0.1 && ared > 0.0) {
            for (size_t i = 0; i < n; i++) base[i] = xtrial[i];
            fbase = ftrial;
            for (size_t k = 0; k < ngens; k++) gbase[k] = gtrial[k];
            if (r > 0.7 && dnorm > 0.5 * Delta) {
                Delta *= 2.0; if (Delta > Delta_max) Delta = Delta_max;
            }
        } else {
            if (Delta <= dend) { ok = true; break; }
            Delta *= 0.5; if (Delta <= 1.5 * dend) Delta = dend;
        }
    }

    {
        double base_viol = fm_cobyqa_viol(gbase, gens, ngens);
        if (base_viol > 1e-8 && have_best) {
            for (size_t i = 0; i < n; i++) base[i] = xbest[i];
            fbase = best_f;
        } else if (base_viol > 1e-6) {
            fm_warn(g_fm_name, "infeas", "could not satisfy constraints to tolerance");
        }
    }
    for (size_t i = 0; i < n; i++) x[i] = base[i];
    if (boxes) fm_project_box(x, n, boxes);
    *fx_out = fbase;
    ok = true;
cleanup:
    free(base); free(gbase); free(gf); free(Hf); free(gck); free(Hck); free(HL);
    free(lam); free(pen); free(xs); free(gp); free(gm); free(gcn);
    free(d); free(xtrial); free(gtrial); free(xbest);
    return ok;
}
