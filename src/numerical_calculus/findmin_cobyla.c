/* findmin_cobyla.c — COBYLA constrained derivative-free local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  COBYLA -- Constrained Optimization BY Linear Approximation           *
 * ------------------------------------------------------------------ *
 * Powell's derivative-free (function-values-only) trust-region method for
 * CONSTRAINED optimization, the analogue of scipy's minimize(method="COBYLA").
 * It is the FIRST derivative-free method in Mathilda to accept general
 * (non-box) constraints -- Powell/NelderMead reject them, and SLSQP/TNC need a
 * gradient -- so it is the method for non-smooth / noisy / black-box objectives
 * WITH constraints.
 *
 * At each iterate it models f and every constraint by an AFFINE (linear)
 * approximation, then solves a trust-region subproblem: (stage 1) minimise the
 * worst constraint violation, then (stage 2) minimise the linear objective model
 * while holding that minimum violation -- Powell's two-stage `trstlp`. Here the
 * linear models are recovered by finite differences on a coordinate cross
 * (exact for the trust radius, and immune to the simplex-degeneracy bookkeeping
 * that dominates a literal cobyla.f port), and the two-stage LP is solved by the
 * same dual active-set QP as SLSQP (`fm_slsqp_activeset`) with a tiny
 * regularisation and an inf-norm trust box -- reaching the same constrained
 * optimum as reference COBYLA. Step acceptance uses Powell's L-infinity
 * exact-penalty merit `Phi = f + mu*max_i max(0, c_i)` with his PARMU update.
 * Equalities `h==0` are handled by splitting into `h<=0` and `-h<=0` (so this is
 * strictly more capable than scipy's inequality-only COBYLA); box bounds enter
 * as ordinary linear inequalities. Machine precision only (WorkingPrecision >
 * MachinePrecision falls back to QuasiNewton). Exposed as Method -> "COBYLA".
 * Reference: Powell 1994. */

/* Two-stage inf-norm trust-region LP for the COBYLA step, solved as a pair of
 * strictly-convex QPs by fm_slsqp_activeset (B = eps*I, so the quadratic term is
 * negligible and each solve is an LP to rounding). The mcon internal constraints
 * are in Mathilda's `c <= 0` feasible convention: row j is the affine model
 * `cval0[j] + Amod[j,:]*d`. Stage 1 finds t* = min over ||d||_inf <= rho of the
 * worst model violation max_j (cval0[j]+Amod[j]*d); stage 2 minimises af*d subject
 * to every model constraint staying <= max(t*,0). Writes the step into d. */
static bool fm_cobyla_subproblem(size_t n, size_t mcon, const double* Amod,
                                 const double* cval0, const double* af,
                                 double rho, double* d) {
    const double eps = 1e-8;
    size_t nv1 = n + 1;                 /* stage-1 variables (d, t) */
    size_t m1  = mcon + 2 * n;          /* constraint rows + 2n trust-box rows */
    double* L1   = (double*)calloc(nv1 * nv1, sizeof(double));
    double* g1   = (double*)calloc(nv1, sizeof(double));
    double* Nrm1 = (double*)calloc(m1 * nv1, sizeof(double));
    double* b1   = (double*)malloc(sizeof(double) * m1);
    int*    ie1  = (int*)calloc(m1, sizeof(int));
    double* mul1 = (double*)malloc(sizeof(double) * m1);
    double* y1   = (double*)malloc(sizeof(double) * nv1);
    double* L2   = (double*)calloc(n * n, sizeof(double));
    double* Nrm2 = (double*)calloc((mcon + 2 * n) * n, sizeof(double));
    double* b2   = (double*)malloc(sizeof(double) * (mcon + 2 * n));
    int*    ie2  = (int*)calloc(mcon + 2 * n, sizeof(int));
    double* mul2 = (double*)malloc(sizeof(double) * (mcon + 2 * n));
    bool ok = false;
    if (!L1 || !g1 || !Nrm1 || !b1 || !ie1 || !mul1 || !y1
        || !L2 || !Nrm2 || !b2 || !ie2 || !mul2) goto done;

    /* Stage 1: min t + eps/2*||(d,t)||^2. */
    for (size_t i = 0; i < nv1; i++) L1[i * nv1 + i] = sqrt(eps);
    g1[n] = 1.0;
    /* model rows: (-Amod_j)*d + 1*t >= cval0_j  (i.e. cval0_j + Amod_j*d <= t) */
    for (size_t j = 0; j < mcon; j++) {
        double* row = &Nrm1[j * nv1];
        for (size_t i = 0; i < n; i++) row[i] = -Amod[j * n + i];
        row[n] = 1.0;
        b1[j] = cval0[j];
    }
    /* trust box: d_i >= -rho and -d_i >= -rho */
    for (size_t i = 0; i < n; i++) {
        double* rlo = &Nrm1[(mcon + 2 * i) * nv1];
        double* rhi = &Nrm1[(mcon + 2 * i + 1) * nv1];
        rlo[i] = 1.0;  b1[mcon + 2 * i]     = -rho;
        rhi[i] = -1.0; b1[mcon + 2 * i + 1] = -rho;
    }
    (void)fm_slsqp_activeset(L1, nv1, g1, Nrm1, b1, ie1, m1, y1, mul1);
    double tstar = y1[n];
    double tcap = (tstar > 0.0) ? tstar : 0.0;

    /* Stage 2: min af*d + eps/2*||d||^2 s.t. cval0_j + Amod_j*d <= tcap. */
    for (size_t i = 0; i < n; i++) L2[i * n + i] = sqrt(eps);
    for (size_t j = 0; j < mcon; j++) {
        double* row = &Nrm2[j * n];
        for (size_t i = 0; i < n; i++) row[i] = -Amod[j * n + i];
        b2[j] = cval0[j] - tcap;
    }
    for (size_t i = 0; i < n; i++) {
        double* rlo = &Nrm2[(mcon + 2 * i) * n];
        double* rhi = &Nrm2[(mcon + 2 * i + 1) * n];
        rlo[i] = 1.0;  b2[mcon + 2 * i]     = -rho;
        rhi[i] = -1.0; b2[mcon + 2 * i + 1] = -rho;
    }
    (void)fm_slsqp_activeset(L2, n, af, Nrm2, b2, ie2, mcon + 2 * n, d, mul2);
    ok = true;
done:
    free(L1); free(g1); free(Nrm1); free(b1); free(ie1); free(mul1); free(y1);
    free(L2); free(Nrm2); free(b2); free(ie2); free(mul2);
    return ok;
}

/* Max constraint violation of the internal (c<=0) list at a point whose raw
 * gens values are gv[ngens] and coordinates p[n]. */
static double fm_cobyla_maxviol(const double* gv, const double* p, size_t n,
                                const size_t* gi_k, const double* gi_sgn, size_t ng,
                                const size_t* bnd_var, const bool* bnd_lo, size_t nb,
                                const FmBox* boxes) {
    (void)n;
    double mv = 0.0;
    for (size_t gi = 0; gi < ng; gi++) {
        double v = gi_sgn[gi] * gv[gi_k[gi]];
        if (v > mv) mv = v;
    }
    for (size_t bi = 0; bi < nb; bi++) {
        size_t i = bnd_var[bi];
        double v = bnd_lo[bi] ? (boxes[i].lo - p[i]) : (p[i] - boxes[i].hi);
        if (v > mv) mv = v;
    }
    return mv;
}

/* Evaluate f and every raw constraint gens[k].expr at p; fills *fout and
 * gv[ngens]. Returns false if any evaluation is non-numeric. */
bool fm_cobyla_eval(Expr* f, const FmGenCon* gens, size_t ngens,
                           FmVarBind* binds, const FmOpts* opts,
                           const double* p, size_t n, double* fout, double* gv) {
    if (!fm_eval_scalar(f, binds, p, n, opts, fout)) return false;
    for (size_t k = 0; k < ngens; k++)
        if (!fm_eval_scalar(gens[k].expr, binds, p, n, opts, &gv[k])) return false;
    return true;
}

bool fm_run_cobyla(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars; (void)g_exprs; (void)mu;

    /* n==1 with no general constraints -> the exact Brent path (like
     * Powell/NelderMead); general constraints at n==1 use the full machinery. */
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

    /* Internal constraint list: each inequality -> 1 entry; each equality ->
     * 2 entries (+h<=0, -h<=0); each active bound -> 1 entry (exact +/-e_i). */
    size_t n_eq = 0;
    for (size_t k = 0; k < ngens; k++) if (gens[k].equality) n_eq++;
    size_t ng = ngens + n_eq;
    size_t nb = 0;
    for (size_t i = 0; i < n; i++) {
        if (boxes && boxes[i].has_lo) nb++;
        if (boxes && boxes[i].has_hi) nb++;
    }
    size_t mcon = ng + nb;
    size_t ng1 = ng ? ng : 1, mc1 = mcon ? mcon : 1, ngv = ngens ? ngens : 1;

    size_t* gi_k   = (size_t*)malloc(sizeof(size_t) * ng1);
    double* gi_sgn = (double*)malloc(sizeof(double) * ng1);
    size_t* bnd_var = (size_t*)malloc(sizeof(size_t) * (nb ? nb : 1));
    bool*   bnd_lo  = (bool*)malloc(sizeof(bool) * (nb ? nb : 1));
    double* base   = (double*)malloc(sizeof(double) * n);
    double* gbase  = (double*)malloc(sizeof(double) * ngv);
    double* xs     = (double*)malloc(sizeof(double) * n);   /* cross sample point */
    double* gs     = (double*)malloc(sizeof(double) * ngv); /* plus-side gens vals */
    double* gs2    = (double*)malloc(sizeof(double) * ngv); /* minus-side gens vals */
    double* af     = (double*)malloc(sizeof(double) * n);
    double* agen   = (double*)malloc(sizeof(double) * ngv * n); /* row k = grad gens[k] */
    double* Amod   = (double*)malloc(sizeof(double) * mc1 * n);
    double* cval0  = (double*)malloc(sizeof(double) * mc1);
    double* dstep  = (double*)malloc(sizeof(double) * n);
    double* xtrial = (double*)malloc(sizeof(double) * n);
    double* gtrial = (double*)malloc(sizeof(double) * ngv);
    double* xbest  = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    double fbase = 0.0;
    bool have_best = false; double best_f = 0.0;
    if (!gi_k || !gi_sgn || !bnd_var || !bnd_lo || !base || !gbase || !xs || !gs
        || !gs2 || !af || !agen || !Amod || !cval0 || !dstep || !xtrial || !gtrial || !xbest)
        goto cleanup;

    { size_t gi = 0;
      for (size_t k = 0; k < ngens; k++) {
          gi_k[gi] = k; gi_sgn[gi] = 1.0; gi++;
          if (gens[k].equality) { gi_k[gi] = k; gi_sgn[gi] = -1.0; gi++; }
      } }
    { size_t bi = 0;
      for (size_t i = 0; i < n; i++) {
          if (boxes && boxes[i].has_lo) { bnd_var[bi] = i; bnd_lo[bi] = true;  bi++; }
          if (boxes && boxes[i].has_hi) { bnd_var[bi] = i; bnd_lo[bi] = false; bi++; }
      } }

    if (boxes) fm_project_box(x, n, boxes);
    for (size_t i = 0; i < n; i++) base[i] = x[i];
    if (!fm_cobyla_eval(f, gens, ngens, binds, opts, base, n, &fbase, gbase)) {
        fm_warn(g_fm_name, "nlnum", "objective/constraint evaluation failed at start point");
        goto cleanup;
    }

    /* rho schedule: rhobeg from the box scale (default 1), rhoend = tol_prec. */
    double rhoend = pow(10.0, -opts->prec_goal_digits);
    double rhobeg = 1.0;
    for (size_t i = 0; i < n; i++)
        if (boxes && boxes[i].has_lo && boxes[i].has_hi) {
            double half = 0.5 * (boxes[i].hi - boxes[i].lo);
            if (half > 0.0 && half < rhobeg) rhobeg = half;
        }
    if (rhobeg < 4.0 * rhoend) rhobeg = 4.0 * rhoend;
    double rho = rhobeg;
    double parmu = 0.0;

    for (int64_t it = 0; it < opts->max_iter; it++) {
        /* Linear models by CENTRAL differences on a coordinate cross of side rho.
         * Central (vs forward) differences are essential for the non-smooth
         * objectives COBYLA targets: at a kink the two-sided slope averages the
         * left and right derivatives (e.g. |t| at t=0 reads 0, correctly flat),
         * where a one-sided difference would report a spurious ±1 and strand the
         * search on the kink. Each side is clamped into the box so every sample
         * stays feasible; a coordinate pinned against a bound falls back to the
         * one available side. */
        for (size_t i = 0; i < n; i++) af[i] = 0.0;
        bool model_ok = true;
        for (size_t j = 0; j < n && model_ok; j++) {
            double sp = rho, sm = rho;   /* plus / minus half-widths */
            if (boxes) {
                if (boxes[j].has_hi) { double up = boxes[j].hi - base[j]; if (up < sp) sp = up; }
                if (boxes[j].has_lo) { double dn = base[j] - boxes[j].lo; if (dn < sm) sm = dn; }
            }
            if (sp < 1e-12 * (1.0 + fabs(base[j]))) sp = 0.0;
            if (sm < 1e-12 * (1.0 + fabs(base[j]))) sm = 0.0;
            if (sp == 0.0 && sm == 0.0) { sp = 1e-8; }   /* degenerate box: nudge */
            double fp = fbase, fm = fbase;
            if (sp > 0.0) {
                for (size_t i = 0; i < n; i++) xs[i] = base[i];
                xs[j] = base[j] + sp;
                if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xs, n, &fp, gs)) { model_ok = false; break; }
            } else { for (size_t k = 0; k < ngens; k++) gs[k] = gbase[k]; }
            if (sm > 0.0) {
                for (size_t i = 0; i < n; i++) xs[i] = base[i];
                xs[j] = base[j] - sm;
                if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xs, n, &fm, gs2)) { model_ok = false; break; }
            } else { for (size_t k = 0; k < ngens; k++) gs2[k] = gbase[k]; }
            double denom = sp + sm;                 /* > 0 by construction */
            af[j] = (fp - fm) / denom;
            for (size_t k = 0; k < ngens; k++) agen[k * n + j] = (gs[k] - gs2[k]) / denom;
        }
        if (!model_ok) {
            /* A non-numeric sample: shrink rho and retry, or stop. */
            if (rho <= rhoend) { ok = true; break; }
            rho *= 0.5; if (rho <= 1.5 * rhoend) rho = rhoend;
            continue;
        }

        /* Assemble the internal constraint models at the base. */
        for (size_t gi = 0; gi < ng; gi++) {
            double sgn = gi_sgn[gi]; size_t k = gi_k[gi];
            double* row = &Amod[gi * n];
            for (size_t i = 0; i < n; i++) row[i] = sgn * agen[k * n + i];
            cval0[gi] = sgn * gbase[k];
        }
        for (size_t bi = 0; bi < nb; bi++) {
            size_t i = bnd_var[bi];
            double* row = &Amod[(ng + bi) * n];
            for (size_t t = 0; t < n; t++) row[t] = 0.0;
            if (bnd_lo[bi]) { row[i] = -1.0; cval0[ng + bi] = boxes[i].lo - base[i]; }
            else            { row[i] =  1.0; cval0[ng + bi] = base[i] - boxes[i].hi; }
        }

        if (!fm_cobyla_subproblem(n, mcon, Amod, cval0, af, rho, dstep)) goto cleanup;
        double dnorm = 0.0;
        for (size_t i = 0; i < n; i++) if (fabs(dstep[i]) > dnorm) dnorm = fabs(dstep[i]);

        /* Short step at this rho -> reduce rho (or converge). */
        if (dnorm < 0.5 * rho) {
            if (rho <= rhoend) { ok = true; break; }
            rho *= 0.5; if (rho <= 1.5 * rhoend) rho = rhoend;
            continue;
        }

        for (size_t i = 0; i < n; i++) xtrial[i] = base[i] + dstep[i];
        if (boxes) fm_project_box(xtrial, n, boxes);
        double ftrial;
        if (!fm_cobyla_eval(f, gens, ngens, binds, opts, xtrial, n, &ftrial, gtrial))
            goto cleanup;
        fm_fire_monitor(opts->step_monitor);

        /* PARMU update (Powell): raise mu so the predicted merit reduction stays
         * positive, based on the objective cost per unit predicted violation
         * reduction. */
        double base_viol = fm_cobyla_maxviol(gbase, base, n, gi_k, gi_sgn, ng,
                                             bnd_var, bnd_lo, nb, boxes);
        double lin_viol = 0.0;
        for (size_t j = 0; j < mcon; j++) {
            double v = cval0[j];
            for (size_t i = 0; i < n; i++) v += Amod[j * n + i] * dstep[i];
            if (v > lin_viol) lin_viol = v;
        }
        if (lin_viol < 0.0) lin_viol = 0.0;
        double prerec = base_viol - lin_viol;
        double afd = 0.0;
        for (size_t i = 0; i < n; i++) afd += af[i] * dstep[i];
        if (prerec > 0.0) {
            double barmu = afd / prerec;
            if (parmu < 1.5 * barmu) parmu = 2.0 * barmu;
            if (parmu > 1e12) parmu = 1e12;
        }

        double trial_viol = fm_cobyla_maxviol(gtrial, xtrial, n, gi_k, gi_sgn, ng,
                                              bnd_var, bnd_lo, nb, boxes);
        double trial_merit = ftrial + parmu * trial_viol;
        double base_merit  = fbase  + parmu * base_viol;

        if (trial_viol <= 1e-8 && (!have_best || ftrial < best_f)) {
            for (size_t i = 0; i < n; i++) xbest[i] = xtrial[i];
            best_f = ftrial; have_best = true;
        }

        if (trial_merit < base_merit - 1e-12 * (1.0 + fabs(base_merit))) {
            /* Accept: move the base to the trial point. */
            for (size_t i = 0; i < n; i++) base[i] = xtrial[i];
            fbase = ftrial;
            for (size_t k = 0; k < ngens; k++) gbase[k] = gtrial[k];
        } else {
            /* No improvement at this rho -> shrink (or converge). */
            if (rho <= rhoend) { ok = true; break; }
            rho *= 0.5; if (rho <= 1.5 * rhoend) rho = rhoend;
        }
    }

    /* Report the base (best-by-merit) iterate, or the best strictly-feasible
     * point if the base ended infeasible. */
    {
        double base_viol = fm_cobyla_maxviol(gbase, base, n, gi_k, gi_sgn, ng,
                                             bnd_var, bnd_lo, nb, boxes);
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
    free(gi_k); free(gi_sgn); free(bnd_var); free(bnd_lo); free(base); free(gbase);
    free(xs); free(gs); free(gs2); free(af); free(agen); free(Amod); free(cval0);
    free(dstep); free(xtrial); free(gtrial); free(xbest);
    return ok;
}
