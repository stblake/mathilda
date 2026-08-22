/* findmin_penalty.c — penalty outer loop for general-constraint local solves.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Penalty outer loop                                                  *
 * ------------------------------------------------------------------ */

bool fm_run_penalty(Expr* f, Expr** vars, size_t n,
                           FmVarBind* binds, FmMethod method,
                           Expr** g_exprs, Expr*** H_exprs,
                           double* x, /* in/out */
                           const FmGenCon* gens, size_t ngens,
                           const FmBox* boxes,
                           const FmOpts* opts,
                           double* fx_out) {
    /* Outer μ schedule: 1 → 10 → ... up to 10^8.
     *
     * The classical penalty-method convergence theorem requires μ → ∞:
     * for problems with an active constraint at the optimum, the unaugmented
     * gradient is non-zero on the constraint surface, so the augmented
     * minimizer sits a 1/μ-sized step inside the infeasible region. Stopping
     * as soon as the iterate becomes feasible (the previous behaviour) can
     * therefore terminate well before the constrained optimum — e.g. with
     *
     *   FindMinimum[{x + y, 3 x + 2 y >= 7 && x >= 0 && y >= 0}, {x, y}]
     *
     * the BFGS inner solver lands at (x ≈ 2.45, 0) once it becomes feasible,
     * while the true minimum is (7/3, 0). Continuing to ramp μ pulls the
     * iterate down the active boundary to the corner.
     *
     * Termination: stop once both (a) the iterate is feasible and (b)
     * inter-round movement is small (i.e. the augmented minimum has
     * stabilised), or the schedule is exhausted. */
    double mu = 1.0;
    double fx = 0.0;
    bool feas = false;
    double* x_prev = (double*)malloc(sizeof(double) * n);
    if (!x_prev) return false;
    /* Best FEASIBLE iterate seen across the whole μ schedule. The high-μ rounds
     * are ill-conditioned — the augmented objective f + μ·Σpenalty is dominated
     * by the penalty term — so for a hard constraint (e.g. a bilinear equality
     * with large-magnitude terms) the final round can DRIFT to a point *more*
     * infeasible than an earlier round already reached, making the polish report
     * the infeasible sentinel on a problem that is in fact feasible (and making
     * a larger MaxIterations, which converges each ill-conditioned round harder,
     * return a worse answer). Remember the lowest-objective feasible point and
     * fall back to it ONLY when the last round ends infeasible, so any run that
     * already ends feasible is bit-for-bit unchanged. */
    double* best_feas_x = (double*)malloc(sizeof(double) * n);
    if (!best_feas_x) { free(x_prev); return false; }
    double best_feas_f = 0.0;
    bool have_feas = false;
    const double feas_eps = 1.0e-8;   /* matches NMinimize's selection tolerance */
    double pen = 0.0;
    const double rel_tol = pow(10.0, -opts->prec_goal_digits);
    /* Augmented-Lagrangian multipliers, one per general constraint, all zero at
     * the first round (so round 0 is exactly the classical quadratic penalty).
     * Installed on the file statics that fm_eval_augmented / fm_eval_aug_gradient
     * consult, and updated by the PHR rule after each inner solve. Saved and
     * restored so a nested solver call is unaffected. */
    double* al_lambda = NULL;
    if (ngens > 0) {
        al_lambda = (double*)calloc(ngens, sizeof(double));
        if (!al_lambda) { free(x_prev); free(best_feas_x); return false; }
    }
    const double* al_saved_lambda = g_fm_al_lambda;
    const FmGenCon* al_saved_gens = g_fm_al_gens;
    g_fm_al_lambda = al_lambda;      /* NULL when ngens==0 → AL branch inert */
    g_fm_al_gens   = gens;
    const double LAM_MAX = 1.0e12;
    for (int round = 0; round < 9; round++) {
        for (size_t i = 0; i < n; i++) x_prev[i] = x[i];
        bool ok;
        switch (method) {
            case FM_METHOD_QUASINEWTON:
                ok = fm_run_bfgs(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_CONJGRAD:
                ok = fm_run_cg(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_NEWTON:
                ok = fm_run_newton(f, vars, n, binds, g_exprs, H_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_LBFGSB:
                ok = fm_run_lbfgsb(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_TNC:
                ok = fm_run_tnc(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            default:
                ok = fm_run_bfgs(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
        }
        if (!ok) {
            g_fm_al_lambda = al_saved_lambda; g_fm_al_gens = al_saved_gens;
            free(al_lambda); free(x_prev); free(best_feas_x); return false;
        }
        if (!fm_eval_penalty(gens, ngens, binds, x, n, opts, &pen)) {
            g_fm_al_lambda = al_saved_lambda; g_fm_al_gens = al_saved_gens;
            free(al_lambda); free(x_prev); free(best_feas_x); return false;
        }
        if (pen <= feas_eps) {
            double f_here;
            if (fm_eval_scalar(f, binds, x, n, opts, &f_here) &&
                (!have_feas || f_here < best_feas_f)) {
                for (size_t i = 0; i < n; i++) best_feas_x[i] = x[i];
                best_feas_f = f_here; have_feas = true;
            }
        }
        /* PHR multiplier update at the post-solve iterate, using this round's μ:
         *   λ_k ← λ_k + 2μ·h_k          (equality)
         *   λ_k ← max(0, λ_k + 2μ·g_k)  (inequality g_k ≤ 0)
         * clamped so a pathological constraint cannot blow a multiplier up. This
         * is what lets the next round reach feasibility at a moderate μ rather
         * than depending on μ→∞. Skipped when there are no general constraints. */
        if (al_lambda) {
            for (size_t k = 0; k < ngens; k++) {
                double c;
                if (!fm_eval_scalar(gens[k].expr, binds, x, n, opts, &c)) continue;
                double nl = al_lambda[k] + 2.0 * mu * c;
                if (!gens[k].equality && nl < 0.0) nl = 0.0;
                if (nl >  LAM_MAX) nl =  LAM_MAX;
                if (nl < -LAM_MAX) nl = -LAM_MAX;
                al_lambda[k] = nl;
            }
        }
        if (pen < 1e-12) feas = true;
        if (feas) {
            /* Stop only if the iterate has stabilised between rounds. */
            double max_step = 0.0, max_x = 0.0;
            for (size_t i = 0; i < n; i++) {
                double ds = fabs(x[i] - x_prev[i]);
                if (ds > max_step) max_step = ds;
                double xa = fabs(x[i]);
                if (xa > max_x) max_x = xa;
            }
            if (max_step < rel_tol * (max_x + 1.0)) break;
        }
        mu *= 10.0;
    }
    /* Drift rescue: the final iterate is infeasible but an earlier round found a
     * feasible point — return that instead. Never triggers when the run already
     * ends feasible (pen ≤ feas_eps), so existing feasible results are unchanged. */
    if (pen > feas_eps && have_feas) {
        for (size_t i = 0; i < n; i++) x[i] = best_feas_x[i];
    }
    g_fm_al_lambda = al_saved_lambda;
    g_fm_al_gens   = al_saved_gens;
    free(al_lambda);
    free(x_prev);
    free(best_feas_x);
    if (!feas) {
        fm_warn(g_fm_name, "infeas", "could not satisfy constraints to tolerance");
    }
    /* Report unaugmented objective value at the final point. */
    if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) return false;
    *fx_out = fx;
    return true;
}
