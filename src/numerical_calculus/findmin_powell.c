/* findmin_powell.c — Powell conjugate-direction derivative-free local solver.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Powell's conjugate-direction method (derivative-free)               *
 * ------------------------------------------------------------------ *
 * Minimises f by sweeping 1-D line minimisations along an evolving set
 * of n directions (initialised to the unit basis), replacing the
 * direction of largest decrease with the averaged cycle step when
 * Powell's parabolic test accepts it.  Function values only -- no
 * gradient -- so it serves the non-smooth / black-box objectives the
 * gradient methods cannot.  The line search reuses the same Brent
 * parabolic-interpolation minimiser as the 1-D "Brent" method,
 * restricted to phi(t) = f(p + t d) via fm_eval_line.  Machine precision
 * only: it rides the compiled objective in g_fm_obj_prog through
 * fm_eval_scalar, and WorkingPrecision > MachinePrecision falls back to
 * QuasiNewton in the driver.  The replace-direction test, the direction
 * cycling and the extrapolated point all match scipy's _minimize_powell,
 * so the two agree on the minimiser to rounding.  Exposed as
 * Method -> "Powell" and its Mathematica alias "PrincipalAxis". */

/* phi(t) = f(p + t d), written into xtmp.  No box projection here:
 * clamping t (fm_bracket_line / fm_brent_line) is what keeps the line
 * feasible; projecting inside the eval would make phi non-smooth and
 * defeat Brent's parabolic model. */
static bool fm_eval_line(Expr* f, FmVarBind* binds, size_t n,
                         const double* p, const double* d, double t,
                         const FmOpts* opts, double* xtmp, double* fval_out) {
    for (size_t i = 0; i < n; i++) xtmp[i] = p[i] + t * d[i];
    return fm_eval_scalar(f, binds, xtmp, n, opts, fval_out);
}

/* mnbrak-style bracketing of phi(t) = f(p + t d) over t in [t_lo, t_hi].
 * Mirrors fm_bracket, but the step scale lives in d (unit-norm basis
 * rows, or a replaced row carrying its own magnitude), so it starts at
 * t = 0 with unit step h = 1 -- scipy's mnbrak default (xa=0, xb=1) --
 * rather than the coordinate-scaled 1e-2 of the 1-D fm_bracket.  When
 * unbounded the interval is [-HUGE_VAL, +HUGE_VAL] and every clamp is a
 * no-op. */
static bool fm_bracket_line(Expr* f, FmVarBind* binds, size_t n,
                            const double* p, const double* d, const FmOpts* opts,
                            double t_lo, double t_hi, double* xtmp,
                            double* a_out, double* b_out, double* c_out) {
    double a = 0.0, b, fa, fb;
    if (!fm_eval_line(f, binds, n, p, d, a, opts, xtmp, &fa)) return false;
    b = a + 1.0;
    if (b > t_hi) b = 0.5 * (a + t_hi);
    if (b < t_lo) b = 0.5 * (a + t_lo);
    if (!fm_eval_line(f, binds, n, p, d, b, opts, xtmp, &fb)) return false;
    if (fb > fa) {
        double t = a; a = b; b = t;
        t = fa; fa = fb; fb = t;
    }
    double c = b + 1.618 * (b - a);
    if (c > t_hi) c = t_hi;
    if (c < t_lo) c = t_lo;
    double fc;
    if (!fm_eval_line(f, binds, n, p, d, c, opts, xtmp, &fc)) return false;
    /* Strict `<` (scipy's mnbrak), NOT the 1-D fm_bracket's `<=`: on a flat
     * coordinate direction (e.g. Beale is constant in x at y==1) `<=` keeps
     * growing the bracket to where floating-point cancellation fakes a
     * spurious minimum at ~1e37; strict `<` stops at the flat region. */
    for (int k = 0; k < 100 && fc < fb; k++) {
        a = b; fa = fb;
        b = c; fb = fc;
        c = b + 1.618 * (b - a);
        if (c >= t_hi) {
            c = t_hi;
            if (!fm_eval_line(f, binds, n, p, d, c, opts, xtmp, &fc)) return false;
            break;
        }
        if (c <= t_lo) {
            c = t_lo;
            if (!fm_eval_line(f, binds, n, p, d, c, opts, xtmp, &fc)) return false;
            break;
        }
        if (!fm_eval_line(f, binds, n, p, d, c, opts, xtmp, &fc)) return false;
    }
    if (a > c) { double t = a; a = c; c = t; }
    *a_out = a; *b_out = b; *c_out = c;
    return true;
}

/* Brent parabolic-interpolation minimiser of phi(t) = f(p + t d) on the
 * bracket [a, c], clamped to [t_lo, t_hi].  One-for-one transliteration
 * of fm_brent_min with abscissa t and line evaluation, omitting the
 * per-step StepMonitor (Powell fires it once per cycle instead) and the
 * non-convergence warning (the outer loop's function-decrease test is the
 * authoritative signal). */
static bool fm_brent_line(Expr* f, FmVarBind* binds, size_t n,
                          const double* p, const double* d, const FmOpts* opts,
                          double a, double b, double c, double t_lo, double t_hi,
                          double* xtmp, double* t_out, double* fx_out) {
    if (a > c) { double t = a; a = c; c = t; }
    double tol = pow(10.0, -opts->prec_goal_digits);
    double tol_acc = pow(10.0, -opts->acc_goal_digits);
    double e = 0.0, dd = 0.0;
    double x, w, v;
    x = w = v = b;
    if (x < a || x > c) x = w = v = 0.5 * (a + c);
    double fx;
    if (!fm_eval_line(f, binds, n, p, d, x, opts, xtmp, &fx)) return false;
    double fw = fx, fv = fx;
    for (int64_t k = 0; k < opts->max_iter; k++) {
        double xm = 0.5 * (a + c);
        double tol1 = tol * fabs(x) + FM_ZEPS;
        double tol2 = 2.0 * tol1;
        if (fabs(x - xm) <= tol2 - 0.5 * (c - a)
            || fabs(fx) < tol_acc * (1.0 + fabs(fx))) {
            *t_out = x; *fx_out = fx; return true;
        }
        double u = x;
        if (fabs(e) > tol1) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double pp = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) pp = -pp;
            q = fabs(q);
            double etemp = e;
            e = dd;
            if (fabs(pp) >= fabs(0.5 * q * etemp)
                || pp <= q * (a - x) || pp >= q * (c - x)) {
                e = (x >= xm) ? (a - x) : (c - x);
                dd = FM_CGOLD * e;
            } else {
                dd = pp / q;
                u = x + dd;
                if (u - a < tol2 || c - u < tol2)
                    dd = (xm - x >= 0.0) ? tol1 : -tol1;
            }
        } else {
            e = (x >= xm) ? (a - x) : (c - x);
            dd = FM_CGOLD * e;
        }
        u = (fabs(dd) >= tol1) ? (x + dd) : (x + ((dd >= 0.0) ? tol1 : -tol1));
        if (u < t_lo) u = t_lo;
        if (u > t_hi) u = t_hi;
        double fu;
        if (!fm_eval_line(f, binds, n, p, d, u, opts, xtmp, &fu)) return false;
        if (fu <= fx) {
            if (u >= x) a = x; else c = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else c = u;
            if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    *t_out = x; *fx_out = fx;
    return true;
}

/* Feasible t-interval for p + t d against the per-variable box, then
 * bracket + Brent-minimise phi(t).  On a strict improvement, advance p by
 * t*.d (re-projected), update *f_cur, and report the realised step t*.
 * Returns false (leaving p and *f_cur untouched, *t_star = 0) when the
 * direction is degenerate, the feasible interval collapses, or an
 * evaluation fails -- the caller simply skips that direction. */
static bool fm_powell_line_min(Expr* f, FmVarBind* binds, size_t n,
                               double* p, const double* d,
                               const FmBox* boxes, const FmOpts* opts,
                               double* xtmp, double* f_cur, double* t_star) {
    *t_star = 0.0;
    double dn = 0.0;
    for (size_t i = 0; i < n; i++) dn += d[i] * d[i];
    if (dn < 1e-60) return false;                 /* ||d|| < 1e-30 */

    double t_lo = -HUGE_VAL, t_hi = HUGE_VAL;
    if (boxes) {
        for (size_t i = 0; i < n; i++) {
            if (fabs(d[i]) < 1e-30) continue;     /* component fixed */
            double tl, tu;
            if (d[i] > 0.0) {
                tl = boxes[i].has_lo ? (boxes[i].lo - p[i]) / d[i] : -HUGE_VAL;
                tu = boxes[i].has_hi ? (boxes[i].hi - p[i]) / d[i] :  HUGE_VAL;
            } else {                              /* sign flip */
                tl = boxes[i].has_hi ? (boxes[i].hi - p[i]) / d[i] : -HUGE_VAL;
                tu = boxes[i].has_lo ? (boxes[i].lo - p[i]) / d[i] :  HUGE_VAL;
            }
            if (tl > t_lo) t_lo = tl;
            if (tu < t_hi) t_hi = tu;
        }
        if (t_hi - t_lo < 1e-15) return false;    /* pinned in a corner along d */
    }

    double ta, tb, tc;
    if (!fm_bracket_line(f, binds, n, p, d, opts, t_lo, t_hi, xtmp, &ta, &tb, &tc))
        return false;
    double tmin, fmin;
    if (!fm_brent_line(f, binds, n, p, d, opts, ta, tb, tc, t_lo, t_hi, xtmp, &tmin, &fmin))
        return false;
    if (fmin < *f_cur) {
        for (size_t i = 0; i < n; i++) p[i] += tmin * d[i];
        if (boxes) fm_project_box(p, n, boxes);
        *f_cur = fmin;
        *t_star = tmin;
        return true;
    }
    return false;
}

bool fm_run_powell(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars; (void)g_exprs; (void)gens; (void)ngens; (void)mu;

    /* n == 1 degenerates to a single Brent line search along e_0; delegate
     * to the exact 1-D path (identical result, cheaper) rather than run the
     * full direction-set machinery over one coordinate. */
    if (n == 1) {
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

    double* direc   = (double*)calloc(n * n, sizeof(double));
    double* x_start = (double*)malloc(sizeof(double) * n);
    double* x_extra = (double*)malloc(sizeof(double) * n);
    double* d_avg   = (double*)malloc(sizeof(double) * n);
    double* xtmp    = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    if (!direc || !x_start || !x_extra || !d_avg || !xtmp) goto cleanup;

    for (size_t i = 0; i < n; i++) direc[i*n + i] = 1.0;
    if (boxes) fm_project_box(x, n, boxes);

    double f_cur;
    if (!fm_eval_scalar(f, binds, x, n, opts, &f_cur)) {
        fm_warn(g_fm_name, "nlnum", "objective evaluation failed at start point");
        goto cleanup;
    }

    double ftol     = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        for (size_t i = 0; i < n; i++) x_start[i] = x[i];
        double f_start = f_cur;
        double delta = 0.0;
        size_t ibig = 0;

        /* One sweep of n line minimisations over the current direction set. */
        for (size_t i = 0; i < n; i++) {
            double f_before = f_cur;
            double tstep;
            fm_powell_line_min(f, binds, n, x, &direc[i*n], boxes, opts,
                               xtmp, &f_cur, &tstep);
            double dec = f_before - f_cur;
            if (dec > delta) { delta = dec; ibig = i; }
        }

        /* Convergence: relative function decrease over a full cycle (scipy's
         * Powell test), or a PrecisionGoal step-size floor. */
        if (2.0 * (f_start - f_cur) <= ftol * (fabs(f_start) + fabs(f_cur)) + 1e-20) {
            ok = true; break;
        }
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x[i] - x_start[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x[i]) > max_x) max_x = fabs(x[i]);
        }
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }

        /* Averaged direction d_avg = x - x_start and extrapolated point
         * x_e = 2x - x_start (== x + d_avg). */
        for (size_t i = 0; i < n; i++) {
            d_avg[i]   = x[i] - x_start[i];
            x_extra[i] = 2.0 * x[i] - x_start[i];
        }
        if (boxes) fm_project_box(x_extra, n, boxes);
        double f_extra;
        if (!fm_eval_scalar(f, binds, x_extra, n, opts, &f_extra)) {
            /* Non-finite extrapolate -> keep the direction set, continue. */
            fm_fire_monitor(opts->step_monitor);
            continue;
        }

        /* Powell's replace-direction test (Numerical Recipes / scipy form). */
        if (f_start > f_extra) {
            double s1 = f_start - f_cur - delta;
            double s2 = f_start - f_extra;
            double t = 2.0 * (f_start + f_extra - 2.0 * f_cur) * s1 * s1
                     - delta * s2 * s2;
            if (t < 0.0) {
                double tstep;
                bool moved = fm_powell_line_min(f, binds, n, x, d_avg, boxes,
                                                opts, xtmp, &f_cur, &tstep);
                /* Linear-dependence guard: rotate in only a non-zero step. */
                if (moved && tstep != 0.0) {
                    for (size_t i = 0; i < n; i++) {
                        direc[ibig*n + i]  = direc[(n-1)*n + i];
                        direc[(n-1)*n + i] = tstep * d_avg[i];
                    }
                }
            }
        }
        fm_fire_monitor(opts->step_monitor);
    }

    *fx_out = f_cur;
    ok = true;
cleanup:
    free(direc); free(x_start); free(x_extra); free(d_avg); free(xtmp);
    return ok;
}
