/*
 * gkadapt.c — globally-adaptive Gauss-Kronrod (QK15) quadrature.  See gkadapt.h.
 *
 * QK15 nodes/weights are the standard QUADPACK constants.  The adaptive driver
 * keeps the subintervals in a binary max-heap keyed by local error so the worst
 * panel is always bisected first; this is the QAG strategy.  When extrapolation
 * is requested the running total after each bisection feeds Wynn's epsilon
 * (shared seqaccel kernel), giving the QAGS acceleration for endpoint
 * singularities.
 */

#include "gkadapt.h"
#include "seqaccel.h"
#include "denint.h"   /* tanh-sinh: the endpoint singularity handler */

#include <float.h>
#include <math.h>
#include <stdlib.h>

/* ---- QK15: 7-point Gauss embedded in a 15-point Kronrod rule ---- */
/* Abscissae of the 15-point Kronrod rule (nonnegative half; xgk[7] = centre). */
static const double XGK[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000
};
/* Weights of the 15-point Kronrod rule.  Full QUADPACK precision: at 15
 * significant digits they summed to 2 - 6e-15, a 3e-15 relative error in
 * every integral (the last digit of a smooth machine integrand). */
static const double WGK[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714
};
/* Weights of the 7-point Gauss rule (wg[3] is the centre weight). */
static const double WG[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327
};

/* One QK15 panel on [a,b].  Returns false if any sample is non-numeric. */
static bool qk15(GkSampleMachine f, void* ctx, double a, double b,
                 double _Complex* result, double* abserr,
                 double* resabs_out, double* resasc_out, long* nev) {
    const double epmach = DBL_EPSILON;
    const double uflow  = DBL_MIN;
    double centr = 0.5 * (a + b);
    double hlgth = 0.5 * (b - a);
    double dhlgth = fabs(hlgth);

    double _Complex fv1[8], fv2[8], fc;
    if (!f(ctx, centr, &fc)) return false;
    (*nev)++;

    double _Complex resg = WG[3] * fc;     /* 7-pt Gauss   */
    double _Complex resk = WGK[7] * fc;    /* 15-pt Kronrod */
    double resabs = cabs(resk);

    /* Gauss/Kronrod shared abscissae (odd xgk indices 1,3,5). */
    for (int j = 0; j < 3; j++) {
        int idx = 2 * j + 1;
        double absc = hlgth * XGK[idx];
        double _Complex f1, f2;
        if (!f(ctx, centr - absc, &f1)) return false;
        if (!f(ctx, centr + absc, &f2)) return false;
        *nev += 2;
        fv1[idx] = f1; fv2[idx] = f2;
        double _Complex fsum = f1 + f2;
        resg += WG[j] * fsum;
        resk += WGK[idx] * fsum;
        resabs += WGK[idx] * (cabs(f1) + cabs(f2));
    }
    /* Kronrod-only abscissae (even xgk indices 0,2,4,6). */
    for (int j = 0; j < 4; j++) {
        int idx = 2 * j;
        double absc = hlgth * XGK[idx];
        double _Complex f1, f2;
        if (!f(ctx, centr - absc, &f1)) return false;
        if (!f(ctx, centr + absc, &f2)) return false;
        *nev += 2;
        fv1[idx] = f1; fv2[idx] = f2;
        double _Complex fsum = f1 + f2;
        resk += WGK[idx] * fsum;
        resabs += WGK[idx] * (cabs(f1) + cabs(f2));
    }

    double _Complex reskh = resk * 0.5;
    double resasc = WGK[7] * cabs(fc - reskh);
    for (int j = 0; j < 7; j++)
        resasc += WGK[j] * (cabs(fv1[j] - reskh) + cabs(fv2[j] - reskh));

    *result = resk * hlgth;
    resabs *= dhlgth;
    resasc *= dhlgth;
    double err = cabs((resk - resg) * hlgth);
    if (resasc != 0.0 && err != 0.0)
        err = resasc * fmin(1.0, pow(200.0 * err / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) {
        double floor_err = (epmach * 50.0) * resabs;
        if (floor_err > err) err = floor_err;
    }
    *abserr = err;
    if (resabs_out) *resabs_out = resabs;
    if (resasc_out) *resasc_out = resasc;
    return true;
}

/* ---- subinterval heap (max-heap on err) ---- */
typedef struct {
    double a, b;
    double _Complex result;
    double err;
    double perr;   /* error estimate of the parent panel (INFINITY at the root) */
} Panel;

typedef struct { Panel* v; int n, cap; } Heap;

static void heap_swap(Heap* h, int i, int j) { Panel t = h->v[i]; h->v[i] = h->v[j]; h->v[j] = t; }

static void heap_up(Heap* h, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->v[p].err >= h->v[i].err) break;
        heap_swap(h, p, i); i = p;
    }
}
static void heap_down(Heap* h, int i) {
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, m = i;
        if (l < h->n && h->v[l].err > h->v[m].err) m = l;
        if (r < h->n && h->v[r].err > h->v[m].err) m = r;
        if (m == i) break;
        heap_swap(h, m, i); i = m;
    }
}
static bool heap_push(Heap* h, Panel p) {
    if (h->n == h->cap) {
        int nc = h->cap ? h->cap * 2 : 32;
        Panel* nv = realloc(h->v, (size_t)nc * sizeof(Panel));
        if (!nv) return false;
        h->v = nv; h->cap = nc;
    }
    h->v[h->n] = p; heap_up(h, h->n); h->n++;
    return true;
}
static Panel heap_pop(Heap* h) { Panel top = h->v[0]; h->n--; h->v[0] = h->v[h->n]; heap_down(h, 0); return top; }

/* ---- endpoint singularity handler ---- */

/* Counting thunk: the tanh-sinh panel rule samples through this so its
 * evaluations are charged to the run's n_eval / max_eval budget. */
typedef struct { GkSampleMachine f; void* ctx; long n; } GkCount;
static bool gk_count_sample(void* v, double x, double _Complex* out) {
    GkCount* c = (GkCount*)v;
    c->n++;
    return c->f(c->ctx, x, out);
}

/* Tanh-sinh levels for one endpoint panel: level L costs ~9 * 2^L samples and
 * the levels are not nested, so 6 bounds a panel at ~1100 samples.  An
 * algebraic or logarithmic endpoint singularity converges by level 3-5. */
#define GK_SING_DE_LEVELS 6

/* Integrate an endpoint panel [pa,pb] by tanh-sinh, whose double-exponential
 * clustering resolves an integrable endpoint singularity that bisection only
 * approaches algebraically (x^-1/2, (1-x)^1/4, Log[x], ...).  The rule is
 * asked for the run's relative tolerance on the panel; the level difference it
 * converges on is the error of the level BEFORE the last, so the returned
 * value is accurate to roughly the square of that tolerance.  Accepted only
 * when it converged and agrees with the Kronrod value to within the panel's
 * own error estimate (+ the new one): a panel the rule cannot resolve, or a
 * near-endpoint feature it misjudges, is left to ordinary bisection. */
static bool gk_sing_panel(GkSampleMachine f, void* ctx, const Panel* p,
                          double reltol, long* nev,
                          double _Complex* val, double* err) {
    GkCount cc = { f, ctx, 0 };
    double _Complex v; double e;
    double rt = reltol > 1e-14 ? reltol : 1e-14;
    bool conv = denint_tanhsinh_machine(gk_count_sample, &cc, p->a, p->b, rt,
                                        GK_SING_DE_LEVELS, &v, &e);
    *nev += cc.n;
    if (!conv || !isfinite(creal(v)) || !isfinite(cimag(v)) || !isfinite(e))
        return false;
    if (cabs(v - p->result) > p->err + e + 50.0 * DBL_EPSILON * cabs(v))
        return false;
    /* roundoff floor, as qk15's */
    double floor_err = 50.0 * DBL_EPSILON * cabs(v);
    *val = v;
    *err = e > floor_err ? e : floor_err;
    return true;
}

GkResult gk_integrate_machine(GkSampleMachine f, void* ctx,
                              double a, double b,
                              double abstol, double reltol,
                              int max_subdiv, long max_eval,
                              bool extrapolate) {
    return gk_integrate_machine_ex(f, ctx, a, b, abstol, reltol, max_subdiv,
                                   max_eval, extrapolate, 0);
}

GkResult gk_integrate_machine_ex(GkSampleMachine f, void* ctx,
                                 double a, double b,
                                 double abstol, double reltol,
                                 int max_subdiv, long max_eval,
                                 bool extrapolate, int sing_depth) {
    GkResult R;
    R.value = 0.0; R.abs_err = 0.0; R.status = GK_NONNUMERIC;
    R.n_eval = 0; R.n_subdiv = 0;

    if (max_subdiv < 0) max_subdiv = 0;
    if (max_eval <= 0) max_eval = (long)(max_subdiv + 1) * 15L + 30L;
    if (abstol < 0.0) abstol = 0.0;
    if (reltol < 0.0) reltol = 0.0;

    Panel p0;
    p0.a = a; p0.b = b; p0.perr = INFINITY;
    if (!qk15(f, ctx, a, b, &p0.result, &p0.err, NULL, NULL, &R.n_eval))
        return R;   /* non-numeric on the whole interval */

    Heap h = { NULL, 0, 0 };
    if (!heap_push(&h, p0)) { R.status = GK_NONNUMERIC; return R; }

    double _Complex total = p0.result;
    double total_err = p0.err;

    /* Endpoint singularity handler (sing_depth > 0): a panel touching a or b
     * that bisection has already narrowed to (b-a) 2^-sing_depth is settled by
     * tanh-sinh (gk_sing_panel) instead of being halved again, and leaves the
     * heap.  One attempt per endpoint: a decline returns that end to plain
     * bisection for the rest of the run. */
    double sing_w = (sing_depth > 0 && sing_depth < 60)
                  ? (b - a) * ldexp(1.0, -sing_depth) : -1.0;
    bool sing_tried_lo = false, sing_tried_hi = false;
    double _Complex settled = 0.0;   /* sum of the settled panels' values */

    /* Extrapolation: feed the running total after each bisection to Wynn. */
    double _Complex* seq = NULL;
    int seq_n = 0, seq_cap = 0;

    R.status = (total_err <= fmax(abstol, reltol * cabs(total))) ? GK_OK : GK_NOCONV;
    while (R.status != GK_OK && R.n_subdiv < max_subdiv && R.n_eval + 30 <= max_eval) {
        double tol = fmax(abstol, reltol * cabs(total));
        if (total_err <= tol) { R.status = GK_OK; break; }

        Panel worst = heap_pop(&h);
        bool at_lo = (worst.a == a), at_hi = (worst.b == b);
        /* Singular signature: still the worst panel after sing_depth halvings,
         * and the last halving bought less than a factor 20 (an analytic panel
         * gains ~2^-20 a halving once resolved; x^alpha only 2^-(1+alpha)).
         * An under-resolved oscillation running into the end (Sin[Exp[Exp[x]]]
         * at the upper limit) passes too, deliberately: tanh-sinh resolves that
         * panel as well, and the run ends at 15 digits where bisecting it to
         * the goal ended at 13, for ~40% more samples. */
        if (sing_w > 0.0 && (worst.b - worst.a) <= sing_w
            && worst.err > 0.05 * worst.perr
            && ((at_lo && !sing_tried_lo) || (at_hi && !sing_tried_hi))) {
            if (at_lo) sing_tried_lo = true; else sing_tried_hi = true;
            double _Complex sv; double se;
            if (gk_sing_panel(f, ctx, &worst, reltol, &R.n_eval, &sv, &se)) {
                /* settled: out of the heap, its value and error kept in the
                 * running totals; the extrapolation sequence restarts, since
                 * the running sum just jumped */
                total += sv - worst.result;
                settled += sv;
                total_err += se - worst.err;
                if (total_err < 0.0) total_err = 0.0;
                seq_n = 0;
                if (h.n == 0) {
                    R.status = (total_err <= fmax(abstol, reltol * cabs(total)))
                             ? GK_OK : GK_NOCONV;
                    break;
                }
                continue;
            }
        }
        double mid = 0.5 * (worst.a + worst.b);
        if (!(mid > worst.a && mid < worst.b)) {
            /* Interval too narrow to bisect (floating-point underflow toward an
             * endpoint singularity): no further refinement is possible.  Stop —
             * the caller's tanh-sinh fallback handles the singular tail. */
            heap_push(&h, worst);
            break;
        }
        Panel left, right;
        left.a = worst.a; left.b = mid;
        right.a = mid; right.b = worst.b;
        left.perr = right.perr = worst.err;
        bool okl = qk15(f, ctx, left.a, left.b, &left.result, &left.err, NULL, NULL, &R.n_eval);
        bool okr = qk15(f, ctx, right.a, right.b, &right.result, &right.err, NULL, NULL, &R.n_eval);
        if (!okl || !okr) {
            /* singular panel: cannot refine here; report best so far */
            heap_push(&h, worst);
            break;
        }
        /* Update running totals incrementally. */
        total += (left.result + right.result) - worst.result;
        total_err += (left.err + right.err) - worst.err;
        if (total_err < 0.0) total_err = 0.0;
        heap_push(&h, left);
        heap_push(&h, right);
        R.n_subdiv++;

        if (extrapolate) {
            if (seq_n == seq_cap) {
                int nc = seq_cap ? seq_cap * 2 : 16;
                if (nc > 64) nc = 64;
                if (seq_n < nc) {
                    double _Complex* ns = realloc(seq, (size_t)nc * sizeof(double _Complex));
                    if (ns) { seq = ns; seq_cap = nc; }
                }
            }
            if (seq_n < seq_cap) seq[seq_n++] = total;
        }
    }

    /* The running total was updated by one add/subtract per bisection, so
     * it carries ~n_subdiv roundings of |total| -- on a cancelling integrand
     * (|f| integrating to ~40x the result) the last 1-2 digits.  Re-sum the
     * final panels with Neumaier's compensated sum: O(panels), and the value
     * is then as accurate as the panel results themselves. */
    if (R.n_subdiv > 0) {
        double sr = creal(settled), si = cimag(settled), cr = 0.0, ci = 0.0;
        for (int i = 0; i < h.n; i++) {
            double xr = creal(h.v[i].result), xi = cimag(h.v[i].result);
            double tr = sr + xr, ti = si + xi;
            cr += (fabs(sr) >= fabs(xr)) ? (sr - tr) + xr : (xr - tr) + sr;
            ci += (fabs(si) >= fabs(xi)) ? (si - ti) + xi : (xi - ti) + si;
            sr = tr; si = ti;
        }
        total = (sr + cr) + (si + ci) * I;
    }
    R.value = total;
    R.abs_err = total_err;

    /* QAGS-style acceleration: if the raw error stalled, try Wynn on the
     * running-total sequence and accept it when it tightens the estimate. */
    if (extrapolate && R.status != GK_OK && seq_n >= 5) {
        double _Complex acc; double step;
        if (seqaccel_wynn_machine(seq, seq_n, 1, &acc, &step)
            && isfinite(creal(acc)) && isfinite(cimag(acc)) && isfinite(step)) {
            double tol = fmax(abstol, reltol * cabs(acc));
            if (step < total_err && cabs(acc - total) < 100.0 * total_err + tol) {
                R.value = acc;
                R.abs_err = fmin(total_err, step);
                if (R.abs_err <= tol) R.status = GK_OK;
            }
        }
    }

    free(seq);
    free(h.v);
    return R;
}

/* ================================================================== *
 *  QAGS (QUADPACK dqagse + dqk21 + dqelg)                             *
 * ================================================================== */

/* QK21 abscissae (nonnegative half; XGK21[10] = centre) and weights; the
 * 10-point Gauss rule uses the odd-indexed abscissae.  QUADPACK constants. */
static const double XGK21[11] = {
    0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
    0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
    0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
    0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
    0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
    0.000000000000000000000000000000000
};
static const double WGK21[11] = {
    0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077208977235800, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821
};
static const double WG10[5] = {
    0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
    0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
    0.295524224714752870173892994651338
};

/* One QK21 panel (dqk21): *resabs = the rule applied to |f|, *resasc = to
 * |f - mean| (the smoothness gauge of the error estimate). */
static bool qk21(GkSampleMachine f, void* ctx, double a, double b,
                 double _Complex* result, double* abserr,
                 double* resabs, double* resasc, long* nev) {
    const double epmach = DBL_EPSILON, uflow = DBL_MIN;
    double centr = 0.5 * (a + b), hlgth = 0.5 * (b - a), dhlgth = fabs(hlgth);
    double _Complex fv1[10], fv2[10], fc;
    if (!f(ctx, centr, &fc)) return false;
    (*nev)++;
    double _Complex resg = 0.0, resk = WGK21[10] * fc;
    double rabs = cabs(resk);
    for (int j = 0; j < 5; j++) {                 /* Gauss and Kronrod nodes */
        int jt = 2 * j + 1;
        double absc = hlgth * XGK21[jt];
        double _Complex f1, f2;
        if (!f(ctx, centr - absc, &f1) || !f(ctx, centr + absc, &f2)) return false;
        *nev += 2;
        fv1[jt] = f1; fv2[jt] = f2;
        resg += WG10[j] * (f1 + f2);
        resk += WGK21[jt] * (f1 + f2);
        rabs += WGK21[jt] * (cabs(f1) + cabs(f2));
    }
    for (int j = 0; j < 5; j++) {                 /* Kronrod-only nodes */
        int jt = 2 * j;
        double absc = hlgth * XGK21[jt];
        double _Complex f1, f2;
        if (!f(ctx, centr - absc, &f1) || !f(ctx, centr + absc, &f2)) return false;
        *nev += 2;
        fv1[jt] = f1; fv2[jt] = f2;
        resk += WGK21[jt] * (f1 + f2);
        rabs += WGK21[jt] * (cabs(f1) + cabs(f2));
    }
    double _Complex reskh = 0.5 * resk;
    double rasc = WGK21[10] * cabs(fc - reskh);
    for (int j = 0; j < 10; j++)
        rasc += WGK21[j] * (cabs(fv1[j] - reskh) + cabs(fv2[j] - reskh));
    *result = resk * hlgth;
    rabs *= dhlgth; rasc *= dhlgth;
    double err = cabs((resk - resg) * hlgth);
    if (rasc != 0.0 && err != 0.0) err = rasc * fmin(1.0, pow(200.0 * err / rasc, 1.5));
    if (rabs > uflow / (50.0 * epmach)) err = fmax(epmach * 50.0 * rabs, err);
    *abserr = err; *resabs = rabs; *resasc = rasc;
    return true;
}

/* dqelg: the epsilon algorithm on epstab[0 .. *n-1] (capacity
 * QAGS_LIMEXP + 2), with QUADPACK's error estimate from the last three
 * extrapolated results (res3la).  *n may be reduced (table truncation). */
#define QAGS_LIMEXP 50
static void qags_qelg(int* n, double _Complex* epstab, double _Complex* result,
                      double* abserr, double _Complex* res3la, int* nres) {
    const double epmach = DBL_EPSILON, oflow = DBL_MAX;
    (*nres)++;
    *abserr = oflow;
    *result = epstab[*n - 1];
    if (*n < 3) { *abserr = fmax(*abserr, 5.0 * epmach * cabs(*result)); return; }
    /* 1-based indices below, as in QUADPACK */
    epstab[*n + 1] = epstab[*n - 1];
    int newelm = (*n - 1) / 2;
    epstab[*n - 1] = oflow;
    int num = *n, k1 = *n;
    for (int i = 1; i <= newelm; i++) {
        int k2 = k1 - 1, k3 = k1 - 2;
        double _Complex res = epstab[k1 + 1];
        double _Complex e0 = epstab[k3 - 1], e1 = epstab[k2 - 1], e2 = res;
        double e1abs = cabs(e1);
        double _Complex delta2 = e2 - e1, delta3 = e1 - e0;
        double err2 = cabs(delta2), err3 = cabs(delta3);
        double tol2 = fmax(cabs(e2), e1abs) * epmach;
        double tol3 = fmax(e1abs, cabs(e0)) * epmach;
        if (err2 <= tol2 && err3 <= tol3) {
            /* e0, e1 and e2 agree to machine accuracy: convergence assumed */
            *result = res;
            *abserr = fmax(err2 + err3, 5.0 * epmach * cabs(*result));
            return;
        }
        double _Complex e3 = epstab[k1 - 1];
        epstab[k1 - 1] = e1;
        double _Complex delta1 = e1 - e3;
        double err1 = cabs(delta1);
        double tol1 = fmax(e1abs, cabs(e3)) * epmach;
        bool irregular = (err1 <= tol1 || err2 <= tol2 || err3 <= tol3);
        double _Complex ss = 0.0;
        if (!irregular) {
            ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
            irregular = !(cabs(ss * e1) > 1e-4);
        }
        if (irregular) { *n = i + i - 1; break; }   /* omit part of the table */
        res = e1 + 1.0 / ss;
        epstab[k1 - 1] = res;
        k1 -= 2;
        double error = err2 + cabs(res - e2) + err3;
        if (error <= *abserr) { *abserr = error; *result = res; }
    }
    /* shift the table */
    if (*n == QAGS_LIMEXP) *n = 2 * (QAGS_LIMEXP / 2) - 1;
    int ib = ((num / 2) * 2 == num) ? 2 : 1;
    for (int i = 1; i <= newelm + 1; i++) { epstab[ib - 1] = epstab[ib + 1]; ib += 2; }
    if (num != *n) {
        int indx = num - *n + 1;
        for (int i = 1; i <= *n; i++) { epstab[i - 1] = epstab[indx - 1]; indx++; }
    }
    if (*nres < 4) {
        res3la[*nres - 1] = *result;
        *abserr = oflow;
    } else {
        *abserr = cabs(*result - res3la[2]) + cabs(*result - res3la[1])
                + cabs(*result - res3la[0]);
        res3la[0] = res3la[1]; res3la[1] = res3la[2]; res3la[2] = *result;
    }
    *abserr = fmax(*abserr, 5.0 * epmach * cabs(*result));
}

/* The QAGS subinterval list, and iord: the interval indices in decreasing
 * order of error (kept fully sorted, which dqpsrt approximates to save work
 * on a long list; limit is small here). */
typedef struct {
    double *alist, *blist, *elist;
    double _Complex* rlist;
    int* iord;
} QagsList;

/* dqpsrt: after the interval at iord[nrmax-1] was bisected (its first half
 * kept at index *maxerr, the second appended at index last-1), restore the
 * ordering and return the interval to bisect next, iord[nrmax-1].  During
 * extrapolation the nrmax-1 entries above are small intervals being skipped;
 * a half whose error exceeds some of them moves above them, and nrmax drops
 * with it, as in QUADPACK. */
static void qags_sort(QagsList* L, int last, int* nrmax, int* maxerr, double* errmax) {
    int mi = *maxerr;
    double e_mi = L->elist[mi];
    while (*nrmax > 1 && e_mi > L->elist[L->iord[*nrmax - 2]]) (*nrmax)--;
    /* remove mi, then insert mi and the new interval in order */
    int cnt = last - 1, w = 0;
    for (int i = 0; i < cnt; i++) if (L->iord[i] != mi) L->iord[w++] = L->iord[i];
    int ins[2] = { mi, last - 1 };
    for (int q = 0; q < 2; q++) {
        int id = ins[q], i = w;
        double e = L->elist[id];
        while (i > 0 && L->elist[L->iord[i - 1]] < e) { L->iord[i] = L->iord[i - 1]; i--; }
        L->iord[i] = id;
        w++;
    }
    *maxerr = L->iord[*nrmax - 1];
    *errmax = L->elist[*maxerr];
}

GkResult gk_qags_machine(GkSampleMachine f, void* ctx, double a, double b,
                         double abstol, double reltol, int limit) {
    const double epmach = DBL_EPSILON, uflow = DBL_MIN, oflow = DBL_MAX;
    GkResult R;
    R.value = 0.0; R.abs_err = 0.0; R.status = GK_NONNUMERIC;
    R.n_eval = 0; R.n_subdiv = 0;
    if (limit <= 0) limit = 50;
    if (abstol < 0.0) abstol = 0.0;
    if (reltol < 0.0) reltol = 0.0;
    if (abstol <= 0.0 && reltol < fmax(50.0 * epmach, 0.5e-28))
        reltol = fmax(50.0 * epmach, 0.5e-28);   /* dqagse's ier = 6 bound */

    double _Complex result, area1, area2;
    double abserr, defabs, resasc, error1, error2, defab1, defab2, rabs1, rabs2;
    if (!qk21(f, ctx, a, b, &result, &abserr, &defabs, &resasc, &R.n_eval))
        return R;
    double dres = cabs(result);
    double errbnd = fmax(abstol, reltol * dres);
    int ier = 0;
    if (abserr <= 100.0 * epmach * defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || (abserr <= errbnd && abserr != resasc) || abserr == 0.0) {
        R.value = result; R.abs_err = abserr;
        R.status = (ier == 0) ? GK_OK : GK_NOCONV;
        return R;
    }

    QagsList L;
    L.alist = malloc((size_t)limit * sizeof(double));
    L.blist = malloc((size_t)limit * sizeof(double));
    L.elist = malloc((size_t)limit * sizeof(double));
    L.rlist = malloc((size_t)limit * sizeof(double _Complex));
    L.iord  = malloc((size_t)limit * sizeof(int));
    if (!L.alist || !L.blist || !L.elist || !L.rlist || !L.iord) {
        free(L.alist); free(L.blist); free(L.elist); free(L.rlist); free(L.iord);
        R.value = result; R.abs_err = abserr; R.status = GK_NOCONV;
        return R;
    }
    L.alist[0] = a; L.blist[0] = b; L.rlist[0] = result; L.elist[0] = abserr;
    L.iord[0] = 0;

    double _Complex rlist2[QAGS_LIMEXP + 2], res3la[3] = { 0.0, 0.0, 0.0 };
    rlist2[0] = result;
    double errmax = abserr, errsum = abserr, erlarg = 0.0, ertest = 0.0;
    double correc = 0.0, small = 0.0, erlast;
    double _Complex area = result;
    int maxerr = 0, nrmax = 1, nres = 0, numrl2 = 2, ktmin = 0, ierro = 0;
    int iroff1 = 0, iroff2 = 0, iroff3 = 0;
    bool extrap = false, noext = false;
    abserr = oflow;
    int ksgn = (dres >= (1.0 - 50.0 * epmach) * defabs) ? 1 : -1;
    int last;
    bool sample_fail = false;
    bool sum_up = false;          /* dqagse's label 115: result = sum of rlist */

    for (last = 2; last <= limit; last++) {
        /* bisect the subinterval with the nrmax-th largest error estimate */
        double a1 = L.alist[maxerr], b1 = 0.5 * (L.alist[maxerr] + L.blist[maxerr]);
        double a2 = b1, b2 = L.blist[maxerr];
        erlast = errmax;
        if (!qk21(f, ctx, a1, b1, &area1, &error1, &rabs1, &defab1, &R.n_eval)
            || !qk21(f, ctx, a2, b2, &area2, &error2, &rabs2, &defab2, &R.n_eval)) {
            /* a sample of a half failed (a point the integrand cannot be
             * evaluated at): the list is intact, keep the best so far */
            sample_fail = true; last--; break;
        }
        R.n_subdiv++;
        double _Complex area12 = area1 + area2;
        double erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - L.rlist[maxerr];
        if (!(defab1 == error1 || defab2 == error2)) {
            if (cabs(L.rlist[maxerr] - area12) <= 1e-5 * cabs(area12)
                && erro12 >= 0.99 * errmax) {
                if (extrap) iroff2++; else iroff1++;
            }
            if (last > 10 && erro12 > errmax) iroff3++;
        }
        L.rlist[maxerr] = area1;
        L.rlist[last - 1] = area2;
        errbnd = fmax(abstol, reltol * cabs(area));
        /* roundoff, the limit, or a too narrow interval */
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        if (last == limit) ier = 1;
        if (fmax(fabs(a1), fabs(b2)) <= (1.0 + 100.0 * epmach) * (fabs(a2) + 1000.0 * uflow))
            ier = 4;
        if (error2 > error1) {
            L.alist[maxerr] = a2; L.alist[last - 1] = a1; L.blist[last - 1] = b1;
            L.rlist[maxerr] = area2; L.rlist[last - 1] = area1;
            L.elist[maxerr] = error2; L.elist[last - 1] = error1;
        } else {
            L.alist[last - 1] = a2; L.blist[maxerr] = b1; L.blist[last - 1] = b2;
            L.elist[maxerr] = error1; L.elist[last - 1] = error2;
        }
        qags_sort(&L, last, &nrmax, &maxerr, &errmax);
        if (errsum <= errbnd) { sum_up = true; break; }
        if (ier != 0) break;
        if (last == 2) {
            small = fabs(b - a) * 0.375;
            erlarg = errsum;
            ertest = errbnd;
            rlist2[1] = area;
            continue;
        }
        if (noext) continue;
        erlarg -= erlast;
        if (fabs(b1 - a1) > small) erlarg += erro12;
        if (!extrap) {
            /* is the interval to be bisected next the smallest one? */
            if (fabs(L.blist[maxerr] - L.alist[maxerr]) > small) continue;
            extrap = true;
            nrmax = 2;
        }
        if (ierro != 3 && erlarg > ertest) {
            /* the smallest interval has the largest error: before
             * extrapolating, bisect the large intervals first */
            bool large = false;
            int jupbnd = (last > 2 + limit / 2) ? limit + 3 - last : last;
            for (int k = nrmax; k <= jupbnd && nrmax <= last; k++) {
                maxerr = L.iord[nrmax - 1];
                errmax = L.elist[maxerr];
                if (fabs(L.blist[maxerr] - L.alist[maxerr]) > small) { large = true; break; }
                nrmax++;
            }
            if (large) continue;
            if (nrmax > last) { nrmax = last; maxerr = L.iord[nrmax - 1]; errmax = L.elist[maxerr]; }
        }
        /* perform extrapolation */
        numrl2++;
        rlist2[numrl2 - 1] = area;
        double _Complex reseps; double abseps;
        qags_qelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ktmin++;
        if (ktmin > 5 && abserr < 1e-3 * errsum) ier = 5;
        if (abseps < abserr) {
            ktmin = 0;
            abserr = abseps;
            result = reseps;
            correc = erlarg;
            ertest = fmax(abstol, reltol * cabs(reseps));
            if (abserr <= ertest) break;
        }
        /* prepare bisection of the smallest interval */
        if (numrl2 == 1) noext = true;
        if (ier == 5) break;
        maxerr = L.iord[0];
        errmax = L.elist[maxerr];
        nrmax = 1;
        extrap = false;
        small *= 0.5;
        erlarg = errsum;
    }
    if (last > limit) last = limit;

    /* final result and error estimate (dqagse's labels 100-130) */
    if (!sum_up) {
        if (abserr == oflow) {
            sum_up = true;
        } else {
            bool test_div = true;
            if (ier + ierro != 0 || sample_fail) {
                if (ierro == 3) abserr += correc;
                if (ier == 0) ier = 3;
                if (cabs(result) != 0.0 && cabs(area) != 0.0) {
                    if (abserr / cabs(result) > errsum / cabs(area)) sum_up = true;
                } else if (abserr > errsum) {
                    sum_up = true;
                } else if (cabs(area) == 0.0) {
                    test_div = false;
                }
            }
            if (!sum_up && test_div
                && !(ksgn == -1 && fmax(cabs(result), cabs(area)) <= defabs * 0.01)) {
                /* test on divergence */
                double _Complex q = result / area;
                double rat = (cimag(q) == 0.0) ? creal(q) : cabs(q);
                if (0.01 > rat || rat > 100.0 || errsum > cabs(area)) ier = 6;
            }
        }
    }
    if (sum_up) {
        /* the sum of the subinterval results, compensated */
        double sr = 0.0, si = 0.0, cr = 0.0, ci = 0.0;
        for (int k = 0; k < last; k++) {
            double xr = creal(L.rlist[k]), xi = cimag(L.rlist[k]);
            double tr = sr + xr, ti = si + xi;
            cr += (fabs(sr) >= fabs(xr)) ? (sr - tr) + xr : (xr - tr) + sr;
            ci += (fabs(si) >= fabs(xi)) ? (si - ti) + xi : (xi - ti) + si;
            sr = tr; si = ti;
        }
        result = (sr + cr) + (si + ci) * I;
        abserr = errsum;
        if (sample_fail && ier == 0) ier = 3;
    }
    if (ier > 2) ier--;
    R.value = result;
    R.abs_err = abserr;
    R.status = (ier == 0) ? GK_OK : GK_NOCONV;
    free(L.alist); free(L.blist); free(L.elist); free(L.rlist); free(L.iord);
    return R;
}
