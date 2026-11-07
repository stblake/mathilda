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
