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

#include <float.h>
#include <math.h>
#include <stdlib.h>

/* ---- QK15: 7-point Gauss embedded in a 15-point Kronrod rule ---- */
/* Abscissae of the 15-point Kronrod rule (nonnegative half; xgk[7] = centre). */
static const double XGK[8] = {
    0.991455371120813, 0.949107912342759, 0.864864423359769,
    0.741531185599394, 0.586087235467691, 0.405845151377397,
    0.207784955007898, 0.000000000000000
};
/* Weights of the 15-point Kronrod rule. */
static const double WGK[8] = {
    0.022935322010529, 0.063092092629979, 0.104790010322250,
    0.140653259715525, 0.169004726639267, 0.190350578064785,
    0.204432940075298, 0.209482141084728
};
/* Weights of the 7-point Gauss rule (wg[3] is the centre weight). */
static const double WG[4] = {
    0.129484966168870, 0.279705391489277, 0.381830050505119,
    0.417959183673469
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

GkResult gk_integrate_machine(GkSampleMachine f, void* ctx,
                              double a, double b,
                              double abstol, double reltol,
                              int max_subdiv, long max_eval,
                              bool extrapolate) {
    GkResult R;
    R.value = 0.0; R.abs_err = 0.0; R.status = GK_NONNUMERIC;
    R.n_eval = 0; R.n_subdiv = 0;

    if (max_subdiv < 0) max_subdiv = 0;
    if (max_eval <= 0) max_eval = (long)(max_subdiv + 1) * 15L + 30L;
    if (abstol < 0.0) abstol = 0.0;
    if (reltol < 0.0) reltol = 0.0;

    Panel p0;
    p0.a = a; p0.b = b;
    if (!qk15(f, ctx, a, b, &p0.result, &p0.err, NULL, NULL, &R.n_eval))
        return R;   /* non-numeric on the whole interval */

    Heap h = { NULL, 0, 0 };
    if (!heap_push(&h, p0)) { R.status = GK_NONNUMERIC; return R; }

    double _Complex total = p0.result;
    double total_err = p0.err;

    /* Extrapolation: feed the running total after each bisection to Wynn. */
    double _Complex* seq = NULL;
    int seq_n = 0, seq_cap = 0;

    R.status = (total_err <= fmax(abstol, reltol * cabs(total))) ? GK_OK : GK_NOCONV;
    while (R.status != GK_OK && R.n_subdiv < max_subdiv && R.n_eval + 30 <= max_eval) {
        double tol = fmax(abstol, reltol * cabs(total));
        if (total_err <= tol) { R.status = GK_OK; break; }

        Panel worst = heap_pop(&h);
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
