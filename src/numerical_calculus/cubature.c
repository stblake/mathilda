/*
 * cubature.c — adaptive Genz-Malik multidimensional cubature (see cubature.h).
 *
 * The degree-7 rule (with embedded degree-5 error estimate) for an axis-aligned
 * box follows Genz & Malik, "An adaptive algorithm for numerical integration
 * over an N-dimensional rectangular region" (J. Comput. Appl. Math. 6, 1980),
 * in the normalisation used by Steven G. Johnson's public-domain `cubature`
 * library.  Per region of half-width h about centre c, the rule evaluates
 * 2d² + 2d + 1 + 2^d points: the centre, two axial shells at λ2 and λ4, the
 * λ4 coordinate pairs, and the 2^d λ5 corners.  The fourth difference along
 * each axis picks the subdivision direction.
 */

#include "cubature.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- Genz-Malik degree-7 / degree-5 rule on one region ---- */

/* Generator distances (squared values are 9/70, 9/10, 9/19). */
#define GM_LAMBDA2 0.3585685828003180919906451539079374954541
#define GM_LAMBDA4 0.9486832980505137995996680633298155601160
#define GM_LAMBDA5 0.6882472016116852977216287342936235251269

/* Evaluate the rule on the box of centre c[], half-width h[].  Writes the
 * degree-7 estimate to *val, |deg7 - deg5| to *err, and the axis of largest
 * fourth difference to *split.  `p` is caller-supplied scratch of length d.
 * Returns false on a non-numeric sample. */
static bool gm_rule(McSampleMachine f, void* ctx,
                    const double* c, const double* h, size_t d, double* p,
                    double _Complex* val, double* err, int* split) {
    const double lambda2 = GM_LAMBDA2, lambda4 = GM_LAMBDA4, lambda5 = GM_LAMBDA5;
    const double weight2 = 980.0 / 6561.0;
    const double weight4 = 200.0 / 19683.0;
    const double weightE2 = 245.0 / 486.0;
    const double weightE4 = 25.0 / 729.0;
    const double ratio = (lambda2 * lambda2) / (lambda4 * lambda4);
    const double dd = (double)d;
    const double weight1 = (12824.0 - 9120.0 * dd + 400.0 * dd * dd) / 19683.0;
    const double weight3 = (1820.0 - 400.0 * dd) / 19683.0;
    const double weight5 = 6859.0 / 19683.0 / (double)(1UL << d);
    const double weightE1 = (729.0 - 950.0 * dd + 50.0 * dd * dd) / 729.0;
    const double weightE3 = (265.0 - 100.0 * dd) / 1458.0;

    double _Complex f0;
    memcpy(p, c, d * sizeof(double));
    if (!f(ctx, p, d, &f0)) return false;

    double _Complex sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;
    double maxdiff = -1.0;
    *split = 0;

    /* Axial shells at λ2 and λ4 (also yields the fourth difference per axis). */
    for (size_t i = 0; i < d; i++) {
        double _Complex fa, fb, fc, fd;
        memcpy(p, c, d * sizeof(double));
        p[i] = c[i] + lambda2 * h[i]; if (!f(ctx, p, d, &fa)) return false;
        p[i] = c[i] - lambda2 * h[i]; if (!f(ctx, p, d, &fb)) return false;
        p[i] = c[i] + lambda4 * h[i]; if (!f(ctx, p, d, &fc)) return false;
        p[i] = c[i] - lambda4 * h[i]; if (!f(ctx, p, d, &fd)) return false;
        sum2 += fa + fb;
        sum3 += fc + fd;
        double diff = cabs((fa + fb - 2.0 * f0) - ratio * (fc + fd - 2.0 * f0));
        if (diff > maxdiff) { maxdiff = diff; *split = (int)i; }
    }

    /* λ4 coordinate pairs (four sign combinations each). */
    for (size_t i = 0; i < d; i++) {
        for (size_t j = i + 1; j < d; j++) {
            double _Complex f1, f2, f3, f4;
            memcpy(p, c, d * sizeof(double));
            p[i] = c[i] + lambda4 * h[i]; p[j] = c[j] + lambda4 * h[j];
            if (!f(ctx, p, d, &f1)) return false;
            p[j] = c[j] - lambda4 * h[j]; if (!f(ctx, p, d, &f2)) return false;
            p[i] = c[i] - lambda4 * h[i]; if (!f(ctx, p, d, &f3)) return false;
            p[j] = c[j] + lambda4 * h[j]; if (!f(ctx, p, d, &f4)) return false;
            sum4 += f1 + f2 + f3 + f4;
        }
    }

    /* 2^d λ5 corners. */
    unsigned long ncorner = 1UL << d;
    for (unsigned long m = 0; m < ncorner; m++) {
        for (size_t i = 0; i < d; i++) {
            double s = ((m >> i) & 1UL) ? -1.0 : 1.0;
            p[i] = c[i] + s * lambda5 * h[i];
        }
        double _Complex fc;
        if (!f(ctx, p, d, &fc)) return false;
        sum5 += fc;
    }

    double vol = 1.0;
    for (size_t i = 0; i < d; i++) vol *= 2.0 * h[i];

    double _Complex res7 = vol * (weight1 * f0 + weight2 * sum2 + weight3 * sum3
                                  + weight4 * sum4 + weight5 * sum5);
    double _Complex res5 = vol * (weightE1 * f0 + weightE2 * sum2
                                  + weightE3 * sum3 + weightE4 * sum4);
    *val = res7;
    *err = cabs(res7 - res5);
    return true;
}

/* Per-region evaluation count for budget bookkeeping. */
static long gm_npoints(size_t d) {
    return 1L + 4L * (long)d + 2L * (long)d * ((long)d - 1L) + (long)(1UL << d);
}

/* ---- region max-heap (keyed by error) ---- */

typedef struct {
    double*       c;     /* centre,     length d (owned) */
    double*       h;     /* half-width, length d (owned) */
    double _Complex val;
    double        err;
    int           split; /* axis to bisect next */
} Region;

typedef struct { Region* v; int n, cap; } Heap;

static bool heap_push(Heap* hp, Region r) {
    if (hp->n == hp->cap) {
        int nc = hp->cap ? hp->cap * 2 : 64;
        Region* nv = realloc(hp->v, (size_t)nc * sizeof(Region));
        if (!nv) return false;
        hp->v = nv; hp->cap = nc;
    }
    int i = hp->n++;
    hp->v[i] = r;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (hp->v[p].err >= hp->v[i].err) break;
        Region t = hp->v[p]; hp->v[p] = hp->v[i]; hp->v[i] = t; i = p;
    }
    return true;
}

static Region heap_pop(Heap* hp) {
    Region top = hp->v[0];
    hp->n--;
    hp->v[0] = hp->v[hp->n];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, m = i;
        if (l < hp->n && hp->v[l].err > hp->v[m].err) m = l;
        if (r < hp->n && hp->v[r].err > hp->v[m].err) m = r;
        if (m == i) break;
        Region t = hp->v[m]; hp->v[m] = hp->v[i]; hp->v[i] = t; i = m;
    }
    return top;
}

static void heap_free(Heap* hp) {
    for (int i = 0; i < hp->n; i++) { free(hp->v[i].c); free(hp->v[i].h); }
    free(hp->v);
}

/* ---- adaptive driver ---- */

CubStatus cub_integrate_machine(McSampleMachine f, void* ctx,
                                const double* a, const double* b, size_t d,
                                double abstol, double reltol, long max_eval,
                                double _Complex* result, double* abserr) {
    if (d < 2 || d > 15) return CUB_NONNUMERIC;
    if (max_eval <= 0) max_eval = 1000000;
    if (abstol < 0.0) abstol = 0.0;
    if (reltol < 0.0) reltol = 0.0;

    double* p = malloc(d * sizeof(double));   /* shared scratch eval point */
    if (!p) return CUB_NONNUMERIC;
    long per = gm_npoints(d);

    Region r0;
    r0.c = malloc(d * sizeof(double));
    r0.h = malloc(d * sizeof(double));
    if (!r0.c || !r0.h) { free(p); free(r0.c); free(r0.h); return CUB_NONNUMERIC; }
    for (size_t i = 0; i < d; i++) {
        r0.c[i] = 0.5 * (a[i] + b[i]);
        r0.h[i] = 0.5 * (b[i] - a[i]);
    }
    if (!gm_rule(f, ctx, r0.c, r0.h, d, p, &r0.val, &r0.err, &r0.split)) {
        free(p); free(r0.c); free(r0.h);
        return CUB_NONNUMERIC;
    }

    double _Complex total = r0.val;
    double total_err = r0.err;
    /* Σ|region value| dominates |total| and stays meaningful when the signed
     * integral cancels to ~0, so the relative test does not chase an
     * unreachable tolerance (e.g. ∫∫ sin(x+y) over [0,π]² = 0). */
    double total_absval = cabs(r0.val);
    long n_eval = per;

    Heap hp = { NULL, 0, 0 };
    if (!heap_push(&hp, r0)) { free(p); free(r0.c); free(r0.h); return CUB_NONNUMERIC; }

    CubStatus status = CUB_NOCONV;
    if (total_err <= fmax(abstol, reltol * total_absval)) status = CUB_OK;

    while (status != CUB_OK && n_eval + 2 * per <= max_eval) {
        if (total_err <= fmax(abstol, reltol * total_absval)) { status = CUB_OK; break; }

        Region worst = heap_pop(&hp);
        int ax = worst.split;

        Region left, right;
        left.c = malloc(d * sizeof(double));  left.h = malloc(d * sizeof(double));
        right.c = malloc(d * sizeof(double)); right.h = malloc(d * sizeof(double));
        if (!left.c || !left.h || !right.c || !right.h) {
            free(left.c); free(left.h); free(right.c); free(right.h);
            free(worst.c); free(worst.h);
            break;   /* out of memory: report best so far */
        }
        memcpy(left.c, worst.c, d * sizeof(double));
        memcpy(left.h, worst.h, d * sizeof(double));
        memcpy(right.c, worst.c, d * sizeof(double));
        memcpy(right.h, worst.h, d * sizeof(double));
        double hh = 0.5 * worst.h[ax];
        left.h[ax] = hh;  left.c[ax]  = worst.c[ax] - hh;
        right.h[ax] = hh; right.c[ax] = worst.c[ax] + hh;

        bool okl = gm_rule(f, ctx, left.c, left.h, d, p, &left.val, &left.err, &left.split);
        bool okr = okl && gm_rule(f, ctx, right.c, right.h, d, p, &right.val, &right.err, &right.split);
        if (!okl || !okr) {
            free(left.c); free(left.h); free(right.c); free(right.h);
            /* re-insert the parent so its contribution is not lost */
            heap_push(&hp, worst);
            status = CUB_NONNUMERIC;
            break;
        }
        n_eval += 2 * per;

        total += (left.val + right.val) - worst.val;
        total_err += (left.err + right.err) - worst.err;
        if (total_err < 0.0) total_err = 0.0;
        total_absval += (cabs(left.val) + cabs(right.val)) - cabs(worst.val);
        if (total_absval < 0.0) total_absval = 0.0;
        free(worst.c); free(worst.h);

        if (!heap_push(&hp, left) || !heap_push(&hp, right)) break;
    }

    if (status != CUB_NONNUMERIC && total_err <= fmax(abstol, reltol * total_absval))
        status = CUB_OK;

    *result = total;
    *abserr = total_err;
    heap_free(&hp);
    free(p);
    return status;
}
