/*
 * denint.c — tanh-sinh (finite) and sinh-sinh (whole-line) double-exponential
 * quadrature.  See denint.h.  Structure mirrors dequad.c: level-by-level step
 * halving, symmetric tail accumulation with a relative-floor cut and a
 * roundoff-noise blow-up guard, machine and MPFR variants.
 */

#include "denint.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEN_TS_TMAX   4.5      /* tanh-sinh |t| cap (weight underflows beyond)  */
#define DEN_SS_TMAX   5.0      /* sinh-sinh |t| cap (x overflows beyond)        */
#define DEN_REL_FLOOR 1e-17    /* tail cut: term < floor * running max          */

/* ================================================================== *
 *  Machine precision                                                  *
 * ================================================================== */

/* tanh-sinh transformed sample on [a,b] at parameter t.  The abscissa is built
 * from the accurate endpoint offset dist = (b-a)/2 · (1 - tanh|u|) so it never
 * collapses onto a singular endpoint while the weight is still significant. */
static bool den_ts_term_machine(DeQuadSampleMachine f, void* ctx,
                                double a, double b, double t,
                                double _Complex* term) {
    double half = 0.5 * (b - a);
    double u = (M_PI / 2.0) * sinh(t);
    double au = fabs(u);
    if (au > 350.0) return false;
    double eu = exp(-2.0 * au);
    double onep = 1.0 + eu;
    double sech2 = 4.0 * eu / (onep * onep);            /* sech^2(u)             */
    double weight = (M_PI / 2.0) * cosh(t) * sech2;     /* x'(t) / half          */
    if (!isfinite(weight) || weight == 0.0) return false;
    double dist = half * (2.0 * eu / onep);             /* half·(1 - tanh|u|)    */
    double x = (t >= 0.0) ? (b - dist) : (a + dist);
    if (x <= a || x >= b) return false;                 /* rounded onto endpoint */
    double _Complex fv;
    if (!f(ctx, x, &fv)) return false;
    double _Complex tm = fv * (weight * half);
    if (!isfinite(creal(tm)) || !isfinite(cimag(tm))) return false;
    *term = tm;
    return true;
}

/* sinh-sinh transformed sample on (-∞,∞):  x = sinh((π/2) sinh t). */
static bool den_ss_term_machine(DeQuadSampleMachine f, void* ctx, double t,
                                double _Complex* term) {
    double u = (M_PI / 2.0) * sinh(t);
    if (fabs(u) > 350.0) return false;
    double x = sinh(u);
    if (!isfinite(x)) return false;
    double weight = (M_PI / 2.0) * cosh(t) * cosh(u);   /* x'(t)                 */
    if (!isfinite(weight)) return false;
    double _Complex fv;
    if (!f(ctx, x, &fv)) return false;
    double _Complex tm = fv * weight;
    if (!isfinite(creal(tm)) || !isfinite(cimag(tm))) return false;
    *term = tm;
    return true;
}

/* Shared machine DE driver, parameterised by the transform's term function. */
typedef bool (*DenTermMachine)(DeQuadSampleMachine, void*, double t, double _Complex*);

/* Wrappers giving the two transforms a uniform (t -> term) signature. */
typedef struct { DeQuadSampleMachine f; void* ctx; double a, b; } DenTsCtx;
static bool den_ts_thunk(DeQuadSampleMachine f, void* ctx, double t, double _Complex* o) {
    DenTsCtx* c = (DenTsCtx*)ctx; (void)f;
    return den_ts_term_machine(c->f, c->ctx, c->a, c->b, t, o);
}
static bool den_ss_thunk(DeQuadSampleMachine f, void* ctx, double t, double _Complex* o) {
    DenTsCtx* c = (DenTsCtx*)ctx; (void)f;
    return den_ss_term_machine(c->f, c->ctx, t, o);
}

static bool den_run_machine(DenTermMachine term, void* tctx, double tmax,
                            double reltol, int max_levels,
                            double _Complex* result, double* abserr) {
    if (reltol <= 0.0) reltol = 1e-13;
    if (max_levels < 1) max_levels = 1;
    double h = 1.0;
    double _Complex prev = 0.0, best = 0.0;
    bool have_prev = false, any = false;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++, h *= 0.5) {
        double _Complex sum = 0.0, c0;
        if (term(NULL, tctx, 0.0, &c0)) { sum = c0; any = true; }
        double maxterm = cabs(sum);
        const double lcenter = de_log2_mag(maxterm);   /* the t = 0 term */
        for (int side = -1; side <= 1; side += 2) {
            /* Floor cut over a t-window, noise cut confirmed before it is
             * taken: an interior zero is not the tail (DeTailGuard, dequad.h). */
            DeTailGuard g; de_tail_guard_init(&g, lcenter);
            double _Complex ckpt = 0.0;
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > tmax) break;
                double _Complex tm;
                if (!term(NULL, tctx, t, &tm)) break;
                any = true;
                int act = de_tail_step(&g, fabs(t), cabs(tm), &maxterm,
                                       DEN_REL_FLOOR);
                if (act == DE_TAIL_ROLLBACK_STOP) { sum = ckpt; break; }
                if (act == DE_TAIL_CHECKPOINT) ckpt = sum;
                sum += tm;
                if (act == DE_TAIL_ADD_STOP) break;
            }
        }
        double _Complex Ival = h * sum;
        best = Ival;
        if (have_prev) {
            double diff = cabs(Ival - prev);
            *abserr = diff;
            if (diff <= reltol * cabs(Ival) + 1e-300) {
                if (!any) break;               /* no numeric sample at all */
                *result = Ival; return true;
            }
        }
        prev = Ival;
        have_prev = true;
    }
    *result = any ? best : NAN + NAN * I;
    return false;
}

bool denint_tanhsinh_machine(DeQuadSampleMachine f, void* ctx, double a, double b,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr) {
    DenTsCtx c = { f, ctx, a, b };
    return den_run_machine(den_ts_thunk, &c, DEN_TS_TMAX, reltol, max_levels, result, abserr);
}

bool denint_sinhsinh_machine(DeQuadSampleMachine f, void* ctx,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr) {
    DenTsCtx c = { f, ctx, 0.0, 0.0 };
    return den_run_machine(den_ss_thunk, &c, DEN_SS_TMAX, reltol, max_levels, result, abserr);
}

/* ================================================================== *
 *  MPFR precision                                                     *
 * ================================================================== */

#ifdef USE_MPFR

typedef struct {
    mpfr_t halfpi, t, s, u, au, eu, onep, tmp, weight, x, fr, fi, half, a, b;
    bool finite;   /* false => tanh-sinh (uses a,b,half); n/a for sinh-sinh */
} DenMpfrScratch;

static void den_scratch_init(DenMpfrScratch* sc, mpfr_prec_t p) {
    mpfr_inits2(p, sc->halfpi, sc->t, sc->s, sc->u, sc->au, sc->eu, sc->onep,
                sc->tmp, sc->weight, sc->x, sc->fr, sc->fi, sc->half,
                sc->a, sc->b, (mpfr_ptr)0);
    mpfr_const_pi(sc->halfpi, MPFR_RNDN);
    mpfr_div_2ui(sc->halfpi, sc->halfpi, 1, MPFR_RNDN);   /* π/2 */
}
static void den_scratch_clear(DenMpfrScratch* sc) {
    mpfr_clears(sc->halfpi, sc->t, sc->s, sc->u, sc->au, sc->eu, sc->onep,
                sc->tmp, sc->weight, sc->x, sc->fr, sc->fi, sc->half,
                sc->a, sc->b, (mpfr_ptr)0);
}

/* tanh-sinh MPFR term at parameter tv. */
static bool den_ts_term_mpfr(DeQuadSampleMPFR f, void* ctx, DenMpfrScratch* sc,
                             double tv, mpfr_t term_re, mpfr_t term_im) {
    mpfr_set_d(sc->t, tv, MPFR_RNDN);
    /* sinh t and cosh t from one exponential (the transform otherwise costs
     * three exps a node, as much as a cheap integrand) */
    mpfr_sinh_cosh(sc->s, sc->x, sc->t, MPFR_RNDN);       /* x holds cosh t */
    mpfr_mul(sc->u, sc->s, sc->halfpi, MPFR_RNDN);        /* u = (π/2) sinh t */
    mpfr_abs(sc->au, sc->u, MPFR_RNDN);
    if (mpfr_get_d(sc->au, MPFR_RNDN) > 1e9) return false;
    mpfr_mul_2ui(sc->tmp, sc->au, 1, MPFR_RNDN);          /* 2|u| */
    mpfr_neg(sc->tmp, sc->tmp, MPFR_RNDN);
    mpfr_exp(sc->eu, sc->tmp, MPFR_RNDN);                 /* eu = exp(-2|u|) */
    mpfr_add_ui(sc->onep, sc->eu, 1, MPFR_RNDN);          /* 1 + eu */
    /* weight = (π/2) cosh t · sech^2(u) = (π/2) cosh t · 4 eu/(1+eu)^2 */
    mpfr_mul(sc->weight, sc->halfpi, sc->x, MPFR_RNDN);
    mpfr_mul_2ui(sc->tmp, sc->eu, 2, MPFR_RNDN);          /* 4 eu */
    mpfr_mul(sc->weight, sc->weight, sc->tmp, MPFR_RNDN);
    mpfr_mul(sc->tmp, sc->onep, sc->onep, MPFR_RNDN);     /* (1+eu)^2 */
    mpfr_div(sc->weight, sc->weight, sc->tmp, MPFR_RNDN);
    if (!mpfr_number_p(sc->weight) || mpfr_zero_p(sc->weight)) return false;
    /* dist = half · 2 eu/(1+eu);  x = (t>=0)? b-dist : a+dist */
    mpfr_mul_2ui(sc->tmp, sc->eu, 1, MPFR_RNDN);          /* 2 eu */
    mpfr_div(sc->tmp, sc->tmp, sc->onep, MPFR_RNDN);
    mpfr_mul(sc->tmp, sc->tmp, sc->half, MPFR_RNDN);      /* dist */
    if (tv >= 0.0) mpfr_sub(sc->x, sc->b, sc->tmp, MPFR_RNDN);
    else           mpfr_add(sc->x, sc->a, sc->tmp, MPFR_RNDN);
    if (mpfr_lessequal_p(sc->x, sc->a) || mpfr_greaterequal_p(sc->x, sc->b)) return false;
    if (!f(ctx, sc->x, sc->fr, sc->fi)) return false;
    if (!mpfr_number_p(sc->fr) || !mpfr_number_p(sc->fi)) return false;
    mpfr_mul(sc->weight, sc->weight, sc->half, MPFR_RNDN);  /* x'(t) = weight·half */
    mpfr_mul(term_re, sc->fr, sc->weight, MPFR_RNDN);
    mpfr_mul(term_im, sc->fi, sc->weight, MPFR_RNDN);
    return true;
}

/* sinh-sinh MPFR term at parameter tv. */
static bool den_ss_term_mpfr(DeQuadSampleMPFR f, void* ctx, DenMpfrScratch* sc,
                             double tv, mpfr_t term_re, mpfr_t term_im) {
    mpfr_set_d(sc->t, tv, MPFR_RNDN);
    mpfr_sinh(sc->s, sc->t, MPFR_RNDN);
    mpfr_mul(sc->u, sc->s, sc->halfpi, MPFR_RNDN);        /* u = (π/2) sinh t */
    mpfr_sinh(sc->x, sc->u, MPFR_RNDN);                   /* x = sinh u */
    if (!mpfr_number_p(sc->x)) return false;
    mpfr_cosh(sc->tmp, sc->t, MPFR_RNDN);
    mpfr_mul(sc->weight, sc->halfpi, sc->tmp, MPFR_RNDN);
    mpfr_cosh(sc->tmp, sc->u, MPFR_RNDN);
    mpfr_mul(sc->weight, sc->weight, sc->tmp, MPFR_RNDN); /* (π/2) cosh t cosh u */
    if (!mpfr_number_p(sc->weight)) return false;
    if (!f(ctx, sc->x, sc->fr, sc->fi)) return false;
    if (!mpfr_number_p(sc->fr) || !mpfr_number_p(sc->fi)) return false;
    mpfr_mul(term_re, sc->fr, sc->weight, MPFR_RNDN);
    mpfr_mul(term_im, sc->fi, sc->weight, MPFR_RNDN);
    return true;
}

/* The MPFR driver nests its levels: level L halves the step and evaluates only
 * the odd multiples of the new step, so the samples of every coarser level are
 * reused (the running sum is the sum over all nodes so far; I = h * sum).  A
 * side's tail is cut on level 0's grid, through the shared guard; a finer level
 * only fills in nodes inside the span the side has already settled on (its new
 * odd nodes are still run through a fresh guard, so a noise rise between them
 * cuts too), which keeps every level a trapezoidal sum over one common range.
 *
 * Convergence: the double-exponential error falls roughly quadratically per
 * halving, so once the level differences are shrinking (d_L < d_{L-1}) the
 * error of level L is estimated as d_L^2 / d_{L-1} rather than d_L, which
 * saves the final confirming level: half the samples of the run. */
static bool den_run_mpfr(DeQuadSampleMPFR f, void* ctx, DenMpfrScratch* sc,
                         bool tanhsinh, long bits, double tmax,
                         double reltol, mpfr_srcptr abstol, int max_levels,
                         mpfr_t out_re, mpfr_t out_im, double* abserr) {
    if (reltol <= 0.0) reltol = 1e-13;
    if (max_levels < 1) max_levels = 1;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    double lfloor = -(double)bits - 4.0;   /* log2 of the relative tail floor */
    double lrel = log2(reltol);
    double labs = abstol ? de_mpfr_log2_abs(abstol) : -INFINITY;

    mpfr_t sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di;
    mpfr_t ck_re, ck_im;   /* noise checkpoint (DeTailGuard) */
    mpfr_inits2(p, sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di,
                ck_re, ck_im, (mpfr_ptr)0);
    mpfr_set_ui(sum_re, 0, MPFR_RNDN);
    mpfr_set_ui(sum_im, 0, MPFR_RNDN);

    double h = 1.0;
    bool have_prev = false, converged = false;
    double lprevd = INFINITY;          /* log2 of the previous level difference */
    int stalled = 0;                   /* consecutive levels without progress  */
    bool any = false;                  /* some sample was numeric              */
    double tcut[2] = { tmax, tmax };   /* settled |t| span per side            */
    double lmax = -INFINITY;           /* log2 of the largest term so far      */
    double lcenter = -INFINITY;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++) {
        if (level == 0) {
            bool ok0 = tanhsinh ? den_ts_term_mpfr(f, ctx, sc, 0.0, tr, ti)
                                : den_ss_term_mpfr(f, ctx, sc, 0.0, tr, ti);
            if (ok0) {
                any = true;
                mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
                mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
                lmax = de_mpfr_log2_l1(tr, ti);
            }
            lcenter = lmax;                            /* the t = 0 term */
        } else {
            h *= 0.5;
        }
        int kstep = (level == 0) ? 1 : 2;   /* finer levels: odd nodes only */
        /* The guard reads this level's terms against this level's maximum, as
         * the unnested driver did: a maximum carried over from a coarser level
         * would make a side still climbing toward a peak the coarse grid already
         * sampled look decayed, then rising -- a confirmed noise cut. */
        double lmax_lvl = lcenter;
        for (int side = -1; side <= 1; side += 2) {
            int si = (side > 0);
            double span = tcut[si];
            double reached = 0.0;
            DeTailGuard g; de_tail_guard_init(&g, lcenter);   /* dequad.h */
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > tmax) break;
                /* an even node inside the settled span is already summed */
                if (kstep == 2 && (k & 1) == 0 && fabs(t) <= span) continue;
                bool ok = tanhsinh ? den_ts_term_mpfr(f, ctx, sc, t, tr, ti)
                                   : den_ss_term_mpfr(f, ctx, sc, t, tr, ti);
                if (!ok) break;
                any = true;
                double lt = de_mpfr_log2_l1(tr, ti);
                int act = de_tail_step_log2(&g, fabs(t), lt, &lmax_lvl, lfloor);
                if (lt > lmax) lmax = lt;
                if (act == DE_TAIL_ROLLBACK_STOP) {
                    mpfr_set(sum_re, ck_re, MPFR_RNDN);
                    mpfr_set(sum_im, ck_im, MPFR_RNDN);
                    break;
                }
                if (act == DE_TAIL_CHECKPOINT) {
                    mpfr_set(ck_re, sum_re, MPFR_RNDN);
                    mpfr_set(ck_im, sum_im, MPFR_RNDN);
                }
                mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
                mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
                reached = fabs(t);
                if (act == DE_TAIL_ADD_STOP) break;
            }
            /* the span the side's sum now covers on this level's grid */
            if (reached > span || level == 0 || reached + h < span) tcut[si] = reached;
        }
        mpfr_mul_d(I_re, sum_re, h, MPFR_RNDN);
        mpfr_mul_d(I_im, sum_im, h, MPFR_RNDN);
        mpfr_set(out_re, I_re, MPFR_RNDN);
        mpfr_set(out_im, I_im, MPFR_RNDN);
        if (have_prev) {
            mpfr_sub(dr, I_re, prev_re, MPFR_RNDN);
            mpfr_sub(di, I_im, prev_im, MPFR_RNDN);
            *abserr = de_mpfr_l1_d(dr, di);
            double ld = de_mpfr_log2_l1(dr, di);
            /* quadratic-convergence error estimate once the differences shrink */
            double lest = (level >= 3 && ld < lprevd) ? 2.0 * ld - lprevd : ld;
            double lI = de_mpfr_log2_l1(I_re, I_im);
            if (ld == -INFINITY || lest <= lrel + lI || lest <= labs
                || de_mpfr_level_converged(dr, di, I_re, I_im, reltol,
                                           lmax + log2(h), bits)) {
                converged = true; break;
            }
            /* Noise floor: an integrand whose values carry fewer digits than
             * the goal (machine numbers inside a WorkingPrecision -> 25 run)
             * stops improving; each further level doubles the cost for
             * nothing, so after two levels without a bit of progress the best
             * estimate is returned unconverged (the caller warns). */
            if (level >= 4 && !(ld < lprevd - 1.0)) {
                if (++stalled >= 2) break;
            } else {
                stalled = 0;
            }
            lprevd = ld;
        }
        mpfr_set(prev_re, I_re, MPFR_RNDN);
        mpfr_set(prev_im, I_im, MPFR_RNDN);
        have_prev = true;
    }

    /* No sample was a number: there is no estimate (an all-failing sum of
     * zero is not the integral 0). */
    if (!any) {
        mpfr_set_nan(out_re); mpfr_set_nan(out_im);
        converged = false;
    }
    mpfr_clears(sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di,
                ck_re, ck_im, (mpfr_ptr)0);
    return converged;
}

bool denint_tanhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          const mpfr_t a, const mpfr_t b,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr) {
    return denint_tanhsinh_mpfr_abs(f, ctx, a, b, bits, reltol, NULL, max_levels,
                                    out_re, out_im, abserr);
}

bool denint_tanhsinh_mpfr_abs(DeQuadSampleMPFR f, void* ctx,
                              const mpfr_t a, const mpfr_t b,
                              long bits, double reltol, mpfr_srcptr abstol,
                              int max_levels,
                              mpfr_t out_re, mpfr_t out_im, double* abserr) {
    DenMpfrScratch sc; den_scratch_init(&sc, (mpfr_prec_t)bits);
    mpfr_set(sc.a, a, MPFR_RNDN);
    mpfr_set(sc.b, b, MPFR_RNDN);
    mpfr_sub(sc.half, b, a, MPFR_RNDN);
    mpfr_div_2ui(sc.half, sc.half, 1, MPFR_RNDN);   /* (b-a)/2 */
    double tmax = 4.0 + 0.25 * log((double)bits);
    bool ok = den_run_mpfr(f, ctx, &sc, true, bits, tmax, reltol, abstol, max_levels,
                           out_re, out_im, abserr);
    den_scratch_clear(&sc);
    return ok;
}

bool denint_sinhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr) {
    DenMpfrScratch sc; den_scratch_init(&sc, (mpfr_prec_t)bits);
    double tmax = 4.5 + 0.25 * log((double)bits);
    bool ok = den_run_mpfr(f, ctx, &sc, false, bits, tmax, reltol, NULL, max_levels,
                           out_re, out_im, abserr);
    den_scratch_clear(&sc);
    return ok;
}

/* ================================================================== *
 *  Clenshaw-Curtis at MPFR precision (smooth pieces of a split range)  *
 * ================================================================== *
 *
 * On an interval where the integrand is analytic and varies on the
 * interval's own scale, a polynomial rule converges geometrically in the node
 * count, where tanh-sinh, built for endpoint singularities, pays ~2 tmax/h
 * nodes regardless: 344 samples for Exp[-u] on [0, 1] at 120 digits, against
 * 33 Clenshaw-Curtis nodes.  The Chebyshev points cos(j pi/N) nest under
 * doubling N, so each refinement evaluates only the N new odd-index nodes;
 * the node and weight tables are built once per (N, precision) in a cache the
 * caller owns, shared by every piece of one split integral. */

struct CcCache {
    long bits;
    int  have;                          /* levels built: N = CC_N0 << l     */
    mpfr_t* x[CC_MAX_LEVEL + 1];        /* cos(j pi / N), j = 0..N          */
    mpfr_t* w[CC_MAX_LEVEL + 1];        /* weights on [-1, 1], j = 0..N     */
};

CcCache* ccquad_cache_new(long bits) {
    CcCache* c = (CcCache*)calloc(1, sizeof *c);
    if (c) c->bits = bits;
    return c;
}

void ccquad_cache_free(CcCache* c) {
    if (!c) return;
    for (int l = 0; l < c->have; l++) {
        long n = (long)CC_N0 << l;
        for (long j = 0; j <= n; j++) { mpfr_clear(c->x[l][j]); mpfr_clear(c->w[l][j]); }
        free(c->x[l]); free(c->w[l]);
    }
    free(c);
}

/* Build level l: nodes cos(j pi/N) and the Clenshaw-Curtis weights
 *   w_j = (c_j/N) (1 - sum_{k=1}^{N/2} b_k/(4k^2 - 1) cos(2 k j pi/N)),
 * c_0 = c_N = 1, c_j = 2 otherwise; b_{N/2} = 1, b_k = 2 otherwise. */
static void cc_build(CcCache* c, int l) {
    long n = (long)CC_N0 << l;
    mpfr_prec_t p = (mpfr_prec_t)c->bits;
    c->x[l] = (mpfr_t*)malloc((size_t)(n + 1) * sizeof(mpfr_t));
    c->w[l] = (mpfr_t*)malloc((size_t)(n + 1) * sizeof(mpfr_t));
    mpfr_t pi, t, acc, coef;
    mpfr_inits2(p, pi, t, acc, coef, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    for (long j = 0; j <= n; j++) {
        mpfr_init2(c->x[l][j], p);
        mpfr_init2(c->w[l][j], p);
        if (l > 0 && (j & 1) == 0) {           /* shared with level l-1 */
            mpfr_set(c->x[l][j], c->x[l - 1][j / 2], MPFR_RNDN);
        } else {
            mpfr_mul_si(t, pi, j, MPFR_RNDN);
            mpfr_div_si(t, t, n, MPFR_RNDN);
            mpfr_cos(c->x[l][j], t, MPFR_RNDN);
        }
    }
    /* cos(m pi/N) for m = 0..2N-1, read off the nodes by symmetry */
    for (long j = 0; j <= n; j++) {
        mpfr_set_ui(acc, 1, MPFR_RNDN);
        for (long k = 1; k <= n / 2; k++) {
            long m = (2 * k * j) % (2 * n);
            mpfr_srcptr cm = (m <= n) ? c->x[l][m] : c->x[l][2 * n - m];
            mpfr_set_si(coef, (k == n / 2) ? 1 : 2, MPFR_RNDN);
            mpfr_div_si(coef, coef, 4 * k * k - 1, MPFR_RNDN);
            mpfr_mul(t, coef, cm, MPFR_RNDN);
            mpfr_sub(acc, acc, t, MPFR_RNDN);
        }
        mpfr_mul_ui(acc, acc, (j == 0 || j == n) ? 1 : 2, MPFR_RNDN);
        mpfr_div_si(c->w[l][j], acc, n, MPFR_RNDN);
    }
    mpfr_clears(pi, t, acc, coef, (mpfr_ptr)0);
    c->have = l + 1;
}

bool ccquad_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a, const mpfr_t b,
                 long bits, double reltol, mpfr_srcptr abstol, CcCache* cache,
                 mpfr_t out_re, mpfr_t out_im, double* abserr) {
    if (!cache || cache->bits != bits) return false;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    long nmax = (long)CC_N0 << CC_MAX_LEVEL;
    mpfr_t* fr = (mpfr_t*)malloc((size_t)(nmax + 1) * sizeof(mpfr_t));
    mpfr_t* fi = (mpfr_t*)malloc((size_t)(nmax + 1) * sizeof(mpfr_t));
    mpfr_t mid, half, x, Ir, Ii, Pr, Pi, dr, di, t;
    mpfr_inits2(p, mid, half, x, Ir, Ii, Pr, Pi, dr, di, t, (mpfr_ptr)0);
    mpfr_add(mid, a, b, MPFR_RNDN);  mpfr_div_2ui(mid, mid, 1, MPFR_RNDN);
    mpfr_sub(half, b, a, MPFR_RNDN); mpfr_div_2ui(half, half, 1, MPFR_RNDN);
    long have_n = 0;            /* f values held for the nodes of N = have_n */
    bool ok = true, conv = false;
    double lrel = log2(reltol);
    double labs = abstol ? de_mpfr_log2_abs(abstol) : -INFINITY;
    double lprevd = INFINITY;
    *abserr = INFINITY;
    for (int l = 0; l <= CC_MAX_LEVEL && ok; l++) {
        long n = (long)CC_N0 << l;
        if (cache->have <= l) cc_build(cache, l);
        /* gather f on the level's nodes: stored at index j * (nmax / n) */
        long stride = nmax / n;
        for (long j = 0; j <= n; j++) {
            if (have_n && (j & 1) == 0) continue;           /* nested node */
            long s = j * stride;
            mpfr_mul(x, cache->x[l][j], half, MPFR_RNDN);
            mpfr_add(x, mid, x, MPFR_RNDN);
            mpfr_init2(fr[s], p); mpfr_init2(fi[s], p);
            if (!f(ctx, x, fr[s], fi[s]) || !mpfr_number_p(fr[s]) || !mpfr_number_p(fi[s])) {
                mpfr_clear(fr[s]); mpfr_clear(fi[s]);
                /* clear what this level already stored, then stop */
                for (long q = 0; q < j; q++) {
                    if (have_n && (q & 1) == 0) continue;
                    mpfr_clear(fr[q * stride]); mpfr_clear(fi[q * stride]);
                }
                ok = false; break;
            }
        }
        if (!ok) break;
        have_n = n;
        mpfr_set_zero(Ir, 1); mpfr_set_zero(Ii, 1);
        for (long j = 0; j <= n; j++) {
            long s = j * stride;
            mpfr_mul(t, cache->w[l][j], fr[s], MPFR_RNDN); mpfr_add(Ir, Ir, t, MPFR_RNDN);
            mpfr_mul(t, cache->w[l][j], fi[s], MPFR_RNDN); mpfr_add(Ii, Ii, t, MPFR_RNDN);
        }
        mpfr_mul(Ir, Ir, half, MPFR_RNDN);
        mpfr_mul(Ii, Ii, half, MPFR_RNDN);
        if (l > 0) {
            mpfr_sub(dr, Ir, Pr, MPFR_RNDN);
            mpfr_sub(di, Ii, Pi, MPFR_RNDN);
            *abserr = de_mpfr_l1_d(dr, di);
            double ld = de_mpfr_log2_l1(dr, di);
            /* geometric convergence: doubling N squares the error, so once the
             * differences shrink the current error is ~ d_l^2 / d_{l-1} */
            double lest = (l >= 2 && ld < lprevd) ? 2.0 * ld - lprevd : ld;
            double lI = de_mpfr_log2_l1(Ir, Ii);
            /* rounding floor: the sum cannot settle below ~2^-bits |I| */
            if (ld == -INFINITY || lest <= lrel + lI || lest <= labs
                || ld <= lI - (double)bits + 8.0) {
                conv = true;
            }
            lprevd = ld;
        }
        mpfr_set(Pr, Ir, MPFR_RNDN);
        mpfr_set(Pi, Ii, MPFR_RNDN);
        if (conv) break;
    }
    if (have_n) {
        long stride = nmax / have_n;
        for (long j = 0; j <= have_n; j++) {
            mpfr_clear(fr[j * stride]); mpfr_clear(fi[j * stride]);
        }
    }
    if (conv) { mpfr_set(out_re, Ir, MPFR_RNDN); mpfr_set(out_im, Ii, MPFR_RNDN); }
    mpfr_clears(mid, half, x, Ir, Ii, Pr, Pi, dr, di, t, (mpfr_ptr)0);
    free(fr); free(fi);
    return conv;
}

#endif /* USE_MPFR */
