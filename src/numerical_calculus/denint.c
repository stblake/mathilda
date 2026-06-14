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
    bool have_prev = false;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++, h *= 0.5) {
        double _Complex sum = 0.0, c0;
        if (term(NULL, tctx, 0.0, &c0)) sum = c0;
        double maxterm = cabs(sum);
        for (int side = -1; side <= 1; side += 2) {
            double prev_m = INFINITY;
            bool decaying = false;
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > tmax) break;
                double _Complex tm;
                if (!term(NULL, tctx, t, &tm)) break;
                double m = cabs(tm);
                if (decaying && k > 3 && m > 4.0 * prev_m) break;  /* tail noise */
                sum += tm;
                if (m > maxterm) maxterm = m;
                if (m < 0.1 * maxterm) decaying = true;
                if (k > 2 && m <= DEN_REL_FLOOR * maxterm) break;
                prev_m = m;
            }
        }
        double _Complex Ival = h * sum;
        best = Ival;
        if (have_prev) {
            double diff = cabs(Ival - prev);
            *abserr = diff;
            if (diff <= reltol * cabs(Ival) + 1e-300) { *result = Ival; return true; }
        }
        prev = Ival;
        have_prev = true;
    }
    *result = best;
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
static double den_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}

/* tanh-sinh MPFR term at parameter tv. */
static bool den_ts_term_mpfr(DeQuadSampleMPFR f, void* ctx, DenMpfrScratch* sc,
                             double tv, mpfr_t term_re, mpfr_t term_im) {
    mpfr_set_d(sc->t, tv, MPFR_RNDN);
    mpfr_sinh(sc->s, sc->t, MPFR_RNDN);
    mpfr_mul(sc->u, sc->s, sc->halfpi, MPFR_RNDN);        /* u = (π/2) sinh t */
    mpfr_abs(sc->au, sc->u, MPFR_RNDN);
    if (mpfr_get_d(sc->au, MPFR_RNDN) > 1e9) return false;
    mpfr_mul_2ui(sc->tmp, sc->au, 1, MPFR_RNDN);          /* 2|u| */
    mpfr_neg(sc->tmp, sc->tmp, MPFR_RNDN);
    mpfr_exp(sc->eu, sc->tmp, MPFR_RNDN);                 /* eu = exp(-2|u|) */
    mpfr_add_ui(sc->onep, sc->eu, 1, MPFR_RNDN);          /* 1 + eu */
    /* weight = (π/2) cosh t · sech^2(u) = (π/2) cosh t · 4 eu/(1+eu)^2 */
    mpfr_cosh(sc->tmp, sc->t, MPFR_RNDN);
    mpfr_mul(sc->weight, sc->halfpi, sc->tmp, MPFR_RNDN);
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

static bool den_run_mpfr(DeQuadSampleMPFR f, void* ctx, DenMpfrScratch* sc,
                         bool tanhsinh, long bits, double tmax,
                         double reltol, int max_levels,
                         mpfr_t out_re, mpfr_t out_im, double* abserr) {
    if (reltol <= 0.0) reltol = 1e-13;
    if (max_levels < 1) max_levels = 1;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    double floor_rel = ldexp(1.0, -(int)bits - 4);

    mpfr_t sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di;
    mpfr_inits2(p, sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di, (mpfr_ptr)0);

    double h = 1.0;
    bool have_prev = false, converged = false;
    *abserr = INFINITY;

    for (int level = 0; level <= max_levels; level++, h *= 0.5) {
        mpfr_set_ui(sum_re, 0, MPFR_RNDN);
        mpfr_set_ui(sum_im, 0, MPFR_RNDN);
        double maxterm = 0.0;
        bool ok0 = tanhsinh ? den_ts_term_mpfr(f, ctx, sc, 0.0, tr, ti)
                            : den_ss_term_mpfr(f, ctx, sc, 0.0, tr, ti);
        if (ok0) {
            mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
            mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
            maxterm = den_l1_d(tr, ti);
        }
        for (int side = -1; side <= 1; side += 2) {
            double prev_m = INFINITY;
            bool decaying = false;
            for (int k = 1; ; k++) {
                double t = side * k * h;
                if (fabs(t) > tmax) break;
                bool ok = tanhsinh ? den_ts_term_mpfr(f, ctx, sc, t, tr, ti)
                                   : den_ss_term_mpfr(f, ctx, sc, t, tr, ti);
                if (!ok) break;
                double m = den_l1_d(tr, ti);
                if (decaying && k > 3 && m > 4.0 * prev_m) break;
                mpfr_add(sum_re, sum_re, tr, MPFR_RNDN);
                mpfr_add(sum_im, sum_im, ti, MPFR_RNDN);
                if (m > maxterm) maxterm = m;
                if (m < 0.1 * maxterm) decaying = true;
                if (k > 2 && m <= floor_rel * maxterm) break;
                prev_m = m;
            }
        }
        mpfr_mul_d(I_re, sum_re, h, MPFR_RNDN);
        mpfr_mul_d(I_im, sum_im, h, MPFR_RNDN);
        mpfr_set(out_re, I_re, MPFR_RNDN);
        mpfr_set(out_im, I_im, MPFR_RNDN);
        if (have_prev) {
            mpfr_sub(dr, I_re, prev_re, MPFR_RNDN);
            mpfr_sub(di, I_im, prev_im, MPFR_RNDN);
            double diff = den_l1_d(dr, di), mag = den_l1_d(I_re, I_im);
            *abserr = diff;
            if (diff <= reltol * mag + 1e-300) { converged = true; break; }
        }
        mpfr_set(prev_re, I_re, MPFR_RNDN);
        mpfr_set(prev_im, I_im, MPFR_RNDN);
        have_prev = true;
    }

    mpfr_clears(sum_re, sum_im, I_re, I_im, prev_re, prev_im, tr, ti, dr, di, (mpfr_ptr)0);
    return converged;
}

bool denint_tanhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          const mpfr_t a, const mpfr_t b,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr) {
    DenMpfrScratch sc; den_scratch_init(&sc, (mpfr_prec_t)bits);
    mpfr_set(sc.a, a, MPFR_RNDN);
    mpfr_set(sc.b, b, MPFR_RNDN);
    mpfr_sub(sc.half, b, a, MPFR_RNDN);
    mpfr_div_2ui(sc.half, sc.half, 1, MPFR_RNDN);   /* (b-a)/2 */
    double tmax = 4.0 + 0.25 * log((double)bits);
    bool ok = den_run_mpfr(f, ctx, &sc, true, bits, tmax, reltol, max_levels,
                           out_re, out_im, abserr);
    den_scratch_clear(&sc);
    return ok;
}

bool denint_sinhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr) {
    DenMpfrScratch sc; den_scratch_init(&sc, (mpfr_prec_t)bits);
    double tmax = 4.5 + 0.25 * log((double)bits);
    bool ok = den_run_mpfr(f, ctx, &sc, false, bits, tmax, reltol, max_levels,
                           out_re, out_im, abserr);
    den_scratch_clear(&sc);
    return ok;
}

#endif /* USE_MPFR */
