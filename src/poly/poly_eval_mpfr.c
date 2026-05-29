/* poly_eval_mpfr.c — see poly_eval_mpfr.h for the API. */

#include "poly_eval_mpfr.h"

#ifdef USE_MPFR

#include <stdlib.h>

void poly_eval_real_mpfr(const mpfr_t* c, int deg, const mpfr_t x,
                         mpfr_t out, mpfr_t* dout) {
    if (deg < 0) {
        mpfr_set_zero(out, 1);
        if (dout) mpfr_set_zero(*dout, 1);
        return;
    }
    mpfr_prec_t bits = mpfr_get_prec(out);
    mpfr_t v, dv, tmp;
    mpfr_init2(v, bits);
    mpfr_init2(dv, bits);
    mpfr_init2(tmp, bits);
    mpfr_set(v, c[deg], MPFR_RNDN);
    mpfr_set_zero(dv, 1);
    for (int i = deg - 1; i >= 0; i--) {
        /* dv = dv * x + v */
        mpfr_mul(tmp, dv, x, MPFR_RNDN);
        mpfr_add(dv, tmp, v, MPFR_RNDN);
        /* v = v * x + c[i] */
        mpfr_mul(tmp, v, x, MPFR_RNDN);
        mpfr_add(v, tmp, c[i], MPFR_RNDN);
    }
    mpfr_set(out, v, MPFR_RNDN);
    if (dout) mpfr_set(*dout, dv, MPFR_RNDN);
    mpfr_clear(v); mpfr_clear(dv); mpfr_clear(tmp);
}

void poly_eval_complex_mpfr(const mpfr_t* c, int deg,
                            const mpfr_t xr, const mpfr_t xi,
                            mpfr_t out_r, mpfr_t out_i,
                            mpfr_t* dr, mpfr_t* di) {
    if (deg < 0) {
        mpfr_set_zero(out_r, 1); mpfr_set_zero(out_i, 1);
        if (dr) mpfr_set_zero(*dr, 1);
        if (di) mpfr_set_zero(*di, 1);
        return;
    }
    mpfr_prec_t bits = mpfr_get_prec(out_r);
    mpfr_t vr, vi, dvr, dvi, t1, t2;
    mpfr_init2(vr, bits);  mpfr_init2(vi, bits);
    mpfr_init2(dvr, bits); mpfr_init2(dvi, bits);
    mpfr_init2(t1, bits);  mpfr_init2(t2, bits);

    mpfr_set(vr, c[deg], MPFR_RNDN);
    mpfr_set_zero(vi, 1);
    mpfr_set_zero(dvr, 1);
    mpfr_set_zero(dvi, 1);

    for (int i = deg - 1; i >= 0; i--) {
        /* (dvr + i*dvi) = (dvr + i*dvi) * (xr + i*xi) + (vr + i*vi)
         * Real part: dvr*xr - dvi*xi + vr
         * Imag part: dvr*xi + dvi*xr + vi                          */
        mpfr_mul(t1, dvr, xr, MPFR_RNDN);
        mpfr_mul(t2, dvi, xi, MPFR_RNDN);
        mpfr_sub(t1, t1, t2, MPFR_RNDN);   /* dvr*xr - dvi*xi      */
        mpfr_add(t1, t1, vr, MPFR_RNDN);   /* + vr                 */
        mpfr_mul(t2, dvr, xi, MPFR_RNDN);
        mpfr_mul(dvr, dvi, xr, MPFR_RNDN); /* reuse dvr as scratch */
        mpfr_add(dvi, t2, dvr, MPFR_RNDN); /* dvr*xi + dvi*xr      */
        mpfr_add(dvi, dvi, vi, MPFR_RNDN); /* + vi                 */
        mpfr_set(dvr, t1, MPFR_RNDN);

        /* (vr + i*vi) = (vr + i*vi) * (xr + i*xi) + c[i]
         * Real: vr*xr - vi*xi + c[i]
         * Imag: vr*xi + vi*xr                                      */
        mpfr_mul(t1, vr, xr, MPFR_RNDN);
        mpfr_mul(t2, vi, xi, MPFR_RNDN);
        mpfr_sub(t1, t1, t2, MPFR_RNDN);
        mpfr_add(t1, t1, c[i], MPFR_RNDN);
        mpfr_mul(t2, vr, xi, MPFR_RNDN);
        mpfr_mul(vr, vi, xr, MPFR_RNDN);   /* reuse vr as scratch  */
        mpfr_add(vi, t2, vr, MPFR_RNDN);
        mpfr_set(vr, t1, MPFR_RNDN);
    }
    mpfr_set(out_r, vr, MPFR_RNDN);
    mpfr_set(out_i, vi, MPFR_RNDN);
    if (dr) mpfr_set(*dr, dvr, MPFR_RNDN);
    if (di) mpfr_set(*di, dvi, MPFR_RNDN);
    mpfr_clear(vr);  mpfr_clear(vi);
    mpfr_clear(dvr); mpfr_clear(dvi);
    mpfr_clear(t1);  mpfr_clear(t2);
}

void zupoly_to_mpfr_coeffs(const ZUPoly* p, mpfr_prec_t bits,
                           mpfr_t** out, int* deg_out) {
    int d = p->deg;
    *deg_out = d;
    if (d < 0) { *out = NULL; return; }
    mpfr_t* arr = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(d + 1));
    for (int i = 0; i <= d; i++) {
        mpfr_init2(arr[i], bits);
        mpfr_set_z(arr[i], p->c[i], MPFR_RNDN);
    }
    *out = arr;
}

void poly_eval_mpfr_free_coeffs(mpfr_t* coeffs, int deg) {
    if (!coeffs || deg < 0) return;
    for (int i = 0; i <= deg; i++) mpfr_clear(coeffs[i]);
    free(coeffs);
}

#endif /* USE_MPFR */
