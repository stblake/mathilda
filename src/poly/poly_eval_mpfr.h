/* poly_eval_mpfr.h
 * ----------------
 * High-precision Horner evaluation of a dense univariate polynomial at
 * a real or complex MPFR point, with optional simultaneous derivative
 * computation. Built for N[Root[..]] but kept polynomial-generic.
 *
 * Coefficient layout is `c[i] = coefficient of x^i` — the same layout
 * used by `ZUPoly` and produced by `zupoly_to_mpfr_coeffs` below.
 *
 * Memory: callers own all mpfr_t cells; this module never allocates
 * (except `zupoly_to_mpfr_coeffs`, which allocates and mpfr_init2's the
 * full coefficient array — release with `poly_eval_mpfr_free_coeffs`).
 */

#ifndef POLY_EVAL_MPFR_H
#define POLY_EVAL_MPFR_H

#ifdef USE_MPFR

#include <stddef.h>
#include <mpfr.h>

#include "zupoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Evaluate p(x) at a real MPFR point.
 *
 *   c, deg : coefficient array with c[i] = coeff of x^i; deg = highest
 *            non-zero exponent (deg >= 0).
 *   x      : evaluation point (read-only).
 *   out    : caller-allocated mpfr_t, receives p(x).
 *   dout   : optional. If non-NULL, receives p'(x); same precision as `out`.
 *
 * All mpfr_t cells are operated at the precision of `out`. The function
 * runs Horner once and, if dout is given, accumulates the derivative in
 * the same sweep (fused (p, p')). */
void poly_eval_real_mpfr(const mpfr_t* c, int deg, const mpfr_t x,
                         mpfr_t out, mpfr_t* dout);

/* Evaluate p(x) at a complex MPFR point z = xr + i*xi.
 *
 *   c, deg : real coefficients (c[i] = coeff of x^i).
 *   xr, xi : real, imag parts of z.
 *   out_r, out_i : receive p(z) = out_r + i*out_i.
 *   dr, di       : optional. If non-NULL, receive p'(z) = dr + i*di.
 *
 * All output cells are operated at the precision of `out_r`. */
void poly_eval_complex_mpfr(const mpfr_t* c, int deg,
                            const mpfr_t xr, const mpfr_t xi,
                            mpfr_t out_r, mpfr_t out_i,
                            mpfr_t* dr, mpfr_t* di);

/* Materialise a ZUPoly's coefficients as a freshly allocated dense MPFR
 * array of length deg+1, each cell mpfr_init2'd to `bits`. *deg_out is
 * the polynomial's degree (>= 0; the caller must check that p is non-zero
 * before calling).
 *
 * Allocation must be released with `poly_eval_mpfr_free_coeffs`. */
void zupoly_to_mpfr_coeffs(const ZUPoly* p, mpfr_prec_t bits,
                           mpfr_t** out, int* deg_out);

/* Free an array returned by zupoly_to_mpfr_coeffs.  Safe to pass NULL. */
void poly_eval_mpfr_free_coeffs(mpfr_t* coeffs, int deg);

#ifdef __cplusplus
}
#endif

#endif /* USE_MPFR */

#endif /* POLY_EVAL_MPFR_H */
