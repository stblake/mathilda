/* nroots_internal.h — shared internal API for the NRoots engines.
 *
 * Defines the working complex-coefficient polynomial type (NrPoly) and the
 * three root-finding engines.  Included only by nroots.c, nroots_aberth.c and
 * nroots_jt.c.  Guarded by USE_MPFR: every numeric path computes in MPFR
 * complex arithmetic (the `ncpx` toolkit from numeric_complex.h), so without
 * MPFR there are no engines and NRoots returns unevaluated.
 *
 * Coefficient layout matches ZUPoly / poly_eval_mpfr:  c[i] = coeff of x^i.
 */
#ifndef MATHILDA_NROOTS_INTERNAL_H
#define MATHILDA_NROOTS_INTERNAL_H

#ifdef USE_MPFR

#include <mpfr.h>
#include "numeric_complex.h"   /* ncpx + ncpx_* */

/* Dense univariate polynomial with complex MPFR coefficients.
 * `c[i]` is the coefficient of x^i; `deg` is the highest index (deg >= 1).
 * `is_real` is true iff every coefficient's imaginary part is exactly zero
 * (lets engines pick real-arithmetic fast paths).  `prec` is the precision
 * (bits) every coefficient cell was initialised at — the working precision. */
typedef struct {
    int           deg;
    ncpx*         c;       /* length deg+1, each ncpx_init'd at `prec`        */
    int           is_real;
    mpfr_prec_t   prec;
} NrPoly;

/* Horner evaluation of p at z, with optional derivative.
 *   val : receives p(z)            (required, ncpx_init'd at >= wp)
 *   der : receives p'(z) or NULL   (optional)
 * `wp` is the scratch working precision for the inner products. */
void nr_poly_eval(const NrPoly* p, const ncpx* z, ncpx* val, ncpx* der,
                  mpfr_prec_t wp);

/* Newton polish of a single root estimate `z` against p, in place.  Runs up
 * to `max_iter` iterations at working precision `wp`; stops when the step is
 * negligible relative to |z|.  Safe on multiple roots (just stalls). */
void nr_newton_polish(const NrPoly* p, ncpx* z, mpfr_prec_t wp, int max_iter);

/* ---- Engines.  Each fills `roots[0..deg-1]` (caller pre-ncpx_init's every
 * cell at `wp`) with all roots counted with multiplicity, in arbitrary order.
 * Returns 0 on success, non-zero on failure (caller may escalate / warn). ---- */

/* Aberth–Ehrlich simultaneous iteration (default). */
int nr_aberth(const NrPoly* p, ncpx* roots, int max_iter, mpfr_prec_t wp);

/* Frobenius companion matrix eigenvalues (real QR; complex via 2n embedding).
 * Implemented in nroots.c. */
int nr_companion(const NrPoly* p, ncpx* roots, mpfr_prec_t wp);

/* Jenkins–Traub three-stage algorithm (RPOLY real / CPOLY complex). */
int nr_jenkinstraub(const NrPoly* p, ncpx* roots, int max_iter, mpfr_prec_t wp);

#endif /* USE_MPFR */

#endif /* MATHILDA_NROOTS_INTERNAL_H */
