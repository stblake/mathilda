/*
 * cubature.h — adaptive multidimensional cubature over an axis-aligned box
 *              (machine precision)
 *
 * Estimates ∫_box f dV with the Genz-Malik degree-7 rule and its embedded
 * degree-5 rule for the local error estimate.  Regions are kept in a max-heap
 * keyed by error; the worst region is bisected along the axis of largest
 * fourth difference until the global error meets the tolerance or the
 * evaluation budget is exhausted.
 *
 * Unlike iterated 1-D quadrature, a single adaptive process drives the whole
 * box, so the cost grows polynomially (not multiplicatively) with dimension —
 * this is the engine NIntegrate uses for low-to-moderate dimensional integrals
 * over constant rectangular bounds.  The rule never samples the box corners
 * (all abscissae are strictly interior), so an integrable singularity sitting
 * on a box corner is handled without ever evaluating the singular point.
 *
 * The integrand is supplied through the same sample callback as the
 * Monte-Carlo engine (McSampleMachine): it takes the d-vector of coordinates
 * and returns false for a non-numeric sample.
 */
#ifndef MATHILDA_CUBATURE_H
#define MATHILDA_CUBATURE_H

#include <complex.h>
#include <stddef.h>

#include "mcint.h"   /* McSampleMachine */

typedef enum {
    CUB_OK        =  0,   /* converged to the requested tolerance        */
    CUB_NOCONV    =  1,   /* budget exhausted; *result is the best so far */
    CUB_NONNUMERIC = -1   /* a sample was non-numeric; *result invalid    */
} CubStatus;

/* Adaptive Genz-Malik cubature over the box [a,b] in d dimensions
 * (2 <= d <= 15).  On CUB_OK / CUB_NOCONV the estimate and its error are
 * written to *result / *abserr. */
CubStatus cub_integrate_machine(McSampleMachine f, void* ctx,
                                const double* a, const double* b, size_t d,
                                double abstol, double reltol, long max_eval,
                                double _Complex* result, double* abserr);

#endif /* MATHILDA_CUBATURE_H */
