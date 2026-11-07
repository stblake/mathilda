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
#include <stdbool.h>
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

/* As cub_integrate_machine, with two refinements.
 *
 * `reltol_deep` (<= reltol) is the relative tolerance once the refinement has
 * gone deep -- a region narrower than CUB_DEEP_REL of the box on some axis.
 * The embedded |deg 7 - deg 5| estimate is far above the true error of a
 * smooth integrand (it estimates the degree-5 rule's error), but tracks it
 * closely at an integrable singularity, where the bisection goes deep: there
 * a tolerance met is the error obtained, and the caller can ask for headroom
 * without paying for it on smooth integrands.
 *
 * When `hot` (length d) is non-NULL and the run did
 * not converge, it receives the point the refinement concentrated on: the
 * vertex of the smallest region nearest the mean centre of the
 * CUB_HOT_REGIONS smallest regions.  An integrable point singularity is where
 * the bisection piles up, and with it at a dyadic point of the box (the 0 of
 * [-Pi, Pi]) that vertex is the singular point exactly.  *hot_ok says whether
 * a point was written; hot_w (length d, or NULL) receives that region's
 * widths, the scale on which the point is known.
 *
 * `*forecast` true on entry (for a caller with a better fallback, the
 * iterated rule; NULL or false: off) ends the run early with CUB_NOCONV when
 * it cannot finish within max_eval, and says on return whether it did: at
 * each doubling of the samples once the refinement is deep (and past
 * CUB_FORECAST_MIN samples), the error's reduction exponent over the last two
 * doublings is extrapolated to the tolerance.  At an integrable vertex
 * singularity the Genz-Malik error falls only algebraically (x8 per doubling
 * in 2-D, x4 in 3-D), so a 12-digit goal would take millions of samples that
 * are then thrown away; a run on course to converge is not affected. */
#define CUB_HOT_REGIONS  8
#define CUB_DEEP_REL     0x1p-12
#define CUB_FORECAST_MIN 20000
CubStatus cub_integrate_machine_hot(McSampleMachine f, void* ctx,
                                    const double* a, const double* b, size_t d,
                                    double abstol, double reltol, double reltol_deep,
                                    long max_eval,
                                    double _Complex* result, double* abserr,
                                    double* hot, double* hot_w, bool* hot_ok,
                                    bool* forecast);

#endif /* MATHILDA_CUBATURE_H */
