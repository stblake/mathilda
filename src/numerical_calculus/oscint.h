/*
 * oscint.h — quadrature for oscillatory integrands (machine precision)
 *
 * For an integrand that oscillates many times across the integration range
 * neither global-adaptive Gauss-Kronrod (which would need a subdivision per
 * half-period) nor the double-exponential rules (which assume monotone decay
 * toward the endpoints) are effective.  This module integrates between the
 * successive zeros of the oscillation and combines the per-half-period panels:
 *
 *   - finite [a,b]: the panels are summed directly (a large but finite count);
 *   - half line (a, ∞): the partial sums over successive panels form an
 *     alternating-style sequence that is accelerated to its limit with Wynn's
 *     epsilon algorithm (shared seqaccel kernel) — the standard "integration
 *     between the zeros + nonlinear sequence acceleration" method for
 *     ∫_a^∞ of a decaying oscillatory function (Bessel, Sin[x]/x, …).
 *
 * Each half-period panel is integrated with the existing adaptive Gauss-Kronrod
 * rule, so smooth-but-oscillatory and mildly singular panels are both handled.
 * This is the engine behind Method -> "LevinRule" and the Automatic fallback
 * when the primary rule fails to converge on an oscillatory integrand.
 */
#ifndef MATHILDA_OSCINT_H
#define MATHILDA_OSCINT_H

#include <complex.h>
#include <stdbool.h>

#include "gkadapt.h"   /* GkSampleMachine */

#ifdef USE_MPFR
#  include <mpfr.h>
#  include "dequad.h"  /* DeQuadSampleMPFR */
#endif

/* Integrate an oscillatory f over [a,b] (infinite=false) or (a, ∞)
 * (infinite=true, b ignored).  Writes *result and *abserr; returns true if the
 * panel sum / extrapolation settled to `reltol`. */
bool osc_integrate_machine(GkSampleMachine f, void* ctx, double a, double b,
                           bool infinite, double reltol, int max_panels,
                           double _Complex* result, double* abserr);

#ifdef USE_MPFR
/* Arbitrary-precision counterpart of the half-line (a, ∞) branch of
 * osc_integrate_machine.  The oscillation geometry (half-period, successive
 * zeros) is located with the cheap machine sampler `fgeo` — a lobe boundary
 * needs only a handful of correct digits — while each half-period panel is
 * integrated at full `bits` precision with tanh-sinh (`fmp`) and the partial
 * sums are extrapolated with the MPFR Wynn epsilon.  (out_re,out_im) must be
 * pre-initialised by the caller.  *abserr receives the final extrapolation-step
 * magnitude as a double.  Returns whether the extrapolation settled to
 * `reltol`. */
bool osc_integrate_mpfr(GkSampleMachine fgeo, void* geoctx,
                        DeQuadSampleMPFR fmp, void* mpctx,
                        double a, long bits, double reltol,
                        int max_panels, int panel_levels,
                        mpfr_t out_re, mpfr_t out_im, double* abserr);
#endif

#endif /* MATHILDA_OSCINT_H */
