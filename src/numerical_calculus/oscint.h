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

/* Integrate an oscillatory f over [a,b] (infinite=false) or (a, ∞)
 * (infinite=true, b ignored).  Writes *result and *abserr; returns true if the
 * panel sum / extrapolation settled to `reltol`. */
bool osc_integrate_machine(GkSampleMachine f, void* ctx, double a, double b,
                           bool infinite, double reltol, int max_panels,
                           double _Complex* result, double* abserr);

#endif /* MATHILDA_OSCINT_H */
