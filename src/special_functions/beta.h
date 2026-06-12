#ifndef BETA_H
#define BETA_H

#include "expr.h"

/* Beta[a, b]            -- Euler beta function  B(a, b) = Gamma(a)Gamma(b)/Gamma(a+b)
 * Beta[z, a, b]         -- incomplete beta      B_z(a, b) = Int_0^z t^(a-1)(1-t)^(b-1) dt
 * Beta[z0, z1, a, b]    -- generalized incomplete = B_{z1}(a, b) - B_{z0}(a, b)
 *
 * Attributes: Listable, NumericFunction, Protected, ReadProtected.
 *
 * Evaluation strategy (see beta.c for detail):
 *   - Beta[a, b] reduces through the existing Gamma builtin, so exact integer /
 *     half-integer / rational and machine / arbitrary-precision / complex
 *     arguments all inherit Gamma's reductions. Non-positive-integer poles map
 *     to ComplexInfinity; the doubly-singular lattice points reduce by the
 *     finite limit of the gamma ratio.
 *   - Beta[z, a, b] reduces through Hypergeometric2F1 via
 *     B_z(a, b) = z^a/a * 2F1(a, 1-b; a+1; z), which terminates to an exact
 *     closed form for positive-integer b and evaluates numerically otherwise.
 *   - Beta[z0, z1, a, b] = Beta[z1, a, b] - Beta[z0, a, b].
 *   - Everything else stays symbolic. */
Expr* builtin_beta(Expr* res);

void beta_init(void);

#endif /* BETA_H */
