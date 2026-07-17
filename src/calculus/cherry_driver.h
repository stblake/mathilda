/* cherry_driver.h — extended-Liouville dispatch for Cherry's special functions.
 *
 * The C0 seam of CHERRY_PLAN.md (§2.3): a single entry point that dispatches the
 * registered special-function forms (ExpIntegralEi / Erf / LogIntegral / PolyLog)
 * and aggregates their decision verdicts.  Both the outermost integrate dispatch
 * (Integrate`RischTranscendental) and the Thm 5.4 tower hook inside the field
 * recursion call this one function, so the Cherry engines fire uniformly on the
 * original integrand AND on peeled tower-monomial coefficients.
 *
 * Today the body is the historical rt_special_case loop (behaviour-preserving);
 * it grows the finite argument generators and the Sigma-decomposition NON-existence
 * aggregation without changing this signature.  Defined in cherry_driver.c.
 */

#ifndef MATHILDA_CHERRY_DRIVER_H
#define MATHILDA_CHERRY_DRIVER_H

#include "expr.h"

/* Solve the extended-Liouville problem  INT gamma dx = v + Sum_i k_i SF_i(a_i)
 * for the special functions whose top-monomial applicability intersects
 * `top_mask` (an RT_SF_TOP_* bitmask from risch_special.h; pass RT_SF_TOP_ANY at
 * the outermost, tower-free dispatch).  Returns a fresh, diff-back-verified
 * antiderivative in kernel form, or NULL when no form applies.  `f`, `x` borrowed. */
Expr* extended_liouville_solve(Expr* f, Expr* x, unsigned top_mask);

#endif /* MATHILDA_CHERRY_DRIVER_H */
