#ifndef MATHILDA_TRIGRAT_H
#define MATHILDA_TRIGRAT_H

#include "expr.h"
#include "simp.h"

/*
 * trigrat.c -- Fast algebraic normal form for rational functions of
 * trigonometric and hyperbolic functions.
 *
 * Simplify's general path (simp_search) is a 2-round / ~14-transform
 * heuristic loop and is too slow on rational functions of Sin/Cos/Sinh/
 * Cosh (and Tan/Cot/Sec/Csc/Tanh/Coth/Sech/Csch after preprocessing).
 * simp_trig_rational implements a polynomial-time algorithm that works
 * in the quotient ring
 *
 *     K[ s_i, c_i, sigma_j, chi_j, other_vars ]
 *        / < s_i^2 + c_i^2 - 1,  chi_j^2 - sigma_j^2 - 1 >
 *
 * for the distinct trig arguments {a_i} (Sin/Cos) and hyperbolic
 * arguments {b_j} (Sinh/Cosh). It performs ideal reduction, denominator
 * rationalisation, multivariate GCD cancellation via the existing
 * Cancel builtin, and substitutes back. The result is returned only if
 * its leaf count is STRICTLY smaller than the input's, so the algorithm
 * cannot regress any case -- when it doesn't apply or doesn't improve,
 * simp_dispatch falls through to simp_search as before.
 *
 * Contract:
 *   - Input is borrowed (not freed). The classifier in simp_dispatch
 *     has already returned SIMP_SHAPE_TRIG for it.
 *   - On success: returns a freshly-allocated Expr*. The caller is
 *     responsible for freeing it.
 *   - On no-op / non-applicable / no-improvement: returns NULL. The
 *     caller continues with simp_search.
 */
Expr* simp_trig_rational(const Expr* input,
                         const AssumeCtx* ctx,
                         const Expr* complexity_func);

#endif /* MATHILDA_TRIGRAT_H */
