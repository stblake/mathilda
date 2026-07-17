/* cherry_ei.h — Cherry's rational exponential integral: ExpIntegralEi + Erf(i).
 *
 * Implements G. W. Cherry, "An Analysis of the Rational Exponential Integral"
 * (SIAM J. Comput. 1989) for the base field C(x): integrands gamma = g E^f with
 * g, f in C(x).  See CHERRY_PLAN.md §3-4.  Emits
 *   INT g E^f dx = y E^f + Sum_i c_i E^(-alpha_i) ExpIntegralEi[f + alpha_i]
 *                        + Sum_j k_j E^(-beta_j) Erfi[r_j/s]
 * where the ei arguments f + alpha_i come from the P1 resultant
 * Res_x(g1, p + alpha q) plus the P2 q-side term, and (when q = s^2 is a perfect
 * square) the erf arguments r_j/s from the completing-square beta finder
 * (r_j^2 = p + beta_j q).  y and all constants are solved as ONE linear system
 * over the rational matching identity
 *   g = y' + f' y + Sum_i c_i f'/(f + alpha_i) + Sum_j k_j (2/Sqrt[Pi]) (r_j/s)'.
 * Every candidate is diff-back verified, so a mis-generation can only decline,
 * never emit a wrong form.  Defined in cherry_ei.c.
 */

#ifndef MATHILDA_CHERRY_EI_H
#define MATHILDA_CHERRY_EI_H

#include "expr.h"

/* g E^f (g, f rational in x)  ->  y E^f + Sum c_i E^(-a_i) ExpIntegralEi[f+a_i].
 * Returns a fresh, diff-back-verified antiderivative, or NULL if the integrand
 * is not a base-field rational exponential of this form, if the ei argument
 * constants are not rational (the algebraic-constant layer is a later phase), or
 * if no constants solve the matching identity. */
Expr* rt_cherry_ei(Expr* f, Expr* x);

/* Cherry Thm 5.4 case b (flat exponential level): f rational in a single kernel E^w
 * whose Laurent expansion has several essential terms Sum_i p_i E^(i w) that the
 * single-shape rt_cherry_ei cannot peel together (e.g. (E^x + E^(2x))/(x-1)).
 * Integrate term-by-term (each p_i E^(i w) via rt_cherry_ei, the t^0 term rationally)
 * and sum; diff-back verified.  Tried after rt_cherry_ei declines; NULL otherwise. */
Expr* rt_cherry_exp_multiterm(Expr* f, Expr* x);

#endif /* MATHILDA_CHERRY_EI_H */
