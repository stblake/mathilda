/* cherry_li.h — Cherry's logarithmic integral: LogIntegral (li) engine.
 *
 * Implements the single-logarithm case of G. W. Cherry, "Integration in Finite
 * Terms with Special Functions: The Logarithmic Integral" (SIAM J. Comput. 1986)
 * for a tower C(x, theta) with theta = Log[w], w a polynomial in x.  The answer
 *   INT gamma dx = v(x, theta) + Sum_k d_k LogIntegral[w^k]
 * has an elementary Laurent part v (polynomial in theta over C(x)) plus a sum of
 * logarithmic integrals whose arguments are POWERS w^k (the degree-1
 * Sigma-decomposition over the single generator w).  Matching is performed in the
 * TOWER — theta is an independent variable with Log[w^k] = k theta by
 * construction and D[LogIntegral[w^k]] = w^(k-1) w'/theta — so the exact rational
 * identity D_tower[answer] == gamma is itself the certificate (a plain Simplify
 * diff-back cannot reduce Log[w^k] = k Log[w]).  Generalises the single-term
 * rt_try_li recognizer to MULTI-li (Cherry d3: x^2/Log[x+1]).  Defined in
 * cherry_li.c.
 */

#ifndef MATHILDA_CHERRY_LI_H
#define MATHILDA_CHERRY_LI_H

#include "expr.h"

/* gamma over C(x, Log[w])  ->  v + Sum_k d_k LogIntegral[w^k], or NULL. */
Expr* rt_cherry_li(Expr* f, Expr* x);

/* The li NON-existence decision (Cherry 1986 Thm 5.4 case a / Thm 4.4).  True iff
 * f is the pure essential form A/Log[w] (w squarefree) whose reduced function
 * Phi = A/w' PROVABLY has no all-equal Sigma-decomposition over the irreducible
 * factors of w — i.e. INT f is not li-elementary (elementary + LogIntegral).
 * The decision property behind Integrate`LiElementaryQ; f, x are borrowed. */
bool rt_cherry_li_nonelem(Expr* f, Expr* x);

#endif /* MATHILDA_CHERRY_LI_H */
