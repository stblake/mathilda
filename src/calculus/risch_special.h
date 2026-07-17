/* risch_special.h — special-function outputs & base cases for the Risch integrator.
 *
 * The non-tower cases that emit higher transcendental functions or delegate to
 * the rational base: the pure-rational case (delegated to
 * Integrate`BronsteinRational), the Erf/Erfi, ExpIntegralEi, LogIntegral, and
 * dilogarithm (PolyLog) recognisers, and the logarithmic-polynomial case.
 * These are tried before the tower/exponential machinery.  Defined in
 * risch_special.c.
 */

#ifndef MATHILDA_RISCH_SPECIAL_H
#define MATHILDA_RISCH_SPECIAL_H

#include "expr.h"

/* Pure rational function of x -> Bronstein rational integrator (NULL if not
 * rational in x). */
Expr* rt_rational_case(Expr* f, Expr* x);

/* Top-monomial applicability mask for a special-function form (the C0 seam that
 * lets the extended-Liouville driver route a peeled tower coefficient to only the
 * relevant engines — CHERRY_PLAN.md §2.2).  A form declares which top monomials it
 * can close; an exponential top (E^w) admits ei/erf, a logarithmic top (Log[u])
 * admits li/dilog.  The outermost dispatch, where no tower has been built, passes
 * RT_SF_TOP_ANY to try every form (preserving the historical order and behaviour). */
#define RT_SF_TOP_EXP  0x1u   /* E^w  top  -> ExpIntegralEi / Erf   */
#define RT_SF_TOP_LOG  0x2u   /* Log[u] top -> LogIntegral / PolyLog */
#define RT_SF_TOP_ANY  (RT_SF_TOP_EXP | RT_SF_TOP_LOG)

/* Erf / Erfi / ExpIntegralEi / LogIntegral / dilogarithm recogniser (tries every
 * form; equivalent to rt_special_case_routed with RT_SF_TOP_ANY). */
Expr* rt_special_case(Expr* f, Expr* x);

/* As rt_special_case, but only try forms whose applicability mask intersects
 * `top_mask` (an RT_SF_TOP_* bitmask).  RT_SF_TOP_ANY reproduces rt_special_case
 * exactly (same forms, same order).  Used by the extended-Liouville tower hook to
 * route a peeled monomial coefficient by its top kind. */
Expr* rt_special_case_routed(Expr* f, Expr* x, unsigned top_mask);

/* Logarithmic-polynomial case (polynomial in a single Log[u]). */
Expr* rt_log_poly_case(Expr* f, Expr* x);

#endif /* MATHILDA_RISCH_SPECIAL_H */
