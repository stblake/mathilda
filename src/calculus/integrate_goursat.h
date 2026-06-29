/* integrate_goursat.h
 *
 * Goursat's algorithm and its cube-/fourth-root generalisations for
 * pseudo-elliptic integrals -- integrands of the form
 *
 *     F(t) / R(t)^p ,   F rational in t,  R a polynomial in t,
 *
 * for p in {1/2, 1/3, 2/3, 1/4, 3/4}.  Such integrands live on a
 * positive-genus algebraic curve y^n = R(t) (n = den p) yet integrate
 * elementarily exactly when an eigenspace criterion holds under a Mobius
 * automorphism that cyclically permutes the roots of R:
 *
 *   p = 1/2  (R cubic or quartic, simple roots):  the Klein four-group V4
 *            of root-pairing involutions (Goursat 1887).  Elementary iff the
 *            trivial-character projection F^(0) vanishes.
 *   p = 1/2  (R cubic with the t^3-1 higher symmetry; Goursat 1887 Section 4):
 *            an order-3 Mobius S fixing one ramification point and cycling the
 *            other three.  Elementary iff F is a non-trivial period-3 character
 *            F(S) = alpha F (alpha = Exp[2 Pi I/3]), i.e. the period-3 trivial
 *            projection F + F(S) + F(S^2) vanishes.  Tried when V4 declines.
 *   p = 1/3, 2/3  (R cubic, or quadratic):  an order-3 Mobius cycle.
 *            Elementary iff the obstructive omega-eigencomponent vanishes
 *            (H1 at 1/3, H0 at 2/3).
 *   p = 1/4, 3/4  (R quartic, roots in harmonic position):  an order-4
 *            Mobius cycle.  Elementary iff the two obstructive eigencomponents
 *            vanish (V1,V2 at 1/4; V0,V1 at 3/4).
 *
 * When the criterion holds the integrand descends to genus-0 curves whose
 * Abelian integrals are rational; the method integrates those recursively and
 * back-substitutes.  Otherwise (obstructed / genuinely elliptic, or roots not
 * obtainable in radicals) it returns NULL so the Integrate cascade continues.
 *
 * Reference: S. Blake, "A Generalisation of Goursat's Algorithm for
 * Integration in Finite Terms" (cube root) and "A Fourth-Root Generalisation
 * of Goursat's Algorithm"; reference implementation GoursatAppendix.wl.
 */
#ifndef INTEGRATE_GOURSAT_H
#define INTEGRATE_GOURSAT_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.
 * Recognises F/R^p, runs the appropriate Goursat reduction, recurses into the
 * full integrator on the rational pieces, and back-substitutes.  Returns a
 * freshly owned antiderivative of `f` with respect to symbol `x`, or NULL when
 * `f` is not of pseudo-elliptic form, the integral is non-elementary
 * (obstructed), the radicand roots are not obtainable in radicals, or a
 * rational reduction does not close.  Does NOT take ownership of `f` or `x`. */
Expr* integrate_goursat_try(Expr* f, Expr* x);

/* `Integrate`GoursatAlgebraic[f, x]` builtin.  Same algorithm as the cascade
 * entry.  Strict: returns NULL (unevaluated) on any non-applicable input. */
Expr* builtin_integrate_goursat(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_goursat_init(void);

#endif /* INTEGRATE_GOURSAT_H */
