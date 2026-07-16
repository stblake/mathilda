/* risch_singleext.h — flat single-extension transcendental Risch cases.
 *
 * The one-kernel integrators over C(x)(theta) for a single monomial theta =
 * Log[u] or E^u (u a rational function of x alone): the exponential
 * Laurent-polynomial case, the fractional (Rothstein–Trager / Hermite /
 * Lazard–Rioboo–Trager) cases, and the non-commensurate exponential-sum case
 * that decouples Sin/Cos·exp products.  These are dispatched before the
 * recursive tower engine (risch_field_integrate.c).  Defined in
 * risch_singleext.c.
 *
 * The Bronstein RdeBoundDegree helper rt_rde_var_bound is defined here but
 * declared in integrate_risch_transcendental.h (it is unit-tested white-box).
 */

#ifndef MATHILDA_RISCH_SINGLEEXT_H
#define MATHILDA_RISCH_SINGLEEXT_H

#include "expr.h"

/* Single-extension case entries (NULL == decline / out of scope). */
Expr* rt_exp_poly_case(Expr* f, Expr* x);   /* Laurent polynomial in E^u */
Expr* rt_frac_case(Expr* f, Expr* x);       /* proper fraction in Log[u] / E^u */
Expr* rt_hyperexp_case(Expr* f, Expr* x);   /* hyperexponential fraction */
Expr* rt_expsum_case(Expr* f, Expr* x);     /* sum of non-commensurate exponentials */
Expr* rt_hermite_case(Expr* f, Expr* x);    /* Hermite (repeated pole) reduction */

/* Kernelize a single-primitive exponential extension (E^(k u) -> rmT^k);
 * sets *u_out to the primitive exponent.  NULL if no single-primitive form. */
Expr* rt_exp_kernelize(Expr* f, Expr* x, Expr** u_out);

#endif /* MATHILDA_RISCH_SINGLEEXT_H */
