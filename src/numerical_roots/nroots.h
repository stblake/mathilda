/* nroots.h — NRoots[lhs == rhs, var, opts]
 *
 * Numerical root-finding for a univariate polynomial equation.  Returns a
 * disjunction of equations  var==r1 || var==r2 || ...  approximating every
 * root (repeating identical equations for roots of multiplicity > 1, per the
 * Wolfram Language).  Works at machine and arbitrary precision, for real and
 * complex coefficients, with three selectable methods: "Aberth" (default),
 * "CompanionMatrix", and "JenkinsTraub".
 *
 * Memory contract: never frees `res`; returns a fresh Expr* or NULL.
 */
#ifndef MATHILDA_NROOTS_H
#define MATHILDA_NROOTS_H

#include "expr.h"

Expr* builtin_nroots(Expr* res);
void  nroots_init(void);

#endif /* MATHILDA_NROOTS_H */
