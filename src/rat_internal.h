/* rat_internal.h — a few rat.c helpers shared with the unified normalizer
 * (src/poly/ratcanon.c).  These were static in rat.c; exposed here so the
 * rewrite (RATCANON_REWRITE_PLAN.md) can reuse them without duplication.
 * Not part of the public builtin surface.
 */
#ifndef MATHILDA_RAT_INTERNAL_H
#define MATHILDA_RAT_INTERNAL_H

#include "expr.h"

/* Split `expr` into numerator and denominator (both freshly owned).  A negative
 * or fractional Power becomes a denominator factor; a Times distributes.  For a
 * non-fraction, num = copy(expr), den = 1. */
void extract_num_den(Expr* expr, Expr** num_out, Expr** den_out);

/* True when `e` reads as negative from its surface form (negative number, or a
 * Times/Rational whose leading factor is negative). */
bool is_superficially_negative(Expr* e);

/* -1 * e, evaluated (borrows e). */
Expr* negate_expr(Expr* e);

#endif /* MATHILDA_RAT_INTERNAL_H */
