/* tonumberfield.h — ToNumberField builtin.
 *
 * Argument forms (WL-faithful):
 *   ToNumberField[a, theta]        express a in the field Q(theta)
 *   ToNumberField[{a1,..}, theta]  express each ai in Q(theta)
 *   ToNumberField[{a1,..}]         express the ai in a common field (Automatic)
 *   ToNumberField[{a1,..}, All]    ... using the smallest common field
 *   ToNumberField[x]               convert x to an explicit AlgebraicNumber
 *
 * Results are AlgebraicNumber objects (or a plain rational when the field is Q).
 * The field/primitive-element math lives in flint_qqbar.c; this file is the WL
 * surface (dispatch + registration). Declines (NULL) without FLINT, and when an
 * argument is not a constant algebraic number or a does not lie in Q(theta).
 */
#ifndef MATHILDA_TONUMBERFIELD_H
#define MATHILDA_TONUMBERFIELD_H

#include "expr.h"

Expr* builtin_tonumberfield(Expr* res);
void tonumberfield_init(void);

#endif /* MATHILDA_TONUMBERFIELD_H */
