/* algebraicnumber.h — AlgebraicNumber[theta, {c0, c1, ..., cn}] builtin.
 *
 * Represents the element c0 + c1 theta + ... + cn theta^n of the rational
 * extension Q(theta). Canonicalised so theta is an algebraic integer and the
 * coefficient list has length equal to the degree of its minimal polynomial;
 * objects that represent a rational number collapse to explicit rational form.
 * The heavy lifting (qqbar conversion, algebraic-integer reduction, power-basis
 * reduction) lives in flint_qqbar.c; this file is the WL surface (arg checking,
 * fixpoint guard, registration). Declines (returns NULL) without FLINT.
 *
 * Ownership: builtin_algebraicnumber follows the standard contract — returns a
 * fresh Expr on success (the evaluator frees res) or NULL to leave res as is.
 */
#ifndef MATHILDA_ALGEBRAICNUMBER_H
#define MATHILDA_ALGEBRAICNUMBER_H

#include "expr.h"

Expr* builtin_algebraicnumber(Expr* res);
void algebraicnumber_init(void);

#endif /* MATHILDA_ALGEBRAICNUMBER_H */
