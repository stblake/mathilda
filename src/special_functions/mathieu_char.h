/* Mathilda -- MathieuCharacteristicA / MathieuCharacteristicB, the Mathieu
 * characteristic values a_r(q), b_r(q). See mathieu_char.c. */
#ifndef MATHILDA_MATHIEU_CHAR_H
#define MATHILDA_MATHIEU_CHAR_H

#include <stdbool.h>
#include "expr.h"

/* The machine helper shared by the scalar path, the NDArray / packed map and
 * Compile[]: a_r(q) (odd == 0) or b_r(q) (odd == 1) for real r and q. False
 * (declines) for b_0, a non-finite input, or |r|, |q| too large (~10^7). */
bool sf_mathieu_char_r(int odd, double r, double q, double* out);
/* Binary-kernel shapes of the above. */
bool sf_mathieu_char_a_r(double r, double q, double* out);
bool sf_mathieu_char_b_r(double r, double q, double* out);

Expr* builtin_mathieucharacteristica(Expr* res);
Expr* builtin_mathieucharacteristicb(Expr* res);
/* Register both heads and their attributes. Called from core_init(). */
void mathieu_char_init(void);

#endif /* MATHILDA_MATHIEU_CHAR_H */
