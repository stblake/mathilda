#ifndef MATHILDA_LINALG_DOT_MPFR_H
#define MATHILDA_LINALG_DOT_MPFR_H

#include "expr.h"

/* a . b for arbitrary-precision vectors / matrices (rank 1 or 2), each output
 * element one fused, correctly rounded MPFR dot product.  One operand must be
 * all MPFR reals and the other MPFR reals or exact numbers (0 included); NULL (the
 * caller takes the generic path) for anything else or a shape mismatch.
 * Borrows both operands.  See dot_mpfr.c. */
Expr* dot_mpfr(const Expr* a, const Expr* b);

#endif
