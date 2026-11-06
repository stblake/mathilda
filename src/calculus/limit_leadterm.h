/*
 * limit_leadterm.h -- Limit at a finite point from the leading term of the
 * expression's expansion there, computed structurally (products, integer
 * powers, compositions F(G(x)) with G vanishing at the point) and by a
 * truncated Series only where a sum must be resolved.  See limit_leadterm.c.
 */
#ifndef LIMIT_LEADTERM_H
#define LIMIT_LEADTERM_H

#include "expr.h"

/* Limit[f, x -> x0] for a finite x0 free of x, when f = c (x - x0)^p (1 + o(1))
 * with a provably nonzero c and p >= 0: c for p = 0, 0 for p > 0.  Owned
 * result, or NULL to decline (no leading term found, p < 0, or a coefficient
 * whose zero test is undecided).  Arguments borrowed. */
Expr* limit_leading_term_value(Expr* f, Expr* x, Expr* x0);

#endif /* LIMIT_LEADTERM_H */
