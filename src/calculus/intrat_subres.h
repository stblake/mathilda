/* intrat_subres.h — native subresultant PRS for the rational-function
 * integrator's Lazard-Rioboo-Trager log part.  See intrat_subres.c. */
#ifndef MATHILDA_INTRAT_SUBRES_H
#define MATHILDA_INTRAT_SUBRES_H

#include "expr.h"

/* The subresultant polynomial remainder sequence {R_0 = A, R_1 = B, R_2, ...,
 * R_k} of A and B in Z[t][x] (Bronstein, Symbolic Integration I, §1.5,
 * algorithm SubResultant), as an Expr List of expanded polynomials in x and t
 * with integer coefficients.  A and B are each scaled by a positive rational
 * constant to clear denominators, which scales every R_j by a constant and so
 * leaves the primitive parts and the x-degrees -- all the LRT log part reads
 * -- unchanged.
 *
 * Returns NULL (the caller keeps its generic path) unless both A and B are
 * polynomials in x and t with exact rational coefficients and deg_x A >=
 * deg_x B >= 1, or when an exact division of the chain fails (it cannot, for
 * well-formed input).  Borrows every argument; the caller owns the result. */
Expr* intrat_subres_chain(const Expr* A, const Expr* B, const Expr* x, const Expr* t);

/* The same subresultant chain with coefficients in Q[t, y1, ..., ym]: every
 * symbol of A and B other than x (parameters, the algebraic-mode generator z)
 * is a polynomial variable.  NULL when A or B is not a polynomial over Q in
 * symbols, past the size bounds, or FLINT is absent. */
Expr* intrat_subres_chain_mv(const Expr* A, const Expr* B, const Expr* x, const Expr* t);

#endif /* MATHILDA_INTRAT_SUBRES_H */
