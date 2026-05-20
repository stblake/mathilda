#ifndef NULLSPACE_H
#define NULLSPACE_H

#include "expr.h"

/* NullSpace[m]                         -- basis for the null space of m
 * NullSpace[m, Method -> "<name>"]     -- explicit RREF method dispatch
 *
 * Returns a List of linearly independent vectors whose linear span is
 * { v : m . v == 0 }.  If the null space is trivial (i.e. m has full
 * column rank) the result is the empty List `{}`.
 *
 * Method names accepted (same grammar as RowReduce / LinearSolve /
 * Inverse via matsol_parse_method_option):
 *   "Automatic"                 -- alias for "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"  -- Bareiss-like fraction-free RREF
 *   "OneStepRowReduction"       -- classical Gauss-Jordan with division
 *   "CofactorExpansion"         -- identity-if-invertible via Det; for
 *                                  singular / rectangular input this
 *                                  falls back to DivFree internally,
 *                                  matching RowReduce's behaviour.
 *
 * For exact integer / rational input the null-space basis vectors are
 * scaled to clear integer denominators (each vector multiplied by the
 * LCM of its entries' integer denominators), so the result is
 * integer-valued whenever the input is integer-valued.  For symbolic
 * input the basis vectors are left in their natural rational form
 * (matching Mathematica's output convention).
 *
 * Vector ordering: basis vectors are returned with the rightmost free
 * column first, then the next-rightmost, and so on, matching
 * Mathematica.
 */
Expr* builtin_nullspace(Expr* res);
void  matnull_init(void);

#endif /* NULLSPACE_H */
