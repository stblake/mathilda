#ifndef MATSOL_H
#define MATSOL_H

#include "expr.h"

/* RowReduce[m]                                  -- default Method -> Automatic.
 * RowReduce[m, Method -> "<name>"]              -- explicit dispatcher.
 *
 * LinearSolve[m, b]                             -- default Method -> Automatic.
 * LinearSolve[m, b, Method -> "<name>"]         -- explicit dispatcher.
 *
 * Supported method names:
 *   "Automatic"                  -- aliased to "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"   -- Bareiss-like fraction-free Gauss-Jordan
 *   "OneStepRowReduction"        -- classical Gauss-Jordan with division
 *   "CofactorExpansion"          -- Cramer's rule (LinearSolve); identity-if-
 *                                   invertible (RowReduce), with fallback to
 *                                   DivisionFreeRowReduction on
 *                                   singular / rectangular input.
 *
 * Method -> Automatic (the symbol) is also accepted, matching Mathematica.
 */
Expr* builtin_rowreduce(Expr* res);
Expr* builtin_linearsolve(Expr* res);
void  matsol_init(void);

#endif /* MATSOL_H */
