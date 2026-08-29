#ifndef JORDANDECOMP_H
#define JORDANDECOMP_H

#include "expr.h"

/* JordanDecomposition[m] -- Jordan canonical form of a square matrix m.
 *
 * Returns a List {s, j} where j is block-diagonal in Jordan blocks and s is
 * the similarity matrix, so that m == s . j . Inverse[s].  Handles exact /
 * rational, free-symbolic, machine-precision (real / complex) and
 * arbitrary-precision (MPFR) matrices, and accepts packed / NDArray input.
 *
 * Standard builtin ownership contract: takes `res`, returns a new Expr* on
 * success (the evaluator frees `res`) or NULL to leave the call unevaluated
 * (does NOT free `res`).  See SPEC.md §4. */
Expr* builtin_jordandecomposition(Expr* res);

void  jordandecomp_init(void);

#endif /* JORDANDECOMP_H */
