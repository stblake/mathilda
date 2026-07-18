/* Mathilda — vector-construction builtins. */
#ifndef VECTORS_H
#define VECTORS_H

#include "expr.h"

/* UnitVector[k]      -> the 2-D unit vector in the k-th direction.
 * UnitVector[n, k]   -> the n-D unit vector: a length-n list with a 1 in
 *                       position k and 0s elsewhere.
 *
 * Components are exact integers by default; the WorkingPrecision option
 * selects a machine-precision or MPFR real representation. */
Expr* builtin_unit_vector(Expr* res);

/* Register the vector builtins in the symbol table. */
void vectors_init(void);

#endif /* VECTORS_H */
