#ifndef SQUAREFREEQ_H
#define SQUAREFREEQ_H

#include "expr.h"

/* SquareFreeQ[expr]
 * SquareFreeQ[expr, vars]
 * SquareFreeQ[..., GaussianIntegers -> True|False|Automatic, Modulus -> p]
 *
 * Tests whether `expr` is square-free as an integer / Gaussian integer /
 * rational / polynomial.  Always returns True or False for any non-error
 * argument set; emits a Mathematica-compatible `SquareFreeQ::argb` or
 * `SquareFreeQ::nonopt` diagnostic and returns NULL (leaving the call
 * unevaluated) for argument errors. */
Expr* builtin_squarefreeq(Expr* res);

void squarefreeq_init(void);

#endif /* SQUAREFREEQ_H */
