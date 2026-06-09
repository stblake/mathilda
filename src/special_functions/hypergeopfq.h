#ifndef HYPERGEOPFQ_H
#define HYPERGEOPFQ_H

#include "expr.h"

/*
 * Mathilda — generalized hypergeometric function HypergeometricPFQ.
 *
 *   HypergeometricPFQ[{a1,...,ap}, {b1,...,bq}, z]
 *     = Sum_{k>=0} ( prod_i Pochhammer(a_i,k) / prod_j Pochhammer(b_j,k) )
 *                  * z^k / k!
 *
 * The builtin evaluates:
 *   - structural rules (z==0 -> 1, common-parameter cancellation, threading
 *     over a list third argument);
 *   - termination to an exact polynomial when an upper parameter is a
 *     non-positive integer (valid for symbolic z);
 *   - machine-, arbitrary-precision (MPFR), and complex numeric values via
 *     direct series summation in the convergent regime.
 *
 * Hypergeometric0F1/1F1/2F1 are thin convenience heads that rewrite to
 * HypergeometricPFQ.
 */
Expr* builtin_hypergeometric_pfq(Expr* res);
Expr* builtin_hypergeometric_0f1(Expr* res);
Expr* builtin_hypergeometric_1f1(Expr* res);
Expr* builtin_hypergeometric_2f1(Expr* res);

void hypergeopfq_init(void);

#endif /* HYPERGEOPFQ_H */
