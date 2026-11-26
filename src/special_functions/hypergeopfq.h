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
 * Hypergeometric0F1/1F1/2F1 evaluate through HypergeometricPFQ (every value,
 * closed form and numeric path is PFQ's) but keep their own head when PFQ
 * cannot reduce the call, as the reference does. Consumers that work in the
 * PFQ form (D, FunctionExpand, Simplify) convert with the helpers below.
 */
Expr* builtin_hypergeometric_pfq(Expr* res);
Expr* builtin_hypergeometric_0f1(Expr* res);
Expr* builtin_hypergeometric_1f1(Expr* res);
Expr* builtin_hypergeometric_2f1(Expr* res);

/* Hypergeometric0F1[b, z] / 1F1[a, b, z] / 2F1[a, b, c, z] as the equivalent
 * (unevaluated) HypergeometricPFQ call; NULL when e is not one of them. */
Expr* hyper_named_to_pfq(const Expr* e);

/* HypergeometricPFQ[{..}, {c}, z] with 0, 1 or 2 upper parameters as the named
 * head (unevaluated); NULL for any other shape. */
Expr* hyper_pfq_to_named(const Expr* e);

void hypergeopfq_init(void);

#endif /* HYPERGEOPFQ_H */
