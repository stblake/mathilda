/* Mathilda -- Hyperfactorial[n] = prod_{k=1}^{n} k^k.
 *
 * Exact for non-negative integer orders (GMP); other orders stay symbolic.
 */
#ifndef HYPERFACTORIAL_H
#define HYPERFACTORIAL_H

#include "expr.h"

Expr* builtin_hyperfactorial(Expr* res);
void hyperfactorial_init(void);

#endif /* HYPERFACTORIAL_H */
