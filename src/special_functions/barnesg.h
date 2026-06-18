/* Mathilda -- BarnesG[z], the Barnes G-function.
 *
 * Exact for integer orders: G(n+1) = prod_{k=1}^{n-1} k! for positive n,
 * G(m) = 0 for non-positive integer m.  Non-integer orders stay symbolic.
 */
#ifndef BARNESG_H
#define BARNESG_H

#include "expr.h"

Expr* builtin_barnesg(Expr* res);
void barnesg_init(void);

#endif /* BARNESG_H */
