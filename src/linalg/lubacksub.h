#ifndef LUBACKSUB_H
#define LUBACKSUB_H

#include "expr.h"

/* LUBackSubstitution[{lu, p, c}, b]: solve m . x == b from the triple
 * LUDecomposition[m] returns (b a vector or a matrix).  See lubacksub.c. */
Expr* builtin_lubacksubstitution(Expr* res);
void  lubacksub_init(void);

#endif /* LUBACKSUB_H */
