#ifndef GATHER_H
#define GATHER_H

#include "expr.h"

/* Gather[list]: collect structurally identical elements of list into sublists,
 * in order of each element's first occurrence. Equivalent to
 * GatherBy[list, Identity]. */
Expr* builtin_gather(Expr* res);

#endif /* GATHER_H */
