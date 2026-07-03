/* Mathilda — RootReduce: canonical (rationalised) form over an algebraic tower.
 * See src/rootreduce.c for the algorithm. */
#ifndef MATHILDA_ROOTREDUCE_H
#define MATHILDA_ROOTREDUCE_H

#include "expr.h"

/* RootReduce[expr]: rationalise radical / root-of-unity denominators and reduce
 * to a canonical representative over Q(params)(radicals). Returns expr unchanged
 * when it carries no algebraic generator or is out of the engine's scope. */
Expr* builtin_rootreduce(Expr* res);

/* Registers the RootReduce builtin. Called from core_init(). */
void rootreduce_init(void);

#endif /* MATHILDA_ROOTREDUCE_H */
