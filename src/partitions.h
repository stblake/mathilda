#ifndef PARTITIONS_H
#define PARTITIONS_H

#include "expr.h"

/* IntegerPartitions[n], [n,kspec], [n,kspec,sspec], [n,kspec,sspec,m].
 * Enumerates the partitions of n in reverse-lexicographic order. See
 * docs/spec/builtins for the full surface. */
Expr* builtin_integerpartitions(Expr* res);

/* Register IntegerPartitions in the symbol table with its attributes. */
void partitions_init(void);

#endif /* PARTITIONS_H */
