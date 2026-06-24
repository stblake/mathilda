#ifndef PARTITIONS_H
#define PARTITIONS_H

#include "expr.h"
#include <gmp.h>

/* IntegerPartitions[n], [n,kspec], [n,kspec,sspec], [n,kspec,sspec,m].
 * Enumerates the partitions of n in reverse-lexicographic order. See
 * docs/spec/builtins for the full surface. */
Expr* builtin_integerpartitions(Expr* res);

/* PartitionsP[n] — the number p(n) of unrestricted partitions of the
 * integer n. builtin_partitionsp dispatches between two engines by size. */
Expr* builtin_partitionsp(Expr* res);

/* Engine 1: Euler pentagonal-number recurrence (exact GMP integers).
 * Computes p(n) for n >= 0 into `out`; returns 0 on success, -1 on OOM. */
int partitionsp_recurrence(unsigned long n, mpz_t out);

/* Engine 2: Hardy-Ramanujan-Rademacher exact formula (MPFR). Computes p(n)
 * for n >= 2 into `out`; returns 0 on success, -1 if it could not converge
 * (or MPFR is unavailable). Exposed for cross-validation against engine 1. */
int partitionsp_hrr(unsigned long n, mpz_t out);

/* Register IntegerPartitions and PartitionsP with their attributes. */
void partitions_init(void);

#endif /* PARTITIONS_H */
