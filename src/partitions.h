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

/* PartitionsQ[n] — the number q(n) of partitions of the integer n into
 * distinct parts (equivalently, into odd parts; OEIS A000009).
 * builtin_partitionsq dispatches between two engines by size. */
Expr* builtin_partitionsq(Expr* res);

/* Engine 1: exact GMP recurrence from prod(1-x^k)prod(1+x^k)=prod(1-x^{2k}).
 * Computes q(n) for n >= 0 into `out`; returns 0 on success, -1 on OOM. */
int partitionsq_recurrence(unsigned long n, mpz_t out);

/* Engine 2: Hardy-Ramanujan-Rademacher / Hagis exact series (MPFR). Computes
 * q(n) for n >= 2 into `out`; returns 0 on success, -1 if it could not
 * converge (or MPFR is unavailable). Exposed for cross-validation. */
int partitionsq_hrr(unsigned long n, mpz_t out);

/* Register IntegerPartitions, PartitionsP and PartitionsQ with attributes. */
void partitions_init(void);

#endif /* PARTITIONS_H */
