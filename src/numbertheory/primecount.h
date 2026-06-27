/* primecount.h -- prime-counting function pi(x) with selectable algorithms.
 *
 * Internal header for the number-theory subsystem.  prime.c (Prime[n] /
 * PrimePi[x] builtins) is the only consumer.  Each algorithm is a recognised
 * method for computing pi(x); having several lets them cross-validate one
 * another (a trivially-correct sieve checks the combinatorial methods, which
 * check each other and the known values of pi(10^k)). */

#ifndef PRIMECOUNT_H
#define PRIMECOUNT_H

#include <stdint.h>

/* Largest x for which the sublinear methods (Lucy / LMO / Deleglise-Rivat)
 * are computed.  Bounds transient memory; larger x is reported out-of-range. */
#define PI_COUNT_MAX 50000000000000LL   /* 5 * 10^13 */

typedef enum {
    PC_AUTOMATIC = 0,  /* pick the best method for the magnitude of x */
    PC_LUCY,           /* Lucy_Hedgehog dynamic programming, O(x^3/4)   */
    PC_SIEVE,          /* segmented sieve of Eratosthenes, O(x log log x) */
    PC_LEGENDRE,       /* Legendre's formula                            */
    PC_MEISSEL,        /* Meissel's formula                             */
    PC_LEHMER,         /* Lehmer's formula                              */
    PC_LMO,            /* Lagarias-Miller-Odlyzko                       */
    PC_DR              /* Deleglise-Rivat                               */
} PrimeCountMethod;

/* pi(x): the number of primes <= x.  Returns -1 when x exceeds the range the
 * chosen method supports (the caller then leaves the call unevaluated).
 * x < 2 yields 0. */
int64_t prime_count(int64_t x, PrimeCountMethod method);

/* Small-prime table (primes < 10^6), built once and shared.  Used by Prime[n]
 * for the direct-index fast path. */
void     primecount_init(void);             /* build the table (idempotent) */
int64_t  primecount_small_table_size(void); /* count of primes in the table */
uint32_t primecount_small_prime(int64_t n); /* nth prime, 1-indexed (n in range) */

#endif
