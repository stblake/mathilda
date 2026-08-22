#ifndef FACINT_H
#define FACINT_H

#include "expr.h"
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>

Expr* builtin_primeq(Expr* res);
Expr* builtin_nextprime(Expr* res);
Expr* builtin_factorinteger(Expr* res);
Expr* builtin_eulerphi(Expr* res);

/* Complete prime factorisation of |n| for callers that need an exact,
 * proven-complete factor list (the Thue solver's maximal-order gate).
 * Fills ascending primes[]/exps[] (primes mpz_init'd on success; caller
 * mpz_clears), sets *count.  Returns true iff COMPLETE (every base a
 * strong probable prime and count <= cap); false => caller must DECLINE. */
bool facint_factor_complete(const mpz_t n, mpz_t* primes, int64_t* exps,
                            int cap, int* count);

void facint_init(void);

#endif
