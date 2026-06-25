#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>

#include "expr.h"

/* Deterministic-stream support for callers (e.g. PossibleZeroQ's
 * Schwartz-Zippel sampler) that need a reproducible draw sequence without
 * disturbing the user's RandomInteger/SeedRandom stream. random_push_seed
 * saves the live PRNG state and reseeds from `seed`; random_pop_seed restores
 * the saved state. Calls nest (small fixed depth) and must be paired. */
void random_push_seed(uint64_t seed);
void random_pop_seed(void);

Expr* builtin_randominteger(Expr* res);
Expr* builtin_randomreal(Expr* res);
Expr* builtin_randomcomplex(Expr* res);
Expr* builtin_randomchoice(Expr* res);
Expr* builtin_randomsample(Expr* res);
Expr* builtin_seedrandom(Expr* res);

void random_init(void);

#endif // RANDOM_H
