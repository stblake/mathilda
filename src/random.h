#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>

#include "expr.h"

/* Deterministic-stream support for callers (e.g. PossibleZeroQ's
 * Schwartz-Zippel sampler) that need a reproducible draw sequence without
 * disturbing the user's RandomInteger/SeedRandom stream. random_push_seed
 * saves the live PRNG state and reseeds from `seed`; random_pop_seed restores
 * the saved state. Calls nest (small fixed depth) and must be paired. */
/* Uniform in [lo, hi] from the GMP stream, for INTERNAL samplers only.
 *
 * Kept separate from RandomInteger on purpose: the zero-test sampler's points
 * are part of a decision procedure (seeded from the expression's hash so a
 * verdict is a pure function of its input), not a user-visible sequence.
 * Speeding up the user-facing generator must not move them -- doing so changed
 * an Integrate result. See the definition in random.c. */
int64_t random_internal_int_range(int64_t lo, int64_t hi);

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
