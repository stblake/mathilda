---
source: src/random.c
---
`builtin_seedrandom` (in `src/random.c`) reseeds the single global GMP **Mersenne Twister** state `g_rand_state` shared by all `Random*` builtins. `SeedRandom[n]` calls `gmp_randseed_ui` for a machine integer or `gmp_randseed` for a bignum seed, making subsequent random draws reproducible; `SeedRandom[]` reseeds from system entropy (`time(NULL) ^ clock()`). Both first call `ensure_rand_init` to lazily construct the state with `gmp_randinit_mt`. Returns `Null`; a non-integer seed returns `NULL` (unevaluated).
