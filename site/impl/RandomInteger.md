---
source: src/random.c
---
**Algorithm.** `builtin_randominteger` (in `src/random.c`) draws uniform integers from a single global GMP random state, `g_rand_state`, lazily initialized by `ensure_rand_init` as a **Mersenne Twister** (`gmp_randinit_mt`) seeded from `time(NULL) ^ clock()`. A range is parsed by `parse_range` into `mpz_t` bounds: a bare `n` means `[0, n]`, `{a, b}` means `[a, b]`. `random_integer_range` computes `range = b - a + 1` and draws `mpz_urandomm(result, g_rand_state, range)` (rejection-free uniform over `[0, range)`), then adds `a`. The result is normalized via `expr_bigint_normalize`, so it demotes to `EXPR_INTEGER` when it fits and stays `EXPR_BIGINT` otherwise — arbitrarily large ranges are supported.

The `RandomInteger[range, n]` and `RandomInteger[range, {n1, n2, ...}]` forms produce a list or nested array via `random_array`, which recurses over the dimension spec drawing one element per leaf.
