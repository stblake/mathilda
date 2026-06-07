---
source: src/random.c
---
**Algorithm.** `builtin_randomreal` (in `src/random.c`) has two paths selected by the requested working precision. The machine path (`randomreal_machine`) draws a uniform `double` in `[0,1)` via `random_uniform_01`, which samples a 53-bit integer with `mpz_urandomm(big, g_rand_state, 2^53)` and divides by `2^53` — i.e. full-mantissa doubles from the shared **Mersenne Twister** state (`gmp_randinit_mt`). `random_real_range` affinely maps it to `[xmin, xmax)`. A bare `x` means `[0, x)`, `{a, b}` means `[a, b)`; bounds are coerced with `expr_to_real`.

When a precision argument requests extended precision, the MPFR path (`randomreal_mpfr` / `random_real_range_mpfr`, guarded by `USE_MPFR`) draws `mpfr_urandomb` at the target bit-width and affinely rescales with `MPFR_RNDN`, preserving exact rational bounds through `get_approx_mpfr`. The `RandomReal[range, n]` and `RandomReal[range, {n1,...}]` forms build lists/arrays via `random_real_array` (or its MPFR counterpart).
