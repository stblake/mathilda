---
source: src/random.c
---
**Algorithm.** `builtin_randomcomplex` (in `src/random.c`) draws a point uniformly from the axis-aligned rectangle spanned by the corners. A bare `z` gives the rectangle `[0,Re z]×[0,Im z]`; `{z1, z2}` gives `[Re z1, Re z2]×[Im z1, Im z2]`. `random_complex_range` draws the real and imaginary parts as two independent uniforms via `random_uniform_01` (the same 53-bit-mantissa Mersenne Twister sampling used by `RandomReal`) and assembles a `Complex[...]`. An extended-precision path (`randomcomplex_mpfr`, guarded by `USE_MPFR`) draws two `mpfr_urandomb` deviates at the target precision. The `RandomComplex[range, n]` / `{n1,...}` forms produce lists or nested arrays via `random_complex_array`.
