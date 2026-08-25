---
source: src/core.c
---
`builtin_mod` handles the two- and three-argument forms. For Integer and BigInt operands it reduces with GMP's `mpz_fdiv_r` — a *floored* remainder, so the result carries the sign of the divisor `n`; a rational path computes `m - n*Floor[m/n]` exactly in `mpq`; MPFR operands compute `m - n*floor(m/n)` at the larger of the two input precisions; and a machine-`double` fallback covers mixed Integer/Real arguments. The three-argument `Mod[m, n, d]` reduces `m - d` and adds `d` back, landing the representative in the half-open range `[d, d+n)`. `Mod` is registered `PROTECTED | NUMERICFUNCTION | LISTABLE` and runs element-wise on a packed machine buffer. Non-numeric arguments return NULL (left symbolic).
