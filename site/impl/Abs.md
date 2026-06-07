---
source: src/complex.c
---
`builtin_abs` handles each numeric kind directly: `mpz_abs` for `EXPR_BIGINT`, sign flips for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR` (via `mpfr_abs`), and `|n|/d` for rationals. For a `Complex[re, im]` literal (or an expression `complex_decompose` splits into numeric real/imag parts) it builds the symbolic modulus `Power[Plus[re^2, im^2], 1/2]`; when MPFR components are present it folds directly through `mpfr_hypot` at the combined working precision instead, which is also numerically stable across disparate magnitudes. Symbolic arguments return `NULL`.
