---
source: src/hyperbolic.c
---
**Algorithm.** `builtin_cosh` (`src/hyperbolic.c`) mirrors the trig cascade. (1) `strip_inverse_call` folds `Cosh[ArcCosh[x]] -> x`. (2) `try_simp_forward_of_inverse_hyp` handles `Cosh[ArcSinh[x]] -> Sqrt[1+x^2]` and `Cosh[ArcTanh[x]] -> 1/Sqrt[1-x^2]`. (3) `even_fold` for evenness `Cosh[-x] -> Cosh[x]`. (4) `hyp_i_fold` rewrites `Cosh[I y] -> Cos[y]`. (5) `Cosh[0] -> 1`; `Cosh[±Infinity] -> Infinity`. (6) Numeric fallback: MPFR via `mpfr_cosh`/`mpfr_complex_cosh`, else `get_approx` + C99 `ccosh` for inexact real/complex inputs. Otherwise `NULL`. There is no exact rational-multiple-of-Pi table for the hyperbolic heads.

**Data structures.** `Expr*` trees built with the shared `make_*` helpers.
