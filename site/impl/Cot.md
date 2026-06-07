---
source: src/trig.c
---
**Algorithm.** `builtin_cot` (`src/trig.c`) runs the same cascade as the other trig heads. (1) `strip_inverse_call` folds `Cot[ArcCot[x]] -> x`. (2) `try_simp_forward_of_inverse` handles `Cot[ArcTan[x]] -> 1/x`. (3) `odd_fold` uses oddness `Cot[-x] -> -Cot[x]`. (4) `trig_i_fold` rewrites `Cot[I y] -> -I Coth[y]`. (5) `Cot[0] -> ComplexInfinity`. (6) For a rational multiple of Pi (via `extract_pi_multiplier`), `exact_cot` reduces n/d mod π into `[0, π/2]` with sign tracking and returns the table value (denominators 1,2,3,4,5,6,10,12; `ComplexInfinity` at multiples of π). (7) Numeric fallback: MPFR via `mpfr_cot`/`mpfr_complex_cot`, else `get_approx` + `1/ctan(c)`. Otherwise `NULL`.

**Data structures.** `Expr*` trees built with the `make_*` helpers; Pi multiples carried as `int64_t n, d`.
