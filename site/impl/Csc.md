---
source: src/trig.c
---
**Algorithm.** `builtin_csc` (`src/trig.c`) applies: (1) `strip_inverse_call` folds `Csc[ArcCsc[x]] -> x`. (2) `odd_fold` for oddness `Csc[-x] -> -Csc[x]`. (3) `trig_i_fold` rewrites `Csc[I y] -> -I Csch[y]`. (4) `Csc[0] -> ComplexInfinity`. (5) For a rational multiple of Pi (`extract_pi_multiplier`), `exact_csc` returns the closed-form surd from the table for denominators 1,2,3,4,5,6,10,12. (6) Numeric fallback: MPFR via `mpfr_csc`/`mpfr_complex_csc`, else `get_approx` + `1/csin(c)` for inexact inputs. Otherwise `NULL`. (Unlike Cos/Tan, Csc has no forward-of-inverse fold step.)

**Data structures.** `Expr*` trees via the `make_*` helpers; Pi multiples as `int64_t n, d`.
