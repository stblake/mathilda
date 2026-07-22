---
source: src/trig.c
---
**Algorithm.** `builtin_arcsin` (`src/trig.c`): (1) `odd_fold` uses oddness `ArcSin[-x] -> -ArcSin[x]`. (2) `trig_i_fold` applies the principal-branch identity `ArcSin[I y] -> I ArcSinh[y]`. (3) Exact inversion via `exact_arcsin`, which brute-forces the forward table: for each denominator d in {1,2,3,4,5,6,10,12} and each n it computes `exact_sin(n,d)` and, if it `expr_eq`-matches the argument, returns `n/d * Pi`. (4) Numeric fallback: MPFR via `mpfr_asin`/`mpfr_complex_asin`, else `get_approx` + C99 `casin` for inexact inputs. On the real branch cut x>1 the imaginary part is negated to match the principal-branch lower-side convention (C99 lands on the upper side). Otherwise `NULL`.

**Data structures.** `Expr*` trees; the exact-inversion search reuses the forward `exact_sin` table rather than a separate inverse table.
