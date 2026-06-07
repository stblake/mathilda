---
source: src/trig.c
---
**Algorithm.** `builtin_arccos` (`src/trig.c`): (1) `arc_pi_minus_fold` applies the reflection `ArcCos[-x] -> Pi - ArcCos[x]` for superficially-negative arguments. (2) `arccos_i_fold` applies the imaginary-axis identity `ArcCos[I y] -> Pi/2 - I ArcSinh[y]`. (3) Exact inversion via `exact_arccos`, which scans n in `[0,d]` for d in {1,2,3,4,5,6,10,12}, computes `exact_cos(n,d)`, and returns `n/d * Pi` on an `expr_eq` match. (4) Numeric fallback: MPFR via `mpfr_acos`/`mpfr_complex_acos`, else `get_approx` + C99 `cacos`; for real x>1 the imaginary part is negated to match Mathematica's `+i*acosh(x)` convention. Otherwise `NULL`.

**Data structures.** `Expr*` trees; exact inversion reuses the forward `exact_cos` table.
