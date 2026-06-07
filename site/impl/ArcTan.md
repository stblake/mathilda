---
source: src/trig.c
---
**Algorithm.** `builtin_arctan` (`src/trig.c`) handles one- and two-argument forms. For `ArcTan[z]`: (1) `odd_fold` for oddness; (2) `trig_i_fold` for `ArcTan[I y] -> I ArcTanh[y]`; (3) exact inversion via `exact_arctan`, scanning n in `[-d/2, d/2]` for d in {1,2,3,4,5,6,10,12} against `exact_tan(n,d)` and returning `n/d * Pi`; (4) numeric fallback MPFR `mpfr_atan`/`mpfr_complex_atan`, else `get_approx` + C99 `catan`. For the two-argument `ArcTan[x, y]` (argument of x + I y): integer inputs are resolved exactly by quadrant (the axes and the four diagonals `±1` map to `0, ±Pi/2, Pi, ±Pi/4, ±3Pi/4`); real inputs use MPFR `mpfr_atan2` or C99 `atan2`. Indeterminate `ArcTan[0,0]` and unhandled cases return `NULL`.

**Data structures.** `Expr*` trees; the single-arg exact path reuses the forward `exact_tan` table.
