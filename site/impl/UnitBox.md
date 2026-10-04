---
source: src/piecewise.c
---
**Algorithm.** `builtin_unitbox` reuses `ustep_class` twice rather than adding a
second classifier: x is in the box iff neither shifted argument `x + 1/2` nor
`1/2 - x` is certified negative. So `UnitBox[x]` is 1 when both one-sided
`UnitStep`-shaped tests pass, 0 when either shifted argument is negative, and
unevaluated when a side cannot be decided. The box is **closed** at both ends
(`UnitBox[1/2] = 1`), matching `UnitStep[0] = 1`, and the result is always the
exact integer 0 or 1 when determined.

**Data structures.** Each element costs two `Expr` allocations and two
`evaluate()` calls (the shifted arguments `x ± 1/2`), traded for reusing
`ustep_class`'s certification logic instead of duplicating it. The ND kernel
(`ndk_UnitBox_i` double→int64, `ndk_UnitBox_ii` int64→int64, `REG_U`) is
narrowing and lives on `pack.c`'s AWARE + `INT64_OK` list. Unlike `UnitStep`,
`Ramp`, `Round` and `IntegerPart`, `UnitBox` does **not** thread over an
`Interval` argument — the interval machinery encloses only monotone functions and
a two-sided box is not one.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers the single-argument
box to `CT_INT` at scalar and rank-1 shapes (`UnitBox[0.5]` is `1`, not `1.`). A
non-real or undecidable argument is left unevaluated.
