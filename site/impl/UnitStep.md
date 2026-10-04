---
source: src/piecewise.c
---
**Algorithm.** `builtin_unitstep` classifies each argument's sign with
`ustep_class` (a numerical-certification test that decides `< 0`, `>= 0`, or
unresolved). `UnitStep[]` is 1; any argument certified negative makes the whole
call 0; arguments certified non-negative contribute a factor of 1 and are
dropped; if every argument is non-negative the result is 1, and if some remain
unresolved the call returns `UnitStep` over just those (returning `NULL`,
unevaluated, when nothing could be resolved). The result is always the exact
integer 0 or 1 when fully determined.

**Data structures.** Plain `Expr` arguments plus a small `int` class array; the
reduced call is rebuilt with `expr_new_function`. The ND kernel is **narrowing**:
`ndk_UnitStep_i` takes a `double` to an `int64` 0/1 and `ndk_UnitStep_ii` takes an
`int64` to an `int64`, with no real-closed or complex arm on purpose (the answer
is always an integer). It is on `pack.c`'s AWARE + `INT64_OK` list, so a packed or
visible numeric `NDArray` narrows to an integer buffer rather than materialising
boxed reals.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it to `CT_INT` at
scalar and rank-1 shapes (a `double` input still yields an integer: `UnitStep[0.5]`
is `1`, not `1.`). The sign test is `UnitStep[0] = 1` (the step is closed at zero,
following `Sign[-0.] = 0`). A non-real or sign-undecidable argument stays symbolic.
