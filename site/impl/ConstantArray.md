---
source: src/list/constant_array.c
---
**Algorithm.** `builtin_constant_array` builds a flat list of `n` copies of `c`
(`ConstantArray[c, n]`), or an `n1 x ... x nk` nested array
(`ConstantArray[c, {n1, ..., nk}]`). It is `Array[]` minus the index
computation: `ca_helper` recurses to the deepest level and returns a fresh
`expr_copy(c)` at each leaf, with no indexed function call built or evaluated.
Every dimension must be a non-negative machine integer (else the call is left
unevaluated); a `0` dimension yields an empty `List` at that level. The optional
third (padding) argument is accepted for arity but not implemented, so that form
is declined.

**Buffer fast path.** When `c` is a machine number (`EXPR_INTEGER` or
`EXPR_REAL`) over a rectangular shape, the whole result is known before anything
is built, so `ndbuild_open` opens a packed `NDArray` (`NDT_INT64` or
`NDT_FLOAT64`) and the single value is written straight into the buffer — no
per-element `Expr` is allocated. `UnitVector[10^6, ...]`-scale construction thus
costs a buffer fill rather than 10⁶ boxed nodes. `ndbuild_open` declines a zero
dimension or sub-threshold size, routing those shapes back to `ca_helper`.

**Data structures / limits.** `Expr**` children for the boxed path; a dense
`int64`/`double` buffer for the packed path (rank under `NDARRAY_MAX_RANK`). A
symbolic or compound `c` always takes the boxed path. `ATTR_PROTECTED`.
