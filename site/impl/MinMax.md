---
source: src/list/minmax.c
---
**Algorithm.** `builtin_minmax` returns `{Min[arg], Max[arg]}`, built by evaluating a
`Min` call and a `Max` call on (refcount-shared copies of) the argument. Delegating
keeps every numeric subtlety — bignums, reals, symbolic extrema, the empty-list
`Infinity`/`-Infinity` — in exactly one place rather than reimplemented here. A
packed buffer, an `NDArray` and an association are all accepted; over an association
the values are used (as `Min`/`Max` already do).

**Data structures.** No bespoke storage: `expr_copy` is a refcount bump, so the same
buffer is handed to both `Min` and `Max`, and each takes its own buffer reduction
(`ndred_min` / `ndred_max`).

**Complexity / limits.** O(n) — effectively two single-pass reductions over the same
data. `MinMax` is on `pack.c`'s `AWARE` list, so a packed/`NDArray` argument stays on
the buffer; without that both halves would materialise 10⁶ boxed nodes (a measured
~735× regression against the buffer path). For a non-list/array/association argument
it declines and the call is left unevaluated.
