---
source: src/list/distance.c
---
**Algorithm.** `builtin_euclidean_distance` computes
`Sqrt[Sum_i Abs[u_i - v_i]^2]` via the shared `dist_builtin(res, p=2, root=true)`.
The core `dist_sum` forms `Sum_i Abs[u_i - v_i]^p` by composing the internal
arithmetic primitives (`internal_subtract`, an Abs helper, `internal_power`,
`internal_plus`) through `eval_and_free`, and `EuclideanDistance` then takes the
square root of that sum. Reusing the evaluator's arithmetic rather than
recomputing it gives three properties for free: exact input stays exact where the
result is rational; complex components use their modulus (`Abs`-then-square, as
Mathematica defines it); and symbolic input survives as a symbolic distance
rather than being rejected.

**Shape and the Abs gap.** `dist_shape` admits only two scalars, or two `List`s
of equal length (a length mismatch, or a list-valued component, declines — these
are vector functions, not matrix-threaded). The internal `dist_abs` decides a
real argument's sign with `list_numeric_sign` and negates if needed, routing
around a bug where `Abs` declines on a bigint-scale rational, which otherwise
left high-precision exact distances unevaluated; complex/symbolic arguments still
go through `internal_abs` for the modulus.

**Complexity / limits.** O(n) arithmetic evaluations over the `n` components;
interpreter-speed, not a buffer path. `ATTR_PROTECTED`. See also
`SquaredEuclideanDistance` (the root-free, exactness-preserving form that
`FindClusters` ranks on), `ManhattanDistance`, and `CosineDistance`.
