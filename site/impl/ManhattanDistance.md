---
source: src/list/distance.c
---
**Algorithm.** `builtin_manhattan_distance` computes `Sum_i Abs[u_i - v_i]` via
the shared `dist_builtin(res, p=1, root=false)`. `dist_sum` forms the sum of
absolute component differences by composing `internal_subtract`, the `dist_abs`
helper and `internal_plus` through `eval_and_free`, so exact input stays exact
(`ManhattanDistance[{1, 2}, {4, 6}]` is `7`), complex components contribute their
modulus, and a symbolic pair comes back symbolically
(`ManhattanDistance[{a}, {b}]` is `Abs[a - b]`, as Mathematica answers).

**Shape.** `dist_shape` admits two scalars or two equal-length `List`s only; a
length mismatch, a list-valued component, or a non-atomic "scalar" (a `Rational`
is stored as `Rational[...]`, an `EXPR_FUNCTION`, so a bare rational scalar is
*not* accepted) declines. `dist_abs` decides a real sign with `list_numeric_sign`
to route around the `Abs` bigint-rational gap.

**Complexity / limits.** O(n) arithmetic evaluations over the components;
interpreter-speed. `ATTR_PROTECTED`. See `EuclideanDistance`,
`SquaredEuclideanDistance`, and `CosineDistance` — all share the `dist_sum` loop
and differ only in `p` and whether a root is taken.
