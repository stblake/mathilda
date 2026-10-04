---
source: src/vectoranal.c
---
**Algorithm.** `builtin_grad` goes through the shared `vecop` front end, which
validates arity (2 or 3), takes a `normalized_copy` of the field `f` and the
coordinate list (materialising a packed `NDArray` to a nested `List` first), and
dispatches on argument count. The 2-argument Cartesian form is a direct
passthrough: `grad_cartesian` builds `D[f, {{x1,...,xn}}]` — the array-derivative
that appends a new innermost tensor slot — and reduces it with a single
`evaluate`. Nothing here re-implements differentiation; the interpreter's
array-`D` does the work, so a scalar becomes a vector, a vector becomes its
Jacobian, and a rank-`k` array gains one rank uniformly.

**Data structures.** Everything is an `Expr` tree assembled with the tiny
`mk_fn*` builders (`mk_d`, `mk_fn2`, `mk_fnN_adopt`) and collapsed by
`eval_and_free`. The 3-argument curvilinear form (`grad_chart`) accepts only a
scalar field and emits `{(1/h_i) D[f, x_i]}` in the orthonormal (physical) basis,
where `resolve_chart`/`chart_scale_factors` supply the Lamé scale factors `h_i`
for one of `"Cartesian"`, `"Polar"` (2-D), `"Cylindrical"` or `"Spherical"`
(3-D). This is a purely symbolic operator — no ND kernel and no `Compile[]`
lowering — since the result is a differentiated expression, not a machine buffer.

**Complexity / limits.** Cost is that of the underlying `D` plus one `evaluate`.
Following Mathilda's "can't evaluate" contract the builtin returns `NULL`
(leaving `Grad[...]` unevaluated) when the coordinate spec is not a list, when
the chart name is unrecognised (emitting `Grad::chart` through the message
funnel), or when a chart form is given a non-scalar field — the gradient of a
vector in a curvilinear chart needs Christoffel symbols, which are deliberately
out of scope.
