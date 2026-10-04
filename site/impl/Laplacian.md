---
source: src/vectoranal.c
---
**Algorithm.** `builtin_laplacian` (through the `vecop` front end) builds the
Cartesian Laplacian as `Sum_i D[f, {x_i, 2}]` and reduces it with one
`eval_and_free`. Because `D` threads over an explicit array `f`, the result
carries `f`'s own dimensions — the scalar Laplacian is applied element-wise to a
vector or tensor field. Nothing re-implements differentiation; the sum of
unmixed second partials is handed straight to the interpreter's `D`.

**Data structures.** `Expr` trees via the `mk_d`/`mk_fnN_adopt` builders;
`normalized_copy` unpacks a packed `NDArray` field to a nested `List` first.
There is no ND kernel or `Compile[]` lowering, as the result is a symbolic
derivative. The 3-argument curvilinear form (`laplacian_chart`) accepts only a
scalar field and emits the Laplace–Beltrami form
`(1/J) Sum_i D[(J/h_i^2) D[f, x_i], x_i]`, with `J = Prod_i h_i` and the Lamé
scale factors `h_i` from `resolve_chart` for `"Cartesian"`, `"Polar"`,
`"Cylindrical"` or `"Spherical"`.

**Complexity / limits.** Cost is `n` second-derivative calls plus one
`evaluate`. It returns `NULL` (leaving `Laplacian[...]` unevaluated) when the
coordinate list is empty (`n < 1`), when a chart form is given a non-scalar
field — the vector Laplacian in a curvilinear basis needs Christoffel symbols
and is out of scope — or when the chart name is unrecognised (emitting
`Laplacian::chart`).
