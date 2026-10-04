---
source: src/vectoranal.c
---
**Algorithm.** `builtin_div` goes through the shared `vecop` front end (arity 2
or 3; `normalized_copy` of the field and the coordinate list, unpacking any
`NDArray`). The 2-argument Cartesian form, `build_div`, contracts the innermost
slot of `f` against the variables: a length-`n` vector becomes
`Sum_i D[f_i, x_i]`, while a tensor (a list of lists) maps `Div` recursively over
its outer structure, yielding a rank-`(k-1)` result. The assembled `Plus`/`List`
expression is reduced with one `eval_and_free`; a scalar, or a vector whose
length does not match `n` (or that is ragged), returns `NULL`.

**Data structures.** `Expr` trees built with the `mk_d`/`mk_fnN_adopt` helpers;
no ND or `Compile[]` wiring, as the output is a symbolic derivative. The
3-argument curvilinear form (`div_chart`) accepts only a length-`n` vector and
emits the orthonormal-basis divergence `(1/J) Sum_i D[(J/h_i) f_i, x_i]`, where
`J = Prod_i h_i` (`mk_jacobian`) and the Lamé factors `h_i` come from
`resolve_chart` for `"Cartesian"`, `"Polar"`, `"Cylindrical"` or `"Spherical"`.

**Complexity / limits.** Dominated by the inner `D` calls and one `evaluate`.
It returns `NULL` (leaving `Div[...]` unevaluated) for a scalar argument, a
shape/length mismatch, an unrecognised chart (emitting `Div::chart`), or a
tensor field in a chart — tensor divergence in a curvilinear basis needs a
metric connection and is out of scope.
