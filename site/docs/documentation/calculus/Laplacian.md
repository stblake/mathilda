# Laplacian

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Laplacian[f, {x1, ..., xn}]`**

gives the Laplacian D\[f,{x1,2}\] + ... + D\[f,{xn,2}\]; for an array f the scalar Laplacian is applied to each component (same dimensions).

**`Laplacian[f, {x1, ..., xn}, chart]`**

gives the Laplacian of a scalar in chart ("Cartesian", "Polar", "Cylindrical", "Spherical").

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Grad[Sin[x^2 + y^2], {x, y}]
Out[1]= {2 x Cos[x^2 + y^2], 2 y Cos[x^2 + y^2]}

In[2]:= Grad[{x y, y z, z x}, {x, y, z}]
Out[2]= {{y, x, 0}, {0, z, y}, {z, 0, x}}

In[3]:= Div[{x^2, y^2, z^2}, {x, y, z}]
Out[3]= 2 x + 2 y + 2 z

In[4]:= Curl[{y, -x}, {x, y}]
Out[4]= -2

In[5]:= Laplacian[x^2 + y^2 + z^2, {x, y, z}]
Out[5]= 6

In[6]:= Div[{r Sin[t], -r Cos[t]}, {r, t}, "Polar"]
Out[6]= 3 Sin[t]

In[7]:= -Grad[k q/r, {r, t, p}, "Spherical"]
Out[7]= {(k q)/r^2, 0, 0}

In[8]:= Laplacian[Sin[r^2], {r, t}, "Polar"] // Simplify
Out[8]= 4 Cos[r^2] - 4 r^2 Sin[r^2]
```

### Applications (3)

Laplacian of the squared radius

```mathematica
In[9]:= Laplacian[x^2 + y^2 + z^2, {x, y, z}]
Out[9]= 6
```

The sum of unmixed second partials

```mathematica
In[10]:= Laplacian[f[x, y], {x, y}]
Out[10]= Derivative[0, 2][f][x, y] + Derivative[2, 0][f][x, y]
```

Laplace-Beltrami form in a polar chart

```mathematica
In[11]:= Laplacian[r^2, {r, t}, "Polar"]
Out[11]= 4
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [Grad](../../calculus/Grad/), [Div](../../calculus/Div/), [Curl](../../calculus/Curl/), [D](../../calculus/D/)

- Source: [`src/vectoranal.c`](https://github.com/stblake/mathilda/blob/main/src/vectoranal.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_vectoranal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vectoranal.c)

## Notes & additional examples

### Notes

`Laplacian[f, {x1, ..., xn}]` is `D[f, {x1, 2}] + ... + D[f, {xn, 2}]`, the sum
of unmixed second partial derivatives. Since `D` threads over an explicit array,
the scalar Laplacian is applied element-wise to a vector or tensor field, and the
result keeps the field's dimensions.

The three-argument form `Laplacian[f, vars, chart]` gives the Laplace–Beltrami
operator on a **scalar** in the orthonormal basis of a chart (`"Cartesian"`,
`"Polar"`, `"Cylindrical"`, `"Spherical"`), as
`(1/J) Sum_i D[(J/h_i^2) D[f, x_i], x_i]` with `J` the product of the Lamé
factors `h_i`. The vector Laplacian in a curvilinear chart needs Christoffel
symbols and is out of scope; such a call, or an unrecognised chart (which warns
`Laplacian::chart`), is returned unevaluated.
