# Grad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Grad[f, {x1, ..., xn}]`**

gives the gradient (D\[f,x1\], ..., D\[f,xn\]) of the scalar f; for an array f a new innermost slot is appended, so a vector field yields its Jacobian.  Equivalent to D\[f, {{x1, ..., xn}}\].

**`Grad[f, {x1, ..., xn}, chart]`**

gives the gradient of a scalar in the orthonormal basis of chart, one of "Cartesian", "Polar", "Cylindrical", "Spherical".

## Examples (12)

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

### Applications (4)

Gradient of a scalar field

```mathematica
In[9]:= Grad[x^2 + y^2, {x, y}]
Out[9]= {2 x, 2 y}
```

An unknown scalar stays symbolic via Derivative

```mathematica
In[10]:= Grad[f[x, y], {x, y}]
Out[10]= {Derivative[1, 0][f][x, y], Derivative[0, 1][f][x, y]}
```

A vector field yields its Jacobian -- one rank higher

```mathematica
In[11]:= Grad[{x y, y z, z x}, {x, y, z}]
Out[11]= {{y, x, 0}, {0, z, y}, {z, 0, x}}
```

Orthonormal gradient in a polar chart

```mathematica
In[12]:= Grad[r^2, {r, t}, "Polar"]
Out[12]= {2 r, 0}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [Div](../../calculus/Div/), [Curl](../../calculus/Curl/), [Laplacian](../../calculus/Laplacian/), [D](../../calculus/D/)

- Source: [`src/vectoranal.c`](https://github.com/stblake/mathilda/blob/main/src/vectoranal.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_vectoranal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vectoranal.c)

## Notes & additional examples

### Notes

`Grad[f, {x1, ..., xn}]` is exactly `D[f, {{x1, ..., xn}}]`: the array-derivative
that appends a new innermost tensor slot. So a scalar becomes a vector, a vector
field becomes its Jacobian, and a rank-`k` array gains one rank uniformly.

The three-argument form `Grad[f, vars, chart]` gives the gradient of a **scalar**
in the orthonormal (physical) basis of an orthogonal coordinate chart, built from
the chart's Lamé scale factors. Supported charts are `"Cartesian"`, `"Polar"`
(2-D), `"Cylindrical"` and `"Spherical"` (3-D). The gradient of a vector field in
a curvilinear chart needs Christoffel symbols and is out of scope: such a call,
an unrecognised chart (which warns `Grad::chart`), or a non-list coordinate spec
is returned unevaluated — Mathilda's honest "can't evaluate" contract.
