# Div

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Div[{f1, ..., fn}, {x1, ..., xn}]`**

gives the divergence D\[f1,x1\] + ... + D\[fn,xn\]; for a rank-k array f it contracts the innermost slot against the variables, yielding a rank-(k-1) result.

**`Div[f, {x1, ..., xn}, chart]`**

gives the divergence of a vector field in the orthonormal basis of chart ("Cartesian", "Polar", "Cylindrical", "Spherical").

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

Divergence of the radial field is the dimension

```mathematica
In[9]:= Div[{x, y, z}, {x, y, z}]
Out[9]= 3
```

Cartesian divergence sums the diagonal partials

```mathematica
In[10]:= Div[{x^2, y^2}, {x, y}]
Out[10]= 2 x + 2 y
```

Divergence of a radial field in cylindrical coordinates

```mathematica
In[11]:= Div[{r^2, 0, 0}, {r, t, z}, "Cylindrical"]
Out[11]= 3 r
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [Grad](../../calculus/Grad/), [Curl](../../calculus/Curl/), [Laplacian](../../calculus/Laplacian/), [D](../../calculus/D/)

- Source: [`src/vectoranal.c`](https://github.com/stblake/mathilda/blob/main/src/vectoranal.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_vectoranal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vectoranal.c)

## Notes & additional examples

### Notes

`Div` contracts the innermost slot of its field against the variables: a
length-`n` vector gives the scalar `D[f1, x1] + ... + D[fn, xn]`, while a
rank-`k` tensor (a list of lists) maps `Div` over its outer structure to give a
rank-`(k-1)` result.

A scalar has no divergence, and a vector whose length does not match the number
of variables is a shape error; both are returned unevaluated. The three-argument
form `Div[f, vars, chart]` gives the divergence of a vector field in the
orthonormal basis of a chart (`"Cartesian"`, `"Polar"`, `"Cylindrical"`,
`"Spherical"`) as `(1/J) Sum_i D[(J/h_i) f_i, x_i]` with `J` the product of the
Lamé factors. Tensor divergence in a curvilinear chart needs a metric connection
and is out of scope; an unrecognised chart warns `Div::chart`.
