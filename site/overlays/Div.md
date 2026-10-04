### Worked examples

```mathematica
In[1]:= Div[{x, y, z}, {x, y, z}]  (* divergence of the radial field is the dimension *)
```

```mathematica
In[1]:= Div[{x^2, y^2}, {x, y}]  (* Cartesian divergence sums the diagonal partials *)
```

```mathematica
In[1]:= Div[{r^2, 0, 0}, {r, t, z}, "Cylindrical"]  (* divergence of a radial field in cylindrical coordinates *)
```

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
