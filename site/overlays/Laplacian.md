### Worked examples

```mathematica
In[1]:= Laplacian[x^2 + y^2 + z^2, {x, y, z}]  (* Laplacian of the squared radius *)
```

```mathematica
In[1]:= Laplacian[f[x, y], {x, y}]  (* the sum of unmixed second partials *)
```

```mathematica
In[1]:= Laplacian[r^2, {r, t}, "Polar"]  (* Laplace-Beltrami form in a polar chart *)
```

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
