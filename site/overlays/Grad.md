### Worked examples

```mathematica
In[1]:= Grad[x^2 + y^2, {x, y}]  (* gradient of a scalar field *)
```

```mathematica
In[1]:= Grad[f[x, y], {x, y}]  (* an unknown scalar stays symbolic via Derivative *)
```

```mathematica
In[1]:= Grad[{x y, y z, z x}, {x, y, z}]  (* a vector field yields its Jacobian -- one rank higher *)
```

```mathematica
In[1]:= Grad[r^2, {r, t}, "Polar"]  (* orthonormal gradient in a polar chart *)
```

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
