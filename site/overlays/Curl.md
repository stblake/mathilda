### Worked examples

```mathematica
In[1]:= Curl[{-y, x}, {x, y}]  (* 2-D curl is the scalar rotation D[f2,x] - D[f1,y] *)
```

```mathematica
In[1]:= Curl[{-y, x, 0}, {x, y, z}]  (* 3-D curl of a rigid-rotation field *)
```

```mathematica
In[1]:= Curl[{0, 0, x^2 + y^2}, {x, y, z}]  (* curl of a z-directed field *)
```

### Notes

The curl is a generalized Levi-Civita contraction
`(1/k!) Sum eps_{a..ij..} d_{x_i} f_{j..}`, where `k` is the depth of the field
`f` and the result has depth `n - k - 1`. Concretely: a 2-D vector gives a scalar
`D[f2, x1] - D[f1, x2]`, a 3-D vector gives the usual vector curl, and a rank-2
tensor gives a scalar.

Because the implementation enumerates all `n!` index permutations, the Cartesian
form is bounded to `2 <= n <= 6`. A scalar field, or a field whose depth exceeds
`n - 1`, has no curl and is returned unevaluated. The three-argument form
`Curl[f, vars, chart]` (dimension 2 or 3 only) gives the orthonormal-basis curl
built from the chart's Lamé factors (`"Cartesian"`, `"Polar"`, `"Cylindrical"`,
`"Spherical"`); an unrecognised chart warns `Curl::chart`.
