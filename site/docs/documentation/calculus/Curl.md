# Curl

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Curl[{f1, f2}, {x1, x2}]`**

gives the scalar curl D\[f2,x1\] - D\[f1,x2\].

**`Curl[{f1, f2, f3}, {x1, x2, x3}]`**

gives the vector curl (D\[f3,x2\]-D\[f2,x3\], D\[f1,x3\]-D\[f3,x1\], D\[f2,x1\]-D\[f1,x2\]).  For an n\*n\*...\*n array the generalized Levi-Civita curl (depth n-k-1) is returned.

**`Curl[f, {x1, ..., xn}, chart]`**

gives the curl of a vector field in the orthonormal basis of chart ("Cartesian", "Polar", "Cylindrical", "Spherical").

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

2-D curl is the scalar rotation D[f2,x] - D[f1,y]

```mathematica
In[9]:= Curl[{-y, x}, {x, y}]
Out[9]= 2
```

3-D curl of a rigid-rotation field

```mathematica
In[10]:= Curl[{-y, x, 0}, {x, y, z}]
Out[10]= {0, 0, 2}
```

Curl of a z-directed field

```mathematica
In[11]:= Curl[{0, 0, x^2 + y^2}, {x, y, z}]
Out[11]= {2 y, -2 x, 0}
```

## Implementation notes

**Algorithm.** `builtin_curl` (through the `vecop` front end) computes the
Cartesian curl as a generalized Levi-Civita contraction:
`(1/k!) Sum eps_{a..ij..} d_{x_i} f_{j..}`, where `k = ncube_depth(f)` is the
field's depth and the result has depth `n - k - 1`. `curl_perm_recur` enumerates
all `n!` permutations; `curl_emit_perm` splits each as `(a_part | i | j_part)`,
reads the leaf `f_{j_part}`, and accumulates `sign * D[f_{j_part}, x_i]`
(`perm_sign` is `(-1)^inversions`) into the cell indexed by `a_part`.
`curl_build_nested` then folds each cell's terms into `Plus[...]/k!` and
assembles the depth-`(n-k-1)` nested list, reduced with one `eval_and_free`.
This covers the three familiar cases uniformly: a 2-D vector gives a scalar, a
3-D vector gives a vector, and a rank-2 tensor gives a scalar.

**Data structures.** A `curl_ctx` holds the field, the variables, and a
per-cell array of growable term buffers (`curl_cell_append` doubles capacity);
`Expr` trees are built with `mk_d`/`mk_neg`/`mk_fnN_adopt`. There is no ND or
`Compile[]` path — the output is a symbolic derivative. The 3-argument
curvilinear form (`curl_chart`, `n = 2` or `3` only) builds the orthonormal-basis
curl from the chart's Lamé factors, e.g. the 3-D component
`(1/(h_j h_k))[D[h_k f_k, x_j] - D[h_j f_j, x_k]]` cyclically.

**Complexity / limits.** The permutation enumeration is `O(n!)`, so the
Cartesian path is bounded to `2 <= n <= 6`; outside that, or when the field is a
scalar (`k < 1`) or `k > n - 1` (negative result depth), it returns `NULL`. The
chart form declines (`NULL`) for `n` other than 2 or 3, a non-vector field, or
an unrecognised chart (emitting `Curl::chart`) — tensor curl in a curvilinear
basis needs a metric and is out of scope.

**Attributes:** `Protected`.

## References

**See also:** [Grad](../../calculus/Grad/), [Div](../../calculus/Div/), [Laplacian](../../calculus/Laplacian/), [D](../../calculus/D/)

- Source: [`src/vectoranal.c`](https://github.com/stblake/mathilda/blob/main/src/vectoranal.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_vectoranal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vectoranal.c)

## Notes & additional examples

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
