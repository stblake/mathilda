# ListGradient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListGradient[f]`**

gives the numerical gradient of the sampled array f by finite differences: second-order central differences in the interior and one-sided differences at the edges, one array per axis. For a vector it returns a vector; for a rank-k array it returns {g1, ..., gk}.

**`ListGradient[f, h] uses spacing h on every axis; ListGradient[f, {s1,`**

<details>
<summary>Notes</summary>

..., sk}\] gives one spacing per axis, each a scalar spacing or a coordinate vector for non-uniform sampling. For a vector, a coordinate list of the same length as f is used directly. Options: Method -\> "Centered" (default), "Forward", or "Backward"; DifferenceOrder -\> p (accuracy order, default 2); WindowLength -\> m (stencil size, default p+1); Axis -\> All | a | {a1, ...} to restrict to particular axes. Integer input yields exact Rationals and symbolic input a symbolic result; arbitrary-precision (MPFR) input keeps its precision. Real and complex machine arrays use the packed/NDArray buffer fast path, and ListGradient\[v\] compiles. ListGradient has the attribute Protected.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ListGradient[{1, 4, 9, 16, 25}]
Out[1]= {3, 4, 6, 8, 9}

In[2]:= ListGradient[{a, b, c, d, e}]
Out[2]= {-a + b, -1/2 a + 1/2 c, -1/2 b + 1/2 d, -1/2 c + 1/2 e, -d + e}

In[3]:= ListGradient[{1, 4, 9, 16, 25}, {0, 1, 3, 6, 10}]
Out[3]= {3, 17/6, 73/30, 193/84, 9/4}
```

### Options (2)

```mathematica
In[4]:= ListGradient[{{1, 2, 6}, {3, 4, 5}}, Axis -> 2]
Out[4]= {{1, 5/2, 4}, {1, 1, 1}}

In[5]:= ListGradient[{0, 1, 8, 27, 64, 125, 216}, DifferenceOrder -> 4]
Out[5]= {0, 3, 12, 27, 48, 75, 108}
```

## Options & behaviour

> **Packed arrays.** `ListGradient` works straight on a `float64`/`float32`
> **or complex** (`complex64`/`complex32`) buffer (packed or visible `NDArray`),
> one axis at a time via a strided axpy-style pass — the finite-difference
> weights are real, so a complex element is just two contiguous doubles carried
> through the same pass. An **integer** buffer takes the ordinary path: the
> central difference of exact integers is a `Rational` (e.g.
> `(f[i+1]-f[i-1])/2`) that no buffer holds, so it materialises and the exact
> `List` path answers. Arbitrary-precision (`MPFR`) and exact
> (`Integer`/`Rational`) arrays likewise compute on the exact/symbolic `List`
> path, retaining full precision.

## Implementation notes

- `Protected`.
- Return shape follows numpy: a rank-1 `f` gives a single vector; a rank-`k` `f`
  gives `{g1, ..., gk}`, one same-shape array per axis (ascending axis order); a
  single requested axis gives one array, not a length-1 list.
- Default (`Method -> "Centered"`, `DifferenceOrder -> 2`) reproduces
  `numpy.gradient` exactly: an interior point is a second-order central
  difference and each extreme endpoint a first-order one-sided difference
  (numpy `edge_order = 1` — the edge accuracy is one order below the interior).
- **spacing** (numpy's `varargs`, as one positional argument):
  - absent → unit spacing on every axis;
  - a scalar (number or symbol) → uniform spacing on every axis;
  - a list of one entry per computed axis, each a scalar spacing or a coordinate
    vector of that axis's length (non-uniform, exact Fornberg weights). For a
    vector, a coordinate list of the same length as `f` is used directly.
- **Options**:
  - `Method -> "Centered"` (default), `"Forward"`, or `"Backward"` — the stencil
    placement. Forward/Backward stay order `p` at every point (shifting inward
    near the far edge); Centered drops the two extreme endpoints to order `p-1`.
  - `DifferenceOrder -> p` (accuracy order, default `2`).
  - `WindowLength -> m` (stencil size; default `Automatic` = `p + 1`). It is an
    alias for the order: an `m`-point stencil is order `m - 1`.
  - `Axis -> All` (default) `| a | {a1, ...}` — restrict to particular axes
    (1-based; negatives count from the end).
- The single kernel is Fornberg's finite-difference-weight recurrence, exact on
  non-uniform grids, so integer input gives exact `Rational`s and symbolic input
  a symbolic result. A machine-`Real` array uses the packed/NDArray fast path,
  and `ListGradient[v]` compiles (rank-1 array form) and auto-compiles.

**Attributes:** `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Real](../../other-advanced/Real/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_list_gradient.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list_gradient.c)
