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

## Examples (12)

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

### Applications (7)

Central differences inside, one-sided at the ends

```mathematica
In[6]:= ListGradient[{1, 2, 4, 7, 11}]
Out[6]= {1, 3/2, 5/2, 7/2, 4}
```

Symbolic input gives a symbolic gradient

```mathematica
In[7]:= ListGradient[{a, b, c, d}]
Out[7]= {-a + b, -1/2 a + 1/2 c, -1/2 b + 1/2 d, -c + d}
```

Machine input rides the float buffer

```mathematica
In[8]:= ListGradient[{1., 2., 4., 7., 11.}]
Out[8]= {1.0, 1.5, 2.5, 3.5, 4.0}
```

Uniform spacing h = 2

```mathematica
In[9]:= ListGradient[{1, 2, 4, 7, 11}, 2]
Out[9]= {1/2, 3/4, 5/4, 7/4, 2}
```

A non-uniform coordinate grid

```mathematica
In[10]:= ListGradient[{0, 1, 4, 9, 16}, {0, 1, 2, 4, 8}]
Out[10]= {1, 2, 17/6, 9/4, 7/4}
```

Rank 2: one array per axis

```mathematica
In[11]:= ListGradient[{{1, 2, 6}, {3, 4, 5}}]
Out[11]= {{{2, 2, -1}, {2, 2, -1}}, {{1, 5/2, 4}, {1, 1, 1}}}
```

Forward differences everywhere

```mathematica
In[12]:= ListGradient[{1, 2, 4, 7, 11}, Method -> "Forward"]
Out[12]= {1/2, 3/2, 5/2, 7/2, 9/2}
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

**Algorithm.** `builtin_list_gradient` is a from-first-principles port of
`numpy.gradient`: the numerical gradient of a sampled rank-*k* array by finite
differences, one same-shape array per differentiated axis (a rank-1 `f` gives
one vector, a rank-*k* `f` gives `{g_1, ..., g_k}`). The default reproduces
numpy exactly — a second-order central difference in the interior and a
first-order one-sided difference at the two extreme endpoints (numpy
`edge_order = 1`). One kernel underlies *every* Method/order/window/grid
combination: Fornberg's finite-difference-weight recurrence, which returns the
first-derivative weights `w[]` of `f'(x0) ~= sum_k w[k] f(z[k])` on arbitrarily
spaced nodes — so a coordinate vector (a non-uniform grid) needs no special
case. `lg_stencil` picks the node window per output index: symmetric in the
interior, shifting inward near an edge, and for the Centered method dropping the
two extreme endpoints to a one-sided `(m-1)`-point stencil; Forward/Backward
stay order `p` everywhere. Options are `Method -> "Centered"|"Forward"|
"Backward"`, `DifferenceOrder -> p` (default 2), `WindowLength -> m`
(default `p+1`), and `Axis -> All | a | {a1, ...}`.

**Data structures.** Two instantiations of the one recurrence.
`fd_weights_double` drives the machine-buffer fast path — a `float64`/`float32`
NDArray, packed or visible; its interior run is a fused tap-outer `axpy` over the
contiguous block (special-cased for the common 2- and 3-point windows) that the
compiler turns into SIMD FMA, with uniform-spacing weight vectors memoised by
`(count, offset)`. A complex machine buffer is carried through the same real
weights over its interleaved `(re, im)` doubles. `fd_weights_expr` builds `Expr`
weights the evaluator reduces — exact Rationals for an integer/rational grid,
symbolic weights for a symbolic grid — so integer and symbolic inputs answer
exactly. The buffer path declines (returns `NULL`) for anything it cannot
represent in a float buffer (int64/bool dtype, symbolic spacing), and
`ndstruct_delist_repack` re-enters the exact/symbolic List path with a
materialised List. On a packed multi-axis call the stacked rank-(*k*+1) output is
built in place to skip the transparency gate's copy; a visible NDArray keeps the
List-of-arrays shape.

**Complexity / limits.** `O(N * m)` for `N` elements and window `m`.
`ListGradient` is on `pack.c`'s `AWARE` list but deliberately **not**
`INT64_OK`: the central difference of exact integers is a Rational
(`(f[i+1]-f[i-1])/2`), which no int64 buffer holds, so an int64 argument
materialises and the exact List path answers. The rank-1 real-array,
default-centered form lowers in `Compile[]` (`compile_ndtables.c`) and
auto-compiles. A differentiated axis of length `< 2` is undefined and the call
declines.

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

- B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. **51** (1988) 699-706.
- Source: [`src/list/list_gradient.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_gradient.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_list_gradient.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list_gradient.c)

## Notes & additional examples

### Notes

`ListGradient[f]` is Mathilda's port of `numpy.gradient`: it estimates the
derivative of a sampled array by finite differences. By default it uses a
second-order central difference in the interior and a first-order one-sided
difference at the two extreme endpoints, so the result has the same shape as `f`.
A vector returns a vector; a rank-*k* array returns `{g_1, ..., g_k}`, one array
per axis.

Spacing is given by the optional second argument: a scalar `h` is uniform on
every axis, `{s_1, ..., s_k}` gives one spec per axis, and a coordinate list the
same length as the axis samples a **non-uniform** grid — handled with no special
case, because the underlying Fornberg weights are exact on arbitrarily spaced
nodes. Options `Method`, `DifferenceOrder`, `WindowLength` and `Axis` select the
one-sided variants, the accuracy order, the stencil size, and which axes to
differentiate.

Exactness follows the input: an integer or rational grid yields exact Rationals
(`(f[i+1]-f[i-1])/2` and friends), a symbolic grid yields symbolic weights, and a
machine-real (or complex) array runs straight on the packed float buffer and
stays packed. Because an integer gradient is Rational and no int64 buffer holds
one, an exact integer input is routed to the List path rather than the buffer.
The rank-1 real-array form also lowers inside `Compile[]`. Each differentiated
axis must have length at least 2.
