# InterpolatingFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InterpolatingFunction[domain, table]`**

represents an approximate function whose values are found by interpolation. domain is {{x1min, x1max}, ...} with one interval per dimension; table is a list of {coord, value} data points on a full tensor grid (coord is a scalar for 1-D, an {x1, ..., xm} list for m-D).

**`InterpolatingFunction[...][x1, ..., xm]`**

gives the interpolated value using tensor-product piecewise- polynomial (default order 3) interpolation. Arguments outside the domain are extrapolated with a warning.

**`Derivative[d1, ..., dm][InterpolatingFunction[...]]`**

gives an InterpolatingFunction for the mixed partial derivative.

<details>
<summary>Notes</summary>

In standard output only the domain is shown; the rest is \<\>.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ifun = InterpolatingFunction[{{0, 5}}, {{0,0},{1,1},{2,3},{3,4},{4,3},{5,0}}]
Out[1]= InterpolatingFunction[{{0, 5}}, <>]

In[2]:= ifun[2.5]
Out[2]= 3.6875

In[3]:= ifun[2]
Out[3]= 3

In[4]:= ifun'[2.5]
Out[4]= 1.04167

In[5]:= Derivative[0, 1][f][1.5, 2.5]
Out[5]= Derivative[0, 1][f][1.5, 2.5]
```

### Applications (6)

```mathematica
In[6]:= f = Interpolation[{1, 4, 9, 16}]
Out[6]= InterpolatingFunction[{{1, 4}}, <>]

In[7]:= f[2.5]
Out[7]= 6.25

In[8]:= g = Interpolation[Table[{x, Sin[x]}, {x, 0., 6., 0.5}]]
Out[8]= InterpolatingFunction[{{0.0, 6.0}}, <>]

In[9]:= g[1.5]
Out[9]= 0.997495

In[10]:= Sin[1.5]
Out[10]= 0.997495

In[11]:= d = Interpolation[{1, 4, 9, 16, 25}]; dd = d'; dd[2.5]
Out[11]= 5.0
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Interpolation evaluate, 20000 points | 18.9 s | 18.1 s | 0.073 s |
| NDS Van der Pol mu=1000 | 4.02 s | 0.234 s | 0.362 s |
| Interpolation over 10^5 array | 1.26 s | 4.57 s | 4.44 s |
| NDS harmonic t=100 | 0.884 s | 0.475 s | 37.2 s |
| NDS Van der Pol mu=100 | 0.647 s | 0.271 s | 0.452 s |
| NDS Lorenz t=5 | 0.337 s | 0.45 s | 14.4 s |

## Implementation notes

**The object.** `InterpolatingFunction[domain, table, ders, orders, method,
periodic]` is the callable normal form produced by `Interpolation`. It carries
the attribute `ATTR_HOLDALL` (set in `interp_init`): without it, every application
`ifun[x]` would force the evaluator to re-traverse and re-evaluate the entire
N-point table on each call (an O(N) cost); holding the constant `domain`/`table`
arguments keeps the object a true normal form so evaluating it is O(1), while the
application's own argument is evaluated by the surrounding head. Application is
reduced in `eval.c` (section 7d), which calls `interp_apply`; `D[ifun, ...]`
folds into a derivative object via `interp_make_derivative` (it just increments
the per-dimension `ders` orders and rewrites the object).

**Evaluation.** `interp_apply` parses `ders`, `orders`, `method`
(`"Spline"`/`"Hermite"`/default), and per-dimension `periodic` flags, coerces the
call coordinates to `double` (`node_to_double`), and dispatches. Exact
node-coincident queries short-circuit to the stored exact value. If the
data/argument carry MPFR precision the MPFR kernels (`interp_eval_mpfr`,
`interp_mpfr.c`) run; otherwise `interp_eval_double`. Per dimension a coordinate
is bracketed by binary search `bracket_interval` (clamping to the end intervals,
extrapolating with a `dmval` warning outside the data range, or reduced modulo the
period for periodic dimensions).

**Kernels.** The default and spline kernels evaluate a separable tensor product
recursively (`eval_dim`): the innermost dimension reads the flat value tensor
`f->V` and outer dimensions recurse, so the m-D interpolant is a nest of 1-D
evaluations. Each 1-D evaluation returns the requested `ders`-th derivative
directly: `newton_deriv_eval` builds Newton divided differences on a sliding
window and differentiates the nested (Horner-style) product form; `spline_eval` /
`spline_eval_periodic` solve the (natural or cyclic) cubic-spline system and
evaluate the chosen derivative of the bracketed cubic. The Hermite/supplied path
(`hermite_tensor_eval`, `build_basis`, `build_T`) interpolates derivative-annotated
data with tensor-product piecewise cubic Hermite, filling unsupplied mixed
partials by finite differences.

**Data structures.** The parsed grid is the `IFun` struct: per-dimension sorted
abscissae `grid[k]` (built by `grid_insert`, indexed by `grid_index`), strides for
row-major flattening, domain bounds, periodicity/period arrays, a node→data-entry
map `entryAt`, and the flat value tensor `V`. Parsing `domain`+`table` into an
`IFun` is memoised by a one-slot module cache `g_grid_cache` keyed by
`expr_eq(domain, table)` + dimensionality + periodicity (`grid_cache_get`), so
repeated applications of the same object reuse the grid (the returned `IFun` is
borrowed, owned by the cache).

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Function](../../functional-programming/Function/), [FullForm](../../expression-information/FullForm/)

- C. de Boor, *A Practical Guide to Splines*, rev. ed. (Springer, 2001).
- E. H. Neville, "Iterative Interpolation", J. Indian Math. Soc. 1934.
- Source: [`src/interp.c`](https://github.com/stblake/mathilda/blob/main/src/interp.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_interp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interp.c)
- Tests: [`tests/test_ndsolve_pde.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_pde.c)

## Notes & additional examples

### Notes

`InterpolatingFunction[domain, table]` is the approximate-function object
returned by `Interpolation`; in standard output only the domain
(`{{xmin, xmax}, ...}`) is shown and the data table is elided as `<>`. Apply it
to a coordinate, `f[x]`, to get the tensor-product piecewise-polynomial value;
arguments outside the domain are extrapolated with a warning.
