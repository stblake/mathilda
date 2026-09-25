# ListInterpolation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListInterpolation[array]`**

constructs an InterpolatingFunction interpolating the values in array, taken to lie on a regular grid at integer positions 1, 2, ... in each direction. array may be nested to any depth; the nesting depth is the number of dimensions.

**`ListInterpolation[array, {{xmin, xmax}, {ymin, ymax}, ...}]`**

places the grid lines equally spaced across the given interval in each direction (one {min, max} pair per dimension).

**`ListInterpolation[array, {{x1, x2, ...}, ...}]`**

uses explicit lists of grid-line positions in each direction.

**`ListInterpolation[array, ..., InterpolationOrder -> n]`**

sets the piecewise-polynomial degree (default 3; 0 constant, 1 linear); Method -\> "Spline" | "Hermite" selects the method and PeriodicInterpolation -\> True builds a periodic interpolant.

<details>
<summary>Notes</summary>

A fast path for List, packed, and NDArray value tensors; works at machine or arbitrary (MPFR) precision, matching the data.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= f = ListInterpolation[{1, 2, 3, 5, 8, 5}]
Out[1]= InterpolatingFunction[{{1, 6}}, <>]

In[2]:= f[2.5]
Out[2]= 2.4375

In[3]:= ListInterpolation[{1, 2, 3, 5, 8, 5}, {{0, 1}}][0.5]
Out[3]= 3.875

In[4]:= ListInterpolation[{{1, 2}, {3, 4}}][1.5, 1.5]
Out[4]= 2.5
```

## Algorithm

interp.c

InterpolatingFunction --- piecewise-polynomial interpolation of tabulated data on a regular (tensor-product) grid, plus the Interpolation[] builder. Modelled on Mathematica's InterpolatingFunction object.

```text
  InterpolatingFunction[domain, table]
  InterpolatingFunction[domain, table, ders]
  InterpolatingFunction[domain, table, ders, orders]
  InterpolatingFunction[domain, table, ders, orders, method]

    domain = {{x1min, x1max}, ...}   -- one interval per dimension; the
             number of intervals m is the dimensionality.
    table  = {{coord, val}, ...}                     -- value-only data, or
             {{coord, val, grad, hess, ...}, ...}     -- derivative-supplied.
             coord is a scalar (1-D value-only) or an {x1,...,xm} list.
             grad = D[f,{vars,1}] (length-m vector), hess = D[f,{vars,2}]
             (m x m matrix), etc.
    ders   = {d1, ..., dm}   -- (optional) derivative-of-interpolant orders.
    orders = {o1, ..., om}   -- (optional) interpolation order per dimension.
    method = "Spline" | "Hermite"   -- (optional) interpolation method.
```

Methods (all evaluate the ders-th mixed derivative so D[ifun[..],..] composes):

```text
  default  : sliding-window Newton divided-difference (order min(3,n-1) or the
             requested InterpolationOrder), per dimension, tensor product.
  "Spline" : natural cubic spline (C2; second derivative 0 at the ends),
             tensor product over the full grid.
  "Hermite": tensor-product piecewise cubic Hermite with node slopes estimated
             by 3-point finite differences.
  supplied : derivative-annotated data is interpolated by tensor-product
             Hermite of per-dimension order k = max(K,1) where K is the highest
             supplied derivative order.  Mixed partials that are not supplied
             are filled by central finite differences across the grid.
```

Precision: machine (double) by default; if the data/argument carry MPFR arbitrary precision the MPFR kernels (interp_mpfr.c) are used instead and an EXPR_MPFR is returned.

Builtin ownership: interp_apply / the Interpolation builtin return a fresh Expr* (or NULL to stay unevaluated); inputs are borrowed.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Interpolation evaluate, 20000 points | 18.9 s | 18.1 s | 0.073 s |
| Interpolation over 10^5 array | 1.26 s | 4.57 s | 4.44 s |
| Interpolation build, 2000 knots | 0.013 s | 0.761 s | 0.092 s |
| Interpolation order 1 (linear) build | -- | 0.764 s | 0.014 s |
| ListInterpolation 2-D 60x60 | -- | 0.823 s | 1.02 s |

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [Interpolation](../../functional-programming/Interpolation/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/), [Rational](../../arithmetic/Rational/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_interp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interp.c)
