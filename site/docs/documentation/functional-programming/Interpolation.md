# Interpolation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Interpolation[data]`**

constructs an InterpolatingFunction that interpolates data, given as {f1, f2, ...} (values at x = 1, 2, ...), {{x1, f1}, ...} (values at given abscissae), or {{{x1, y1, ...}, f1}, ...} (an m-D tensor grid).

**`Interpolation[data, x]`**

builds the interpolating function and evaluates it at x (a number, or a coordinate list in m-D).

**`Interpolation[{{{x1,...}, f1, df1, ddf1, ...}, ...}]`**

reproduces supplied derivatives at the nodes (df = gradient, ddf = Hessian, ...) by tensor-product Hermite interpolation.

**`Interpolation[data, InterpolationOrder -> n]`**

uses piecewise-polynomial pieces of degree n (default 3; 0 gives a piecewise-constant and 1 a piecewise-linear interpolant).

**`Interpolation[data, Method -> m]`**

selects "Spline" (natural/cyclic cubic spline) or "Hermite" (piecewise cubic Hermite with estimated slopes).

**`Interpolation[data, PeriodicInterpolation -> True]`**

builds a periodic interpolant (period = the data span; the data must repeat its first sample at the last). A per-dimension {True, False} list selects periodicity per axis.

<details>
<summary>Notes</summary>

Vector- or array-valued samples (f\_i a list) are interpolated component-wise and return an array of the same shape. Works at machine or arbitrary (MPFR) precision, matching the data.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= f = Interpolation[{1, 2, 3, 5, 8, 5}]
Out[1]= InterpolatingFunction[{{1, 6}}, <>]

In[2]:= f[2.5]
Out[2]= 2.4375
```

F = x^3 with f and f' supplied: cubic reproduced exactly

```mathematica
In[3]:= c = Interpolation[{{{0}, 0, 0}, {{1}, 1, 3}, {{2}, 8, 12}, {{3}, 27, 27}}]; {c[1.5], c'[1.5]}
Out[3]= {3.375, 6.75}
```

Vector-valued: each component interpolated independently

```mathematica
In[4]:= fv = Interpolation[{{{0.}, {1., 2.}}, {{1.2}, {3., 4.}}, {{2.1}, {5., 4.}}, {{3.}, {0., 4.}}}]; fv[1.5]
Out[4]= {4.03175, 4.07143}
```

### Options (4)

```mathematica
In[5]:= Interpolation[{1, 5, 7, 2, 3, 1}, InterpolationOrder -> 1][2.5]
Out[5]= 6.0
```

X^2

```mathematica
In[6]:= Interpolation[{1, 4, 9, 16, 25}, Method -> "Hermite"][2.5]
Out[6]= 6.25
```

High-precision data -> high-precision result

```mathematica
In[7]:= Interpolation[N[{1, 2, 3, 5, 8, 5}, 30], Method -> "Spline"][N[5/2, 30]]
Out[7]= 2.473086124401913875598086124405
```

Periodic: f[x] wraps with period = data span (5)

```mathematica
In[8]:= fp = Interpolation[Table[{x, N[Sin[2 Pi x/5]]}, {x, 0, 5}], PeriodicInterpolation -> True]; {fp[0.5], fp[5.5], fp[-4.5]}
Out[8]= {0.557674, 0.557674, 0.557674}
```

### Applications (7)

```mathematica
In[9]:= f = Interpolation[{1, 4, 9, 16}]
Out[9]= InterpolatingFunction[{{1, 4}}, <>]

In[10]:= f[2]
Out[10]= 4

In[11]:= f[2.5]
Out[11]= 6.25

In[12]:= g = Interpolation[Table[{x, Sin[x]}, {x, 0., 6., 0.5}]]
Out[12]= InterpolatingFunction[{{0.0, 6.0}}, <>]

In[13]:= g[1.5]
Out[13]= 0.997495

In[14]:= p = Interpolation[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}, InterpolationOrder -> 2]
Out[14]= InterpolatingFunction[{{0, 3}}, <>]

In[15]:= p[1.5]
Out[15]= 2.25
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

**Algorithm.** `builtin_interpolation` is the *builder*: it parses tabulated data
into an `InterpolatingFunction[domain, table, ...]` normal-form object (it does
not itself evaluate the interpolant — that is the callable object's job). It
recognises three data forms: a bare list of values (form 1, where abscissae are
synthesised as `1,2,3,...`), `{{x, y}, ...}` (1-D value pairs), and
`{{{x1,...,xm}, y}, ...}` for `m`-dimensional value-only data, or `{{coord, val,
grad, hess, ...}, ...}` for derivative-supplied data (the number of trailing
tensors is `Ksupplied = L - 2`). Options `InterpolationOrder -> o`, `Method ->
"Spline"|"Hermite"`, and `PeriodicInterpolation -> True|False|{...}` are read from
`Rule`/`RuleDelayed` arguments; a lone non-option argument is taken as an
immediate evaluation point.

It then constructs the table of `{coord, val, ...}` entries (synthesising integer
coordinates for form 1) and a per-dimension `domain = {{min,max}, ...}` computed
from the coordinate extrema, preferring the *exact* boundary `Expr`s
(`dminE`/`dmaxE`) when available. The object is emitted with the minimal arity
needed: just `{domain, table}` by default, or with explicit `ders` (all zero),
`orders`, a method slot, and a periodicity list when any non-default option is
present. If an evaluation point was supplied, it immediately calls `interp_apply`
on the freshly built object and returns the value instead of the object.

**Data structures.** Coordinates are pulled to `double` via `node_to_double` for
extent/grid bookkeeping (with a 64-dimension cap), but the stored table keeps the
original exact `Expr` nodes (`expr_copy`). The actual interpolation grid is not
built here — it is constructed lazily and cached when the object is first applied
(see `InterpolatingFunction`).

**Limits.** Requires `>= 2` points; data must fill a full tensor-product grid
(enforced later by `build_grid`). The default method is sliding-window Newton
divided-difference (order `min(3, n-1)` unless `InterpolationOrder` overrides);
`"Spline"` selects a natural/periodic cubic spline and `"Hermite"` a
tensor-product piecewise cubic Hermite.

**Attributes:** `Protected`.

## References

**See also:** [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/)

- C. de Boor, *A Practical Guide to Splines*, rev. ed. (Springer, 2001).
- Source: [`src/interp.c`](https://github.com/stblake/mathilda/blob/main/src/interp.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_interp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interp.c)

## Notes & additional examples

### Notes

`Interpolation[data]` builds an `InterpolatingFunction` (default piecewise
cubic) over the given samples; here `{1, 4, 9, 16}` are the values at
`x = 1, 2, 3, 4`, so evaluating recovers `x^2`. The returned object prints with
only its domain shown and is applied like an ordinary function, `f[x]`.
