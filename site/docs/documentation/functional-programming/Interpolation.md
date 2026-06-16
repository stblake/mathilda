# Interpolation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Interpolation[data]
    constructs an InterpolatingFunction that interpolates data, given
    as {f1, f2, ...} (values at x = 1, 2, ...), {{x1, f1}, ...} (values
    at given abscissae), or {{{x1, y1, ...}, f1}, ...} (an m-D tensor
    grid).
Interpolation[data, x]
    builds the interpolating function and evaluates it at x (a number,
    or a coordinate list in m-D).
Interpolation[{{{x1,...}, f1, df1, ddf1, ...}, ...}]
    reproduces supplied derivatives at the nodes (df = gradient, ddf =
    Hessian, ...) by tensor-product Hermite interpolation.
Interpolation[data, InterpolationOrder -> n]
    uses piecewise-polynomial pieces of degree n (default 3; 0 gives a
    piecewise-constant and 1 a piecewise-linear interpolant).
Interpolation[data, Method -> m]
    selects "Spline" (natural/cyclic cubic spline) or "Hermite"
    (piecewise cubic Hermite with estimated slopes).
Interpolation[data, PeriodicInterpolation -> True]
    builds a periodic interpolant (period = the data span; the data must
    repeat its first sample at the last). A per-dimension {True, False}
    list selects periodicity per axis.
Vector- or array-valued samples (f_i a list) are interpolated
component-wise and return an array of the same shape.
Works at machine or arbitrary (MPFR) precision, matching the data.
```

## Examples

_No verified examples yet for this function._

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

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- C. de Boor, *A Practical Guide to Splines*, rev. ed. (Springer, 2001).
- Source: [`src/interp.c`](https://github.com/stblake/mathilda/blob/main/src/interp.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= f = Interpolation[{1, 4, 9, 16}]
Out[1]= InterpolatingFunction[{{1, 4}}, <>]

In[2]:= f[2]
Out[2]= 4

In[3]:= f[2.5]
Out[3]= 6.25
```

```mathematica
In[1]:= g = Interpolation[Table[{x, Sin[x]}, {x, 0., 6., 0.5}]]
Out[1]= InterpolatingFunction[{{0.0, 6.0}}, <>]

In[2]:= g[1.5]
Out[2]= 0.997495
```

```mathematica
In[1]:= p = Interpolation[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}, InterpolationOrder -> 2]
Out[1]= InterpolatingFunction[{{0, 3}}, <>]

In[2]:= p[1.5]
Out[2]= 2.25
```

### Notes

`Interpolation[data]` builds an `InterpolatingFunction` (default piecewise
cubic) over the given samples; here `{1, 4, 9, 16}` are the values at
`x = 1, 2, 3, 4`, so evaluating recovers `x^2`. The returned object prints with
only its domain shown and is applied like an ordinary function, `f[x]`.
