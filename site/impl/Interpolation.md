---
references:
  - "C. de Boor, *A Practical Guide to Splines*, rev. ed. (Springer, 2001)."
source: src/interp.c
---
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
