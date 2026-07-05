# Functional Programming

## Distribute
Implements the distributive law for operators `f` and `g`.
- `Distribute[expr]`: distributes the head of `expr` over `Plus` appearing in any of its arguments.
- `Distribute[expr, g]`: distributes over `g`.
- `Distribute[expr, g, f]`: performs the distribution only if the head of `expr` is `f`.
- `Distribute[expr, g, f, gp, fp]`: gives `gp` and `fp` in place of `g` and `f` respectively in the result of the distribution.

**Features**:
- `Protected`.
- `Distribute` explicitly constructs the complete result of a distribution; `Expand`, on the other hand, builds up results iteratively, simplifying at each stage.
- For pure products, `Distribute` gives the same results as `Expand`.

```mathematica
In[1]:= Distribute[(a + b) . (x + y + z)]
Out[1]= a . x + a . y + a . z + b . x + b . y + b . z

In[2]:= Distribute[f[a + b, c + d + e]]
Out[2]= f[a, c] + f[a, d] + f[a, e] + f[b, c] + f[b, d] + f[b, e]

In[3]:= Distribute[(a + b + c) (u + v), Plus, Times]
Out[3]= a u + a v + b u + b v + c u + c v

In[4]:= Distribute[{{a, b}, {x, y, z}, {s, t}}, List]
Out[4]= {{a, x, s}, {a, x, t}, {a, y, s}, {a, y, t}, {a, z, s}, {a, z, t}, {b, x, s}, {b, x, t}, {b, y, s}, {b, y, t}, {b, z, s}, {b, z, t}}

In[5]:= Distribute[{{}, {a}}, {{}, {b}}, {{}, {c}}, List, List, List, Join]
Out[5]= {{}, {c}, {b}, {b, c}, {a}, {a, c}, {a, b}, {a, b, c}}
```

## Pure Functions (Function, &)
`Function` is a pure (or "anonymous") function, analogous to $\lambda$ in LISP.
- `Function[body]` or `body &`: A pure function where formal parameters are `#` (or `#1`), `#2`, etc.
- `Function[x, body]` or `Function[{x1, x2, ...}, body]`: A pure function with named local formal parameters.
- `Function[params, body, attrs]`: A pure function treated as having the evaluation attributes `attrs`.
- `Function[Null, body, attrs]`: Slot form (`#`, `##`, etc.) with evaluation attributes.
- `Slot` (`#`): Represents the first argument. `#n` represents the $n$-th argument.
- `SlotSequence` (`##`): Represents all arguments sequence. `##n` represents all arguments starting from the $n$-th.

**Features**:
- **Lexical parameter binding**: named parameters are substituted into the body before evaluation (not bound via the global symbol table). This means held references to a parameter (for example `Unevaluated[x]` inside the body) see the substituted expression rather than the raw symbol.
- **No Hold by default**: the default `Function` has no hold attributes, so its arguments are evaluated before substitution, matching Mathematica. `(Hold[#]&)[1+2]` gives `Hold[3]`, not `Hold[1 + 2]`.
- **3-argument form attributes**: `Function[params, body, attrs]` can assign any of the standard evaluator attributes. Recognised attributes include `HoldFirst`, `HoldRest`, `HoldAll`, `HoldAllComplete`, `Listable`, `Flat`, `Orderless`, `OneIdentity`, `NumericFunction`, `SequenceHold`, `NHoldRest`. `attrs` may be a single attribute symbol or a list.
- `Slot` parameters are positionally mapped to arguments provided. Remaining arguments are ignored.
- Pure functions can be assigned to variables (e.g., `f = #^2 &`) and applied over lists (e.g., `f /@ {1, 2, 3}`).
- Properly scopes nested `Function` expressions: named inner variables shadow outer ones, and unnamed inner functions establish a new scope for `#` slots.

```mathematica
In[1]:= Function[u, 3 + u][x]
Out[1]= 3 + x

In[2]:= (#1^2 + #2^4)&[x, y]
Out[2]= x^2 + y^4

In[3]:= f[X, ##, Y, ##]&[a, b, c]
Out[3]= f[X, a, b, c, Y, a, b, c]

In[4]:= Function[{x}, Length[Unevaluated[x]], {HoldAll}][1+1+1]
Out[4]= 3

In[5]:= Function[{x}, Length[Unevaluated[x]]][1+1+1]
Out[5]= 0

In[6]:= (Hold[#]&)[1+2]
Out[6]= Hold[3]

In[7]:= Function[{a, b}, {Length[Unevaluated[a]], Length[Unevaluated[b]]}, HoldFirst][1+2+3, 4+5+6]
Out[7]= {3, 0}
```

## InterpolatingFunction
`InterpolatingFunction[domain, table]` represents an approximate function
whose values are found by interpolation, and — like `Function` — is evaluated
by applying it to arguments. Any number of real arguments (dimensions) is
supported.

- `domain`: one interval per dimension, `{{x1min, x1max}, ..., {xmmin, xmmax}}`.
  The number of intervals `m` is the dimensionality.
- `table`: the data points, `{{coord1, val1}, {coord2, val2}, ...}`. For `m = 1`
  each `coord` is a scalar `x`; for `m ≥ 2` each `coord` is a list
  `{x1, ..., xm}`. The points must fill the Cartesian product of the
  per-dimension grids (a full tensor grid). Points may be given in any order;
  they are sorted by abscissa internally.
- `InterpolatingFunction[...][x1, ..., xm]` returns the interpolated value.
- `InterpolatingFunction[domain, table, {d1, ..., dm}]` is the optional
  derivative-annotated form (see **Derivatives**).

**Method**: tensor-product piecewise-polynomial interpolation of default
order 3 (cubic) in every dimension. For a query point the bracketing interval
is found in each dimension and a sliding window of `order + 1` consecutive grid
points, centred on that interval, is selected; the value is the tensor product
of the per-dimension Newton divided-difference polynomials. The order drops to
`min(3, n_k - 1)` when a dimension has fewer than four grid points (two →
linear, three → quadratic). The windowing reproduces Mathematica's default
piecewise interpolant.

**Derivatives**: `Derivative[d1, ..., dm][InterpolatingFunction[...]]` (and so
`ifun'`, `D[ifun[x], x]`, `D[f[x,y], {x,2}]`, etc.) returns another
`InterpolatingFunction` carrying the accumulated derivative orders — exactly as
in Mathematica. When applied it evaluates the mixed partial
`∂^{d1}_{x1} ... ∂^{dm}_{xm}` of the windowed interpolant. Orders are additive
across repeated differentiation, and an order beyond the local polynomial
degree yields `0`.

**Semantics**:
- An *exact* argument tuple (Integer/Rational) coinciding with a grid node
  returns that node's stored value exactly (value queries only); all other
  queries return a machine real.
- An argument outside `domain` is extrapolated from the nearest window, and an
  `InterpolatingFunction::dmval` warning is issued.
- A symbolic argument, a wrong number of arguments (≠ `m`), or a malformed
  object (fewer than two grid points in a dimension, an incomplete tensor grid,
  non-numeric data) leaves the application unevaluated.
- In standard output only the `domain` is shown; the rest is abbreviated as
  `<>`. `FullForm` reveals the full structure.

Attribute: `Protected`.

```mathematica
In[1]:= ifun = InterpolatingFunction[{{0, 5}}, {{0,0},{1,1},{2,3},{3,4},{4,3},{5,0}}]
Out[1]= InterpolatingFunction[{{0, 5}}, <>]

In[2]:= ifun[2.5]
Out[2]= 3.6875

In[3]:= ifun[2]
Out[3]= 3

In[4]:= ifun'[2.5]
Out[4]= 1.04167

In[5]:= f = InterpolatingFunction[{{0, 4}, {0, 4}}, {{{0,0},0}, ..., {{4,4},80}}];
        f[1.5, 2.5]              (* tensor-product 2-D interpolation *)
Out[5]= 17.875

In[6]:= Derivative[0, 1][f][1.5, 2.5]
Out[6]= 18.75
```

## Interpolation
`Interpolation[data]` constructs an [`InterpolatingFunction`](#interpolatingfunction)
object that agrees with `data` at every supplied point. The data may be given as:

- `{f1, f2, ...}` — function values at `x = 1, 2, ...`.
- `{{x1, f1}, {x2, f2}, ...}` — values at the given abscissae.
- `{{{x1, y1, ...}, f1}, ...}` — values on a full m-dimensional tensor grid.
- `{{{x1, ...}, f1, df1, ddf1, ...}, ...}` — values **and supplied derivatives**
  at each node; `df` is the gradient (`D[f,{vars,1}]`, a scalar in 1-D, a length-m
  vector in m-D), `ddf` the Hessian (`D[f,{vars,2}]`), and so on.

`Interpolation[data, x]` builds the function and immediately evaluates it at the
point `x` (a number, or an `{x, y, ...}` coordinate list in m-D).

**Methods.** Value-only data uses sliding-window Newton interpolation by default.

- `Method -> "Spline"` — natural cubic spline (C², second derivative zero at the
  ends; reduces to linear for two points).
- `Method -> "Hermite"` — piecewise cubic Hermite with node slopes estimated by
  3-point finite differences (reproduces quadratics exactly).
- Derivative-supplied data is interpolated by **tensor-product Hermite** of
  per-dimension order `k = max(K, 1)`, where `K` is the highest supplied
  derivative order. Mixed partials not present in the data (e.g. `f_xy` from
  gradient-only data) are filled by central finite differences across the grid.
  This reproduces a cubic from `{f, f'}` and a quintic from `{f, f', f''}`
  exactly, and likewise for separable polynomials in m-D.

**Options.**

- `InterpolationOrder -> n` — degree of the piecewise-polynomial pieces
  (default `3`; `0` piecewise-constant, `1` piecewise-linear). Carried on the
  object and preserved through differentiation.
- `PeriodicInterpolation -> True` — build a periodic interpolant. The period in
  each dimension is the data span, and the data must repeat its first sample at
  the last (`f(x_max) = f(x_min)`). Out-of-range arguments wrap (no `dmval`
  warning); the seam is handled by wrap-around windows for the default/Hermite
  methods and a **cyclic** cubic spline for `Method -> "Spline"`. A per-dimension
  list `{True, False, …}` selects periodicity per axis.

**Vector/array-valued data.** A sample `f_i` may itself be a list or array (with
explicit coordinates, e.g. `{{x}, {v1, v2}}`). Each component is interpolated
independently and the result is an array of the same shape. This composes with
every method, with supplied derivatives (whose tensors carry the value-array axes
first, then the derivative axes), in n-D, and at MPFR precision.

The domain is inferred from the data (the min/max abscissa in each dimension),
preserving the data's number type (integer abscissae give an integer domain).

**Precision.** Interpolation is carried out at machine precision for machine data
and at the data's MPFR precision when the data or the evaluation point carry
arbitrary precision; the result is an `EXPR_MPFR` real in that case. All methods
and data forms support both precisions. The interpolants are mathematically
standard (passing through every data point and honouring every supplied
derivative); interior values need not match Wolfram's B-spline output bit-for-bit.

Attribute: `Protected`.

```mathematica
In[1]:= f = Interpolation[{1, 2, 3, 5, 8, 5}]
Out[1]= InterpolatingFunction[{{1, 6}}, <>]

In[2]:= f[2.5]
Out[2]= 2.4375

In[3]:= Interpolation[{1, 5, 7, 2, 3, 1}, InterpolationOrder -> 1][2.5]
Out[3]= 6.

In[4]:= Interpolation[{1, 4, 9, 16, 25}, Method -> "Hermite"][2.5]  (* x^2 *)
Out[4]= 6.25

In[5]:= (* f = x^3 with f and f' supplied: cubic reproduced exactly *)
        c = Interpolation[{{{0}, 0, 0}, {{1}, 1, 3}, {{2}, 8, 12}, {{3}, 27, 27}}];
        {c[1.5], c'[1.5]}
Out[5]= {3.375, 6.75}

In[6]:= (* high-precision data -> high-precision result *)
        Interpolation[N[{1, 2, 3, 5, 8, 5}, 30], Method -> "Spline"][N[5/2, 30]]
Out[6]= 2.473086124401913875598086124405

In[7]:= (* vector-valued: each component interpolated independently *)
        fv = Interpolation[{{{0.}, {1., 2.}}, {{1.2}, {3., 4.}}, {{2.1}, {5., 4.}}, {{3.}, {0., 4.}}}];
        fv[1.5]
Out[7]= {4.03175, 4.07143}

In[8]:= (* periodic: f[x] wraps with period = data span (5) *)
        fp = Interpolation[Table[{x, N[Sin[2 Pi x/5]]}, {x, 0, 5}], PeriodicInterpolation -> True];
        {fp[0.5], fp[5.5], fp[-4.5]}
Out[8]= {0.557674, 0.557674, 0.557674}
```

## Map (/@)
- `f /@ expr` or `Map[f, expr]`

## Apply (@@, @@@)
- `f @@ expr`: Level 0.
- `f @@@ expr`: Level 1.

## MapAll (//@)
- `f //@ expr`: Recursive map.

## MapAt
Applies a function to selected parts of an expression, identified by position specifications of the same form as those returned by `Position`.

- `MapAt[f, expr, n]`: applies `f` to the element at position `n` in `expr`. If `n` is negative, the position is counted from the end. `n = 0` targets the head of `expr`.
- `MapAt[f, expr, {i, j, ...}]`: applies `f` to the part of `expr` at position `{i, j, ...}` (equivalently `expr[[i, j, ...]]`).
- `MapAt[f, expr, {{i1, j1, ...}, {i2, j2, ...}, ...}]`: applies `f` to each of the listed parts of `expr`. If the same position appears more than once, `f` is applied repeatedly to that part.

**Features**:
- `Protected`.
- Path components may be integers (positive or negative), `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- `MapAt[f, expr, {}]` applies `f` to the entire expression itself.

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, 2]
Out[1]= {a, f[b], c, d}

In[2]:= MapAt[f, {a, b, c, d}, {{1}, {4}}]
Out[2]= {f[a], b, c, f[d]}

In[3]:= MapAt[f, {{a, b, c}, {d, e}}, {2, 1}]
Out[3]= {{a, b, c}, {f[d], e}}

In[4]:= MapAt[f, {{a, b, c}, {d, e}}, {All, 2}]
Out[4]= {{a, f[b], c}, {d, f[e]}}

In[5]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, -3]
Out[5]= {{a, b, c}, h[{d, e}], f, g}

In[6]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {2, 1}]
Out[6]= {{a, b, c}, {h[d], e}, f, g}

In[7]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{2}, {1}}]
Out[7]= {h[{a, b, c}], h[{d, e}], f, g}

In[8]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{1, 1}, {2, 2}, {3}}]
Out[8]= {{h[a], b, c}, {d, h[e]}, h[f], g}

In[9]:= MapAt[f, {1, 2, 3, 4, 5, 6}, 3 ;; 4]
Out[9]= {1, 2, f[3], f[4], 5, 6}

In[10]:= MapAt[f, a + b + c + d, 2]
Out[10]= a + c + d + f[b]

In[11]:= MapAt[f, x^2 + y^2, {{1, 1}, {2, 1}}]
Out[11]= f[x]^2 + f[y]^2

In[12]:= MapAt[f, {a, b, c}, 0]
Out[12]= f[List][a, b, c]
```

## Nest
Applies a function repeatedly.
- `Nest[f, expr, n]`: gives an expression with `f` applied `n` times to `expr`.

**Features**:
- `Protected`.
- `n` must be a non-negative integer; `Nest[f, expr, 0]` returns `expr` unchanged.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.

**Examples**:
```
In[1]:= Nest[f, x, 3]
Out[1]= f[f[f[x]]]

In[2]:= Nest[(1 + #)^2 &, 1, 3]
Out[2]= 676

In[3]:= Nest[(1 + #)^2 &, x, 5]
Out[3]= (1 + (1 + (1 + (1 + (1 + x)^2)^2)^2)^2)^2

In[4]:= Nest[Sqrt, 100.0, 4]
Out[4]= 1.33352

In[5]:= Nest[1/(1 + #) &, x, 5]
Out[5]= 1/(1 + 1/(1 + 1/(1 + 1/(1 + 1/(1 + x)))))

In[6]:= Nest[x^# &, x, 6]
Out[6]= x^x^x^x^x^x^x

In[7]:= Nest[#(1 + 0.05) &, 1000, 10]
Out[7]= 1628.89

In[8]:= Nest[(# + 2/#)/2 &, 1.0, 5]
Out[8]= 1.41421

In[9]:= Nest[{{1, 1}, {1, 0}} . # &, {0, 1}, 10]
Out[9]= {55, 34}
```

## NestList
Applies a function repeatedly, collecting every intermediate result.
- `NestList[f, expr, n]`: gives a list of the results of applying `f` to `expr` 0 through `n` times.

**Features**:
- `Protected`.
- Returns a list of length `n + 1` whose first element is `expr` and whose `(k+1)`-th element is `f` applied `k` times to `expr`.
- `n` must be a non-negative integer; `NestList[f, expr, 0]` returns `{expr}`.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.
- `Last[NestList[f, expr, n]]` is equivalent to `Nest[f, expr, n]`.

**Examples**:
```
In[1]:= NestList[f, x, 4]
Out[1]= {x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]}

In[2]:= NestList[Cos, 1.0, 10]
Out[2]= {1., 0.540302, 0.857553, 0.65429, 0.79348, 0.701369, 0.76396, 0.722102, 0.750418, 0.731404, 0.744237}

In[3]:= NestList[(1 + #)^2 &, x, 3]
Out[3]= {x, (1 + x)^2, (1 + (1 + x)^2)^2, (1 + (1 + (1 + x)^2)^2)^2}

In[4]:= NestList[Sqrt, 100.0, 4]
Out[4]= {100., 10., 3.16228, 1.77828, 1.33352}

In[5]:= NestList[2 # &, 1, 10]
Out[5]= {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024}

In[6]:= NestList[# + 1 &, 0, 10]
Out[6]= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

In[7]:= NestList[#^2 &, 2, 6]
Out[7]= {2, 4, 16, 256, 65536, 4294967296, 18446744073709551616}

In[8]:= NestList[#(1 + 0.05) &, 1000, 10]
Out[8]= {1000, 1050., 1102.5, 1157.63, 1215.51, 1276.28, 1340.1, 1407.1, 1477.46, 1551.33, 1628.89}

In[9]:= NestList[(# + 2/#)/2 &, 1.0, 5]
Out[9]= {1., 1.5, 1.41667, 1.41422, 1.41421, 1.41421}

In[10]:= NestList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 100, 20]
Out[10]= {100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2, 1}

In[11]:= NestList[Mod[59 #, 101] &, 1, 15]
Out[11]= {1, 59, 47, 46, 88, 41, 96, 8, 68, 73, 65, 98, 25, 61, 64, 39}
```

## NestWhile
Iteratively applies a function while a predicate continues to yield `True`.
- `NestWhile[f, expr, test]`: Starts with `expr` and keeps applying `f` until `test` no longer yields `True`.
- `NestWhile[f, expr, test, m]`: Supplies the most recent `m` results (not wrapped in a list) as arguments to `test`.
- `NestWhile[f, expr, test, All]`: Supplies every result so far as arguments to `test`.
- `NestWhile[f, expr, test, {mmin, mmax}]`: Defers the first test until `mmin` results exist, then supplies up to `mmax` most-recent results.
- `NestWhile[f, expr, test, m, max]`: Caps the number of `f` applications at `max` (may be `Infinity`).
- `NestWhile[f, expr, test, m, max, n]`: After the loop terminates, applies `f` an additional `n` times.
- `NestWhile[f, expr, test, m, max, -n]`: Returns the result produced `n` applications before the loop ended (i.e. `Part[NestWhileList[...], -n-1]`).

**Features**:
- `Protected`.
- If `test[expr]` does not yield `True` initially, the unchanged `expr` is returned.
- Results passed to `test` are in generation order with the most recent last, so e.g. `# > 1 &` inspects the oldest when more than one result is supplied.
- `NestWhile[f, expr, UnsameQ, 2]` is equivalent to `FixedPoint[f, expr]`.
- `NestWhile[f, expr, UnsameQ, All]` continues until any prior value reappears.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhile` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Examples**:
```
In[1]:= NestWhile[#/2 &, 123456, EvenQ]
Out[1]= 1929

In[2]:= NestWhile[Log, 100., # > 0 &]
Out[2]= -0.859384

In[3]:= NestWhile[Floor[#/2] &, 10, UnsameQ, 2]
Out[3]= 0

In[4]:= NestWhile[#/2 &, 123456, EvenQ, 1, 4]
Out[4]= 7716

In[5]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[5]= 1

In[6]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[6]= 0

In[7]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[7]= 2

In[8]:= NestWhile[# + 1 &, 888, !PrimeQ[#] &]
Out[8]= 907

In[9]:= NestWhile[# + 1 &, 888, !PrimeQ[#1] || !PrimeQ[#3] &, 3]
Out[9]= 1021

In[10]:= NestWhile[Mod[# + 3, 7] &, 0, UnsameQ, All]
Out[10]= 0

In[11]:= NestWhile[If[EvenQ[#], #/2, 3 # + 1] &, 27, # != 1 &]
Out[11]= 1
```

## NestWhileList
Like `NestWhile`, but returns the full list of intermediate results.
- `NestWhileList[f, expr, test]`: Generates `{expr, f[expr], f[f[expr]], ...}`, continuing while `test` applied to the most recent result yields `True`.
- `NestWhileList[f, expr, test, m]`: Supplies the most recent `m` results (not wrapped in a list) as arguments to `test`.
- `NestWhileList[f, expr, test, All]`: Supplies every result so far as arguments to `test`.
- `NestWhileList[f, expr, test, {mmin, mmax}]`: Defers the first test until `mmin` results exist, then supplies up to `mmax` most-recent results.
- `NestWhileList[f, expr, test, m, max]`: Caps the number of `f` applications at `max` (may be `Infinity`).
- `NestWhileList[f, expr, test, m, max, n]`: Appends `n` additional applications of `f` to the list after the loop terminates.
- `NestWhileList[f, expr, test, m, max, -n]`: Drops the last `n` elements from the list.

**Features**:
- `Protected`.
- Results are listed in generation order, including the final element on which `test` yielded a non-`True` value (or the last element produced when `max` iterations were reached).
- If `test[expr]` does not yield `True` initially, the result is just `{expr}`.
- `NestWhileList[f, expr, UnsameQ, 2]` is equivalent to `FixedPointList[f, expr]`.
- `NestWhileList[f, expr, test, All]` is equivalent to `NestWhileList[f, expr, test, {1, Infinity}]`.
- `NestWhileList[f, expr, UnsameQ, All]` continues until a previously-seen value reappears, and the repeat is included as the last element of the list.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhileList` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Examples**:
```
In[1]:= NestWhileList[#/2 &, 123456, EvenQ]
Out[1]= {123456, 61728, 30864, 15432, 7716, 3858, 1929}

In[2]:= NestWhileList[Log, 100., # > 0 &]
Out[2]= {100.0, 4.60517, 1.52718, 0.423423, -0.859384}

In[3]:= NestWhileList[Floor[#/2] &, 20, UnsameQ, 2, 4]
Out[3]= {20, 10, 5, 2, 1}

In[4]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[4]= {20, 10, 5, 2, 1}

In[5]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[5]= {20, 10, 5, 2, 1, 0}

In[6]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[6]= {20, 10, 5, 2}

In[7]:= NestWhileList[# + 1 &, 899, !PrimeQ[#] &]
Out[7]= {899, 900, 901, 902, 903, 904, 905, 906, 907}

In[8]:= NestWhileList[Mod[2 #, 19] &, 2, # != 1 &]
Out[8]= {2, 4, 8, 16, 13, 7, 14, 9, 18, 17, 15, 11, 3, 6, 12, 5, 10, 1}

In[9]:= NestWhileList[Mod[5 #, 7] &, 4, Unequal, All]
Out[9]= {4, 6, 2, 3, 1, 5, 4}

In[10]:= NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All]
Out[10]= {400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2}

In[11]:= NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All, Infinity, -1]
Out[11]= {400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1}
```

## FixedPointList
Generates the list of successive iterates of a function until a fixed point is reached.
- `FixedPointList[f, expr]`: Generates `{expr, f[expr], f[f[expr]], ...}`, continuing until two consecutive results are `SameQ`. The returned list always begins with `expr`, and its last two elements are always equal (the fixed point appears twice).
- `FixedPointList[f, expr, n]`: Stops after at most `n` applications of `f`. If `n` is reached before convergence, the last two elements need not be equal.
- `FixedPointList[f, expr, SameTest -> s]`: Uses the binary predicate `s` instead of `SameQ` to test successive pairs of results. Iteration stops when `s[prev, current]` evaluates to `True`.
- `FixedPointList[f, expr, n, SameTest -> s]`: Combines a bounded iteration count with a custom equivalence test.

**Features**:
- `Protected`.
- `FixedPointList[f, expr]` is equivalent to `NestWhileList[f, expr, UnsameQ, 2]`.
- `n` (when given) must be a non-negative integer or `Infinity`. Malformed argument specs leave `FixedPointList` unevaluated.

```mathematica
In[1]:= FixedPointList[1 + Floor[#/2] &, 1000]
Out[1]= {1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2}

In[2]:= 1 + Floor[Last[%]/2]
Out[2]= 2

In[3]:= FixedPointList[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]
Out[3]= {{28, 21}, {21, 7}, {7, 0}, {7, 0}}

In[4]:= GCD[28, 21]
Out[4]= 7

In[5]:= FixedPointList[1 + Floor[#/2] &, 1000, 5]
Out[5]= {1000, 501, 251, 126, 64, 33}

In[6]:= FixedPointList[(# + 2/#)/2 &, 1.0]
Out[6]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}

In[7]:= FixedPointList[(# + 2/#)/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.01 &)]
Out[7]= {1.0, 1.5, 1.41667, 1.41422}
```

## Fold
Successively applies a binary function to an accumulating seed and the elements of a list.
- `Fold[f, x, list]`: Returns the last element of `FoldList[f, x, list]`, namely `f[...f[f[f[x, list[[1]]], list[[2]]], list[[3]]]..., list[[n]]]`.
- `Fold[f, list]`: Equivalent to `Fold[f, First[list], Rest[list]]`.

**Features**:
- `Protected`.
- The head of the third argument need not be `List` (any compound expression is accepted).
- `Fold[f, x, {}]` returns `x` (the function is never applied); `Fold[f, {a}]` returns `a`.
- `Fold[f, {}]` remains unevaluated (no seed, no elements).
- Each intermediate application is evaluated before the next one.

```mathematica
In[1]:= Fold[f, x, {a, b, c, d}]
Out[1]= f[f[f[f[x, a], b], c], d]

In[2]:= Fold[List, x, {a, b, c, d}]
Out[2]= {{{{x, a}, b}, c}, d}

In[3]:= Fold[Times, 1, {a, b, c, d}]
Out[3]= a b c d

In[4]:= Fold[f, {a, b, c, d}]
Out[4]= f[f[f[a, b], c], d]

In[5]:= Fold[{2 #1, 3 #2} &, x, {a, b, c, d}]
Out[5]= {{{{16 x, 24 a}, 12 b}, 6 c}, 3 d}

In[6]:= Fold[f, x, p[a, b, c, d]]
Out[6]= f[f[f[f[x, a], b], c], d]

(* Horner-form evaluation: 2^3 + 2 + 1 at x = 2 *)
In[7]:= Fold[2 #1 + #2 &, 0, {1, 0, 1, 1}]
Out[7]= 11

(* Form a number from digits *)
In[8]:= Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= 451678

(* Form a continued fraction *)
In[9]:= Fold[1/(#2 + #1) &, x, Reverse[{a, b, c, d}]]
Out[9]= 1/(a + 1/(b + 1/(c + 1/(d + x))))

(* Factorial via Times *)
In[10]:= Fold[Times, 1, Range[5]]
Out[10]= 120

(* When the pure function ignores its second argument, Fold coincides with Nest *)
In[11]:= Fold[f[#1] &, x, Range[5]]
Out[11]= f[f[f[f[f[x]]]]]

(* Edge cases *)
In[12]:= Fold[f, x, {}]
Out[12]= x

In[13]:= Fold[f, {a}]
Out[13]= a

In[14]:= Fold[f, {}]
Out[14]= Fold[f, {}]
```

## FoldList
Produces the list of intermediate fold values.
- `FoldList[f, x, list]`: Gives `{x, f[x, list[[1]]], f[f[x, list[[1]]], list[[2]]], ...}`.
- `FoldList[f, list]`: Gives `{list[[1]], f[list[[1]], list[[2]]], ...}`.

**Features**:
- `Protected`.
- With a length-`n` list, `FoldList` returns a list of length `n + 1` (or `n` in the no-seed form).
- The head of the third argument is preserved in the output: `FoldList[f, x, p[a, b]]` gives `p[x, f[x, a], f[f[x, a], b]]`.
- `FoldList[f, {}]` returns `{}` (an empty list with the input head); `FoldList[f, x, {}]` returns `{x}`.
- `Fold[f, x, list]` is equivalent to `Last[FoldList[f, x, list]]`.

```mathematica
In[1]:= FoldList[f, x, {a, b, c, d}]
Out[1]= {x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]}

In[2]:= FoldList[f, {a, b, c, d}]
Out[2]= {a, f[a, b], f[f[a, b], c], f[f[f[a, b], c], d]}

(* Cumulative sums *)
In[3]:= FoldList[Plus, 0, Range[5]]
Out[3]= {0, 1, 3, 6, 10, 15}

(* Head preservation *)
In[4]:= FoldList[f, x, p[a, b, c, d]]
Out[4]= p[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]]

(* Fold to the right *)
In[5]:= FoldList[g[#2, #1] &, x, {a, b, c, d}]
Out[5]= {x, g[a, x], g[b, g[a, x]], g[c, g[b, g[a, x]]], g[d, g[c, g[b, g[a, x]]]]}

(* Successive factorials *)
In[6]:= FoldList[Times, 1, Range[10]]
Out[6]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800}

(* Build up a continued fraction *)
In[7]:= FoldList[1/(#2 + #1) &, x, Reverse[{a, b, c}]]
Out[7]= {x, 1/(c + x), 1/(b + 1/(c + x)), 1/(a + 1/(b + 1/(c + x)))}

(* Build up a number from digits *)
In[8]:= FoldList[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= {0, 4, 45, 451, 4516, 45167, 451678}

(* Horner-form polynomial evaluation at x = 2, coefficients {1, 0, 1, 1} *)
In[9]:= FoldList[2 #1 + #2 &, 0, {1, 0, 1, 1}]
Out[9]= {0, 1, 2, 5, 11}

(* Edge cases *)
In[10]:= FoldList[f, x, {}]
Out[10]= {x}

In[11]:= FoldList[f, {}]
Out[11]= {}

In[12]:= FoldList[f, p[]]
Out[12]= p[]
```

## Through
Distributes operators that appear inside the heads of expressions.
- `Through[expr]`: Distributes the top-level head.
- `Through[expr, h]`: Performs the transformation wherever `h` occurs in the head of `expr`.

**Features**:
- `Protected`.

```mathematica
In[1]:= Through[{f, g, h}[x]]
Out[1]= {f[x], g[x], h[x]}

In[2]:= Through[(f + g)[x, y]]
Out[2]= f[x, y] + g[x, y]
```

## Thread
Threads a function over the arguments that are lists (or have a given head).
- `Thread[f[args]]`: Threads `f` over any lists that appear in `args`.
- `Thread[f[args], h]`: Threads `f` over any objects with head `h` that appear in `args`.
- `Thread[f[args], h, n]`: Threads `f` over objects with head `h` that appear in the positions specified by `n`.

**Features**:
- `Protected`.
- Functions with attribute `Listable` are automatically threaded over lists.
- All the elements in the specified args whose heads are `h` must be of the same length; otherwise the expression is returned unchanged.
- Arguments that do not have head `h` are copied as many times as there are elements in the arguments that do have head `h`.
- The position specifier `n` uses the standard sequence specification:

  | Spec | Meaning |
  |------|---------|
  | `All` | all elements |
  | `None` | no elements |
  | `n` | elements 1 through `n` |
  | `-n` | last `n` elements |
  | `{n}` | element `n` only |
  | `{m, n}` | elements `m` through `n` inclusive |
  | `{m, n, s}` | elements `m` through `n` in steps of `s` |

```mathematica
In[1]:= Thread[f[{a, b, c}]]
Out[1]= {f[a], f[b], f[c]}

In[2]:= Thread[f[{a, b, c}, x]]
Out[2]= {f[a, x], f[b, x], f[c, x]}

In[3]:= Thread[f[{a, b, c}, {x, y, z}]]
Out[3]= {f[a, x], f[b, y], f[c, z]}

(* Convert equations for lists to lists of equations *)
In[4]:= Thread[{a, b, c} == {x, y, z}]
Out[4]= {a == x, b == y, c == z}

(* Apply a function to both sides of an equation *)
In[5]:= Thread[Log[x == y], Equal]
Out[5]= Log[x] == Log[y]

(* Thread over all arguments (the default) *)
In[6]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List]
Out[6]= {f[a, r, u, x], f[b, s, v, y]}

(* Do not thread at all *)
In[7]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, None]
Out[7]= f[{a, b}, {r, s}, {u, v}, {x, y}]

(* Thread over the first two arguments only *)
In[8]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, 2]
Out[8]= {f[a, r, {u, v}, {x, y}], f[b, s, {u, v}, {x, y}]}

(* Thread over the last two arguments only *)
In[9]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, -2]
Out[9]= {f[{a, b}, {r, s}, u, x], f[{a, b}, {r, s}, v, y]}

(* Thread over arguments 2 through 4 *)
In[10]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {2, 4}]
Out[10]= {f[{a, b}, r, u, x], f[{a, b}, s, v, y]}

(* Thread over every other argument *)
In[11]:= Thread[f[{a, b}, {r, s}, {u, v}, {x, y}], List, {1, -1, 2}]
Out[11]= {f[a, {r, s}, u, {x, y}], f[b, {r, s}, v, {x, y}]}

(* By default, does not thread over heads other than List *)
In[12]:= Thread[f[a + b, c + d]]
Out[12]= f[a + b, c + d]

(* Thread with respect to Plus *)
In[13]:= Thread[f[a + b, c + d], Plus]
Out[13]= f[a, c] + f[b, d]

(* Elements that are not lists are repeated *)
In[14]:= Thread[f[{a, b, c}, h, {x, y, z}]]
Out[14]= {f[a, h, x], f[b, h, y], f[c, h, z]}
```

## Select
Filters elements from an expression matching a criterion.
- `Select[list, crit]`: Returns an expression with the same head as `list`, containing only those elements `e` for which `crit[e]` evaluates to `True`.
- `Select[list, crit, n]`: Returns only the first `n` matching elements.

```mathematica
In[1]:= Select[{1, 2, 4, 7, 6, 2}, EvenQ]
Out[1]= {2, 4, 6, 2}

In[2]:= Select[{1, 2, 4, 7, 6, 2}, # > 2 &, 1]
Out[2]= {4}
```


## AllTrue, AnyTrue, NoneTrue
Test a predicate across the elements of a list (or the values of an association).
- `AllTrue[list, test]`: `True` if `test[e]` is `True` for every element (and for an empty list).
- `AnyTrue[list, test]`: `True` if `test[e]` is `True` for some element.
- `NoneTrue[list, test]`: `True` if `test[e]` is `True` for no element.

Left unevaluated when a test result is neither `True` nor `False`.

```mathematica
In[1]:= AllTrue[{2, 4, 6}, EvenQ]
Out[1]= True

In[2]:= AnyTrue[{1, 3, 4}, EvenQ]
Out[2]= True

In[3]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[3]= True
```
