# Moment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Moment[data, r]`**

gives the r-th raw (power) moment of data, (1/n) Sum\[x\_i^r\].

**`Moment[data, {r_1, ..., r_m}]`**

gives the multivariate raw moment of data. For a matrix or array the moment is taken columnwise over the first axis.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Moment[{1, 2, 3, 4}, 2]
Out[1]= 15/2

In[2]:= Moment[{1., 2., 3., 4.}, 2]
Out[2]= 7.5

In[3]:= Moment[{Pi, E, 2}, 1]
Out[3]= 1/3 (2 + E + Pi)

In[4]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, 3]
Out[4]= {51, 96}

In[5]:= Simplify[Moment[{{a, b}, {c, d}}, {1, 2}]]
Out[5]= 1/2 (a b^2 + c d^2)
```

### Applications (7)

The raw second moment of an integer sample is exact

```mathematica
In[6]:= Moment[{2, 4, 4, 4, 5, 5, 7, 9}, 2]
Out[6]= 29
```

A machine-real sample gives a machine-real moment

```mathematica
In[7]:= Moment[{1., 2., 3., 4.}, 3]
Out[7]= 25.0
```

The first raw moment is just the mean

```mathematica
In[8]:= Moment[{1, 2, 3, 4}, 1]
Out[8]= 5/2
```

And the zeroth raw moment is 1 for any data

```mathematica
In[9]:= Moment[{1, 2, 3, 4}, 0]
Out[9]= 1
```

On a matrix the moment is taken columnwise over the first axis

```mathematica
In[10]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, 2]
Out[10]= {35/3, 56/3}
```

Symbolic data stays symbolic

```mathematica
In[11]:= Moment[{a, b, c}, 2]
Out[11]= 1/3 (a^2 + b^2 + c^2)
```

A list order gives the multivariate mixed raw moment

```mathematica
In[12]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, {1, 2}]
Out[12]= 232/3
```

## Algorithm

moment.c -- Moment[] (raw / power moment). Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

```text
Moment[data, r]              — the r-th raw (power) moment,
                               mu_r = (1/n) Sum[x_i^r]. For a matrix / array the
                               reduction is columnwise over the first axis
                               (equivalently ArrayReduce[Moment[#,r]&, x, 1]).
Moment[data, {r1, ..., rm}]  — the multivariate mixed raw moment,
                               (1/n) Sum_i Product_j x[[i,j]]^r_j,
                               summing the first axis and taking a product over
                               the second (its length must equal Length[{r1,...}]).
```

The raw moment is CentralMoment without the mean subtraction. Because there is no mean to subtract, Mean[data^r] threads correctly for a vector, a matrix (columnwise), AND a higher-rank array in a single expression — Power is Listable so data^r threads elementwise at every rank, and the outer Mean collapses the first axis by n. So (unlike CentralMoment, whose data - Mean[data] would thread row-wise) the scalar-order case needs no separate columnwise routine. Numeric real vectors take a tight C loop (and packed / NDArray inputs a machine-buffer fast path via ndred_moment); every other case — exact, symbolic, matrix/array, multivariate — is built as an expression and handed to the evaluator, which already knows how to be exact or symbolic.

## Implementation notes

**Algorithm.** `builtin_moment` takes `Moment[data, order]`. An `Association`
argument is handled over its values (`assoc_apply_over_values`); an `NDArray` /
packed-array argument takes the buffer fast path `ndred_moment`. Otherwise `data`
must be a non-empty `List` and the dispatch is by the shape of `order`:

1. **List order `{r1, …, rm}`** (`mom_multivariate`) — the multivariate mixed raw
   moment. For each first-axis slice `sub` it forms `Times @@ sub^rvec` (the
   `Listable` `Power` raises block `j` to `r_j` along the second axis, and
   `Times @@` takes the product over that axis), then `Mean` over the slices. The
   slice's second-axis length must equal `Length[rvec]`.
2. **Scalar order over a flat vector** — fast real path: when every element is a
   real number, at least one is inexact, and `order` is a non-negative machine
   integer, a tight C loop accumulates `Σ xᵢ^r` with integer power-by-squaring
   (`mom_ipow`) and returns `acc / n` as a machine `Real`.
3. **Everything else** (`mom_symbolic`) — `Mean[data^r]`. Because the raw moment
   has no mean to subtract, `Power` being `Listable` makes `data^r` thread
   element-wise at *every* rank — so a single expression gives the vector's scalar
   moment, the matrix's columnwise vector, or the array's columnwise array, and no
   separate columnwise routine is needed (unlike `CentralMoment`). This reuses
   `Mean`'s exact-rational, symbolic, and MPFR paths.

`Moment[data, 1]` is `Mean[data]` and `Moment[data, 0]` is `1`, both falling out
of the same definition.

**Data structures.** `Expr` trees driven through `eval_and_free`; the multivariate
case holds an `Expr**` array of `n` per-slice terms before wrapping them in a
`List` for the outer `Mean`. The fast path is a single `double` accumulator.

**Complexity / limits.** `O(n log r)` for the machine-double loop (power by
squaring per element); the symbolic path is `O(n)` tree builds plus the
evaluator's `Mean` cost. Registered as `packed_aware` in `src/pack.c` but **not**
`INT64_OK` — an integer sample's raw moment is a `Rational` no `int64` slot holds,
so an integer buffer degrades to the exact `List` path, exactly like `Variance`.
Lowers inside `Compile[]` at a real vector with integer order, so it participates
in auto-compilation.

- `NHoldAll`, `Protected`.
- The raw moment is `CentralMoment` without the mean subtraction; `Moment[data, 1]` is `Mean[data]`, and `Moment[data, 0]` is `1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[Moment[#, r]&, x, 1]`); because there is no mean to subtract, `Mean[data^r]` threads correctly at every rank.
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `NHoldAll`, `Protected`.

## References

**See also:** [CentralMoment](../../statistics/CentralMoment/), [NDArray](../../linear-algebra/NDArray/), [Rational](../../arithmetic/Rational/), [List](../../other-advanced/List/), [Variance](../../data-structures/Variance/)

- M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), ch. 3 (moments).
- Source: [`src/stats/moment.c`](https://github.com/stblake/mathilda/blob/main/src/stats/moment.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

The raw (power) moment `Moment[data, r]` is `CentralMoment[data, r]` without the
mean subtraction: `μ_r = (1/n) Σ xᵢ^r`. Because there is no mean to subtract,
`Mean[data^r]` threads correctly at every rank, so the matrix and higher-rank
columnwise forms need no special handling — the `Listable` `Power` does the work.

A real `NDArray` / packed buffer takes the kernel `ndred_moment`, and the head
lowers inside `Compile[]` for a real vector with integer order. An integer buffer
falls back to the exact `List` path, because an integer sample's raw moment is a
`Rational` no machine slot can hold.

`Moment` carries the `NHoldAll` attribute (as in Mathematica), so `N` does not
thread into a symbolic order or symbolic data.
