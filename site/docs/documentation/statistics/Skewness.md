# Skewness

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Skewness[data]`**

gives the coefficient of skewness (a measure of asymmetry) of data, equivalent to CentralMoment\[data, 3\] / CentralMoment\[data, 2\]^(3/2). For a matrix it is taken columnwise.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Skewness[{1, 2, 3, 10}]
Out[1]= 18/25 Sqrt[2]

In[2]:= Skewness[{1., 2., 3., 4., 5.}]
Out[2]= 0.0
```

### Applications (3)

A right-skewed integer sample; the exact value is a radical

```mathematica
In[3]:= Skewness[{1, 2, 3, 10}]
Out[3]= 18/25 Sqrt[2]
```

A symmetric sample has zero skewness

```mathematica
In[4]:= Skewness[{1, 2, 3, 4, 5}]
Out[4]= 0
```

A machine-real sample gives a machine-real coefficient

```mathematica
In[5]:= Skewness[{2., 4., 4., 4., 5., 5., 7., 9.}]
Out[5]= 0.65625
```

## Algorithm

skewness.c -- Skewness[]. Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

Skewness[data] -- the coefficient of skewness, a measure of asymmetry. Equivalent to CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2). For a matrix or array it is taken columnwise (the CentralMoment ratio threads). The shared body lives in stats_common.c (stats_standardized_moment), which also routes NDArray / packed inputs to the buffer kernel ndred_skewness.

## Implementation notes

**Algorithm.** `Skewness[data]` is the coefficient of skewness, a measure of
asymmetry, equal to `CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)`.
`builtin_skewness` is a one-line call into the shared body
`stats_standardized_moment(res, 3)` in `src/stats/stats_common.c`, which `Kurtosis`
reuses with `p = 4`. That body:

1. handles an `Association` over its values (`assoc_apply_over_values`);
2. routes an `NDArray` / packed argument to the buffer kernel `ndred_skewness`;
3. otherwise forms `mp = CentralMoment[data, 3]` and `m2 = CentralMoment[data, 2]`
   and, if *either* comes back still wearing a `CentralMoment` head (the data did
   not reduce — a bare symbol, an empty list), returns `NULL` so the caller's head
   stays unevaluated;
4. returns `mp / m2^(3/2)`.

Because `Power` and `Divide` thread, a matrix's columnwise vector of central
moments yields a columnwise vector of skewnesses. Exact input gives exact output,
which is a radical in general (e.g. a `Sqrt[2]` multiple).

**Data structures.** `Expr` trees built and evaluated through `eval_and_free`; the
two central moments are the only intermediates. The buffer path works on the
packed `double` sample inside `ndred_skewness`.

**Complexity / limits.** Two `CentralMoment` reductions, each `O(n)`, plus a
constant-cost combine. `packed_aware` in `src/pack.c` but not `INT64_OK` — an
integer sample's skewness is a radical no machine slot holds, so an integer buffer
degrades to the exact `List` path. Lowers inside `Compile[]` (a `ND_REDS`
delegate) and participates in auto-compilation. A symmetric sample gives exactly
`0`.

- `Protected`.
- Equivalent to `CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)`.
- For a matrix, gives the columnwise skewnesses.
- Handles numerical and symbolic data; exact input gives exact output (a radical in general).
- Fast path on `NDArray`/packed real buffers (`ndred_skewness`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

- M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), §3.31 (the standardized third moment).
- Source: [`src/stats/skewness.c`](https://github.com/stblake/mathilda/blob/main/src/stats/skewness.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`Skewness[data]` is the standardized third central moment,
`CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)` — a dimensionless measure
of asymmetry. A positive value signals a longer right tail, a negative value a
longer left tail, and a symmetric distribution gives exactly `0`. Exact input
yields exact output, which is a radical in general.

For a matrix the coefficient is taken columnwise (the `CentralMoment` ratio
threads). A real `NDArray` / packed buffer takes the kernel `ndred_skewness` and
the head lowers inside `Compile[]`; an integer sample's skewness is a radical, so
an integer buffer degrades to the exact `List` path.
