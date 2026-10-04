# Kurtosis

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Kurtosis[data]`**

gives the coefficient of kurtosis (peak/tail vs flank concentration) of data, equivalent to CentralMoment\[data, 4\] / CentralMoment\[data, 2\]^2. For a matrix it is taken columnwise.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Kurtosis[{1, 2, 3, 4, 5}]
Out[1]= 17/10

In[2]:= Kurtosis[{1, 2, 4, 8}]
Out[2]= 25141/13225
```

### Applications (3)

The Pearson kurtosis of an evenly spaced integer sample

```mathematica
In[3]:= Kurtosis[{1, 2, 3, 4, 5}]
Out[3]= 17/10
```

Exact input gives an exact rational

```mathematica
In[4]:= Kurtosis[{1, 2, 4, 8}]
Out[4]= 25141/13225
```

A machine-real sample gives a machine-real coefficient

```mathematica
In[5]:= Kurtosis[{2., 4., 4., 4., 5., 5., 7., 9.}]
Out[5]= 2.78125
```

## Algorithm

kurtosis.c -- Kurtosis[]. Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

Kurtosis[data] -- the coefficient of kurtosis, a measure of peak/tail vs flank concentration. Equivalent to CentralMoment[data, 4] / CentralMoment[data, 2]^2 (Pearson kurtosis, not the excess form). For a matrix or array it is taken columnwise (the CentralMoment ratio threads). The shared body lives in stats_common.c (stats_standardized_moment), which also routes NDArray / packed inputs to the buffer kernel ndred_kurtosis.

## Implementation notes

**Algorithm.** `Kurtosis[data]` is the coefficient of kurtosis — a measure of how
concentrated the distribution is in the peak and tails versus the flanks — equal
to `CentralMoment[data, 4] / CentralMoment[data, 2]^2`. This is **Pearson**
kurtosis, not the excess form (which subtracts `3`). `builtin_kurtosis` is a
one-line call into the shared body `stats_standardized_moment(res, 4)` in
`src/stats/stats_common.c`, the same routine `Skewness` uses with `p = 3`:

1. an `Association` is handled over its values;
2. an `NDArray` / packed argument routes to the buffer kernel `ndred_kurtosis`;
3. otherwise `mp = CentralMoment[data, 4]` and `m2 = CentralMoment[data, 2]` are
   formed, and if either stays a `CentralMoment` head (data not reducible) the
   result is `NULL` and the caller's head is left unevaluated;
4. the result is `mp / m2^(4/2) = mp / m2²`.

`Power` and `Divide` thread, so a matrix yields the columnwise vector of kurtoses.
Exact input gives exact output (a `Rational` here, since both moments enter with
even powers).

**Data structures.** `Expr` trees through `eval_and_free`; the two central moments
are the only intermediates. The buffer path works on the packed `double` sample
inside `ndred_kurtosis`.

**Complexity / limits.** Two `CentralMoment` reductions (`O(n)` each) plus a
constant combine. `packed_aware` in `src/pack.c` but not `INT64_OK` (the integer
answer is a `Rational`, so an integer buffer degrades to the exact `List` path).
Lowers inside `Compile[]` (a `ND_REDS` delegate) and participates in
auto-compilation.

- `Protected`.
- Equivalent to `CentralMoment[data, 4] / CentralMoment[data, 2]^2` (Pearson kurtosis, not the excess form).
- For a matrix, gives the columnwise kurtoses.
- Handles numerical and symbolic data; exact input gives exact output.
- Fast path on `NDArray`/packed real buffers (`ndred_kurtosis`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

- M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), §3.31 (the standardized fourth moment).
- Source: [`src/stats/kurtosis.c`](https://github.com/stblake/mathilda/blob/main/src/stats/kurtosis.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`Kurtosis[data]` is the standardized fourth central moment,
`CentralMoment[data, 4] / CentralMoment[data, 2]^2`. This is **Pearson** kurtosis,
not the excess form — subtract `3` to obtain the excess kurtosis (which is `0` for a
normal distribution). The exact value is a `Rational`, since both moments enter with
even powers.

For a matrix the coefficient is taken columnwise. A real `NDArray` / packed buffer
takes the kernel `ndred_kurtosis` and the head lowers inside `Compile[]`; an integer
buffer degrades to the exact `List` path.
