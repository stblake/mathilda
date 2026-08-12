# Moment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Moment[data, r]`**

gives the r-th raw (power) moment of data, (1/n) Sum\[x\_i^r\].

**`Moment[data, {r_1, ..., r_m}]`**

gives the multivariate raw moment of data. For a matrix or array the moment is taken columnwise over the first axis.

## Examples (5)

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

- `NHoldAll`, `Protected`.
- The raw moment is `CentralMoment` without the mean subtraction; `Moment[data, 1]` is `Mean[data]`, and `Moment[data, 0]` is `1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[Moment[#, r]&, x, 1]`); because there is no mean to subtract, `Mean[data^r]` threads correctly at every rank.
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `NHoldAll`, `Protected`.

## See also

[CentralMoment](../../statistics/CentralMoment/), [NDArray](../../linear-algebra/NDArray/), [Rational](../../arithmetic/Rational/), [List](../../other-advanced/List/), [Variance](../../data-structures/Variance/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
