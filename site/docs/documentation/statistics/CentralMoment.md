# CentralMoment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CentralMoment[data, r]`**

gives the r-th central moment (moment about the mean) of data, (1/n) Sum\[(x\_i - Mean\[data\])^r\].

**`CentralMoment[data, {r_1, ..., r_m}]`**

gives the multivariate central moment of data. For a matrix or array the moment is taken columnwise over the first axis.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= CentralMoment[{1, 2, 3, 4}, 4]
Out[1]= 41/16

In[2]:= CentralMoment[{1., 2., 3., 4.}, 2]
Out[2]= 1.25

In[3]:= CentralMoment[{{1, 2}, {3, 4}, {5, 6}}, 2]
Out[3]= {8/3, 8/3}

In[4]:= Simplify[CentralMoment[{{a, b}, {c, d}}, {2, 2}]]
Out[4]= 1/16 (a - c)^2 (b - d)^2
```

## Algorithm

central_moment.c -- CentralMoment[]. Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

```text
CentralMoment[data, r]              — the r-th moment about the mean,
                                      mu~_r = (1/n) Sum[(x_i - mu_1)^r], where
                                      mu_1 = Mean[data]. For a matrix / array the
                                      reduction is columnwise over the first axis
                                      (equivalently ArrayReduce[CentralMoment[#,r]&, x, 1]).
CentralMoment[data, {r1, ..., rm}]  — the multivariate mixed central moment,
                                      (1/n) Sum_i Product_j (x[[i,j]] - mu_1[[j]])^r_j,
                                      summing the first axis and taking a product over
                                      the second (its length must equal Length[{r1,...}]).
```

The design mirrors Variance (a central moment is Variance without the n/(n-1) bias correction): divide by n (not n-1), raise to the power r (not square), n >= 1 suffices, and there is no Conjugate — a central moment is (x-mu)^r, not

```text
|x-mu|^2. Numeric real vectors take a tight C loop (and packed / NDArray inputs
```

a machine-buffer fast path via ndred_central_moment); every other case — exact, symbolic, matrix/array, multivariate — is built as an expression and handed to the evaluator, which already knows how to be exact or symbolic.

## Implementation notes

- `Protected`.
- A central moment is `Variance` without the $n/(n-1)$ bias correction: it divides by `n` (not `n-1`), raises to the power `r` (not a square), and needs only `n >= 1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[CentralMoment[#, r]&, x, 1]`).
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_central_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `Protected`.

## See also

[Variance](../../data-structures/Variance/), [NDArray](../../linear-algebra/NDArray/), [Rational](../../arithmetic/Rational/), [List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
