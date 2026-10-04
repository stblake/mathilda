# CentralMoment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CentralMoment[data, r]`**

gives the r-th central moment (moment about the mean) of data, (1/n) Sum\[(x\_i - Mean\[data\])^r\].

**`CentralMoment[data, {r_1, ..., r_m}]`**

gives the multivariate central moment of data. For a matrix or array the moment is taken columnwise over the first axis.

## Examples (11)

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

### Applications (7)

The second central moment — Variance without the bias correction

```mathematica
In[5]:= CentralMoment[{2, 4, 4, 4, 5, 5, 7, 9}, 2]
Out[5]= 4
```

Compare: Variance divides by n - 1, not n

```mathematica
In[6]:= Variance[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[6]= 32/7
```

A machine-real sample gives a machine-real moment

```mathematica
In[7]:= CentralMoment[{1., 2., 3., 4., 5.}, 2]
Out[7]= 2.0
```

The first central moment is always zero

```mathematica
In[8]:= CentralMoment[{2, 4, 4, 4, 5, 5, 7, 9}, 1]
Out[8]= 0
```

An odd central moment vanishes on symmetric data

```mathematica
In[9]:= CentralMoment[{1, 2, 3, 4}, 3]
Out[9]= 0
```

On a matrix the moment is columnwise over the first axis

```mathematica
In[10]:= CentralMoment[{{1, 2}, {3, 4}, {5, 6}}, 2]
Out[10]= {8/3, 8/3}
```

A list order gives the multivariate mixed central moment

```mathematica
In[11]:= CentralMoment[{{0, 1}, {2, 5}, {4, 3}}, {1, 1}]
Out[11]= 4/3
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

**Algorithm.** `builtin_central_moment` takes `CentralMoment[data, order]`, the
`order`-th moment about the mean, `μ̃_r = (1/n) Σ (xᵢ − Mean[data])^r`. An
`Association` is handled over its values; an `NDArray` / packed argument takes the
buffer fast path `ndred_central_moment`. Otherwise `data` must be a non-empty
`List` and dispatch is by shape:

1. **List order `{r1, …, rm}`** (`cm_multivariate`) — the multivariate mixed
   central moment. With `mu = Mean[data]` (rank `k−1`), each first-axis slice
   contributes `Times @@ (sub − mu)^rvec`, then `Mean` over the slices. The slice's
   second-axis length must equal `Length[rvec]`.
2. **Scalar order over an array of depth ≥ 2** (`cm_columnwise`) — reduces the
   first axis columnwise: `mu = Mean[data]`, and each slice gives `(sub − mu)^r`,
   `Mean`-averaged over the slices. The subtraction is done *per slice* because the
   obvious `data − Mean[data]` would thread row-wise, not column-wise.
3. **Scalar order over a flat vector** — fast real path when every element is a
   real number, at least one is inexact, and `order` is a non-negative machine
   integer: one pass computes the mean, a second accumulates `Σ (xᵢ − mean)^r` by
   integer power-by-squaring (`cm_ipow`), returning `acc / n` as a machine `Real`.
4. **Everything else** (`cm_vector_symbolic`) — `Mean[(data − Mean[data])^r]`,
   where the `Listable` `Plus` and `Power` thread over the vector and the outer
   `Mean` divides by `n`, reusing `Mean`'s exact and symbolic paths.

The design mirrors `Variance` with the `n/(n−1)` bias correction removed: it
divides by `n` (not `n−1`), raises to `r` (not a square), needs only `n ≥ 1`, and
applies no `Conjugate` — a central moment is `(x − μ)^r`, not `|x − μ|²`.
`CentralMoment[data, 1]` is `0`.

**Data structures.** `Expr` trees via `eval_and_free`; the columnwise and
multivariate cases hold an `Expr**` of `n` per-slice terms. The fast path uses two
`double` scalar passes.

**Complexity / limits.** `O(n log r)` machine-double, two passes; the symbolic path
is `O(n)` tree builds plus `Mean`. `packed_aware` in `src/pack.c` but not
`INT64_OK` (an integer sample's central moment is a `Rational`, so an integer
buffer degrades to the exact `List` path). Lowers inside `Compile[]` at a real
vector with integer order; participates in auto-compilation.

- `Protected`.
- A central moment is `Variance` without the $n/(n-1)$ bias correction: it divides by `n` (not `n-1`), raises to the power `r` (not a square), and needs only `n >= 1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[CentralMoment[#, r]&, x, 1]`).
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_central_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

**Attributes:** `Protected`.

## References

**See also:** [Variance](../../data-structures/Variance/), [NDArray](../../linear-algebra/NDArray/), [Rational](../../arithmetic/Rational/), [List](../../other-advanced/List/)

- M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), ch. 3 (moments about the mean).
- Source: [`src/stats/central_moment.c`](https://github.com/stblake/mathilda/blob/main/src/stats/central_moment.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`CentralMoment[data, r] = (1/n) Σ (xᵢ − Mean[data])^r`. It is `Variance` with the
bias correction removed: it divides by `n` rather than `n − 1`, raises to the power
`r` rather than squaring, needs only `n ≥ 1`, and applies no `Conjugate`. So
`CentralMoment[data, 2]` and `Variance[data]` differ by exactly the factor
`(n − 1)/n`.

For a matrix or array the centering is done per first-axis slice, because the
obvious `data − Mean[data]` would thread row-wise rather than column-wise.

A real `NDArray` / packed buffer takes the kernel `ndred_central_moment` and the
head lowers inside `Compile[]`; an integer buffer degrades to the exact `List`
path, since the exact answer is a `Rational`.
