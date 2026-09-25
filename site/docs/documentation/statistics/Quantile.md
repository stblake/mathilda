# Quantile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Quantile[data,q]`**

gives the q-th quantile estimate of the elements in data.

**`Quantile[data,{q1,q2,...}]`**

gives a list of quantile estimates.

**`Quantile[data,q,{{a,b},{c,d}}]`**

uses the quantile definition specified by parameters a, b, c, d.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Quantile[{1, 2, 3, 4}, 1/2]
Out[1]= 2

In[2]:= Quantile[{1, 2, 3, 4}, {1/4, 3/4}]
Out[2]= {1, 3}

In[3]:= Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]
Out[3]= 5/2

In[4]:= Quantile[{3, 1, 4, 2}, 1/2]
Out[4]= 2

In[5]:= Quantile[{{1, 2}, {3, 4}}, 1/2]
Out[5]= {1, 2}
```

## Algorithm

quantile.c -- Quantile[].

```text
Quantile[data, q]                 -- Wolfram default parameters {{0,0},{1,0}}:
                                     left-continuous, sorted[[Ceiling[n q]]].
Quantile[data, {q1, q2, ...}]     -- list of quantiles, one result per q.
Quantile[data, q, {{a,b},{c,d}}]  -- the general parameterized definition
                                     (same form Quartiles accepts; Quartiles
                                     is this with q = {1/4,1/2,3/4} and
                                     parameters {{1/2,0},{0,1}}).
```

The interpolation engine is shared with Quartiles: stats_quantile_point in stats_common.c. Matrix input recurses columnwise; a visible NDArray argument is materialised to the exact List path via pack_unpack (correctness-first -- no ndreduce kernel yet; the audit baselines carry the declared reason). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

- `Protected`.
- Default parameters are `{{0, 0}, {1, 0}}` (Wolfram's Type-1 / left-continuous inverse CDF), so `Quantile[{1, 2, 3, 4}, 1/2]` is `2` while `Median[{1, 2, 3, 4}]` is `5/2` — the two deliberately differ on even-length data.
- `Quartiles[data]` is `Quantile[data, {1/4, 1/2, 3/4}, {{1/2, 0}, {0, 1}}]`.
- Exact input gives exact output; sorting uses the canonical `Sort` order.
- For $w \in [0, 1]$ the interpolation is evaluated as the convex combination $(1-w)x_{(j)} + w\,x_{(j+1)}$, not as $x_{(j)} + w\,(x_{(j+1)} - x_{(j)})$. The two are equal in exact arithmetic, but the difference form overflows on `Real` data whose neighbours straddle zero near the machine range: `Quantile[{-1.0*10^308, 1.0*10^308}, 1/2, {{1/2, 0}, {0, 1}}]` is `0.0`, not `Infinity`. Outside $[0, 1]$ — which `{{a,b},{c,d}}` permits, though no standard quantile type uses it — the difference form is used instead, because there the convex form is the one that can produce `NaN`. At $w = 0$ or $w = 1$ the element is selected outright rather than computed.
- For `MatrixQ` data the quantile is computed per column.
- `Quantile` requires REAL numeric data and numeric $q \in [0, 1]$ (`Quantile::q100` otherwise); symbolic arguments stay unevaluated. A `Complex[re, im]` element is rejected with `Quantile::rectn` when `im` is nonzero — including an already-evaluated complex such as `2 + I` (`Complex[2, 1]`), which carries no literal `I` to search for — and accepted when `im` is zero, since `Complex[x, 0]` at MPFR precision is a real number. Not yet caught: a complex value nested under a numeric head, such as `Sqrt[2 + I]`.
- An `NDArray` argument is materialised to the exact `List` path (no buffer fast path yet).

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [Real](../../other-advanced/Real/), [MatrixQ](../../expression-information/MatrixQ/), [I](../../mathematical-constants/I/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)
