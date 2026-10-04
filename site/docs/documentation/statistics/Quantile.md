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

## Examples (11)

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

### Applications (6)

The default Type-1 quantile selects an element, not an average

```mathematica
In[6]:= Quantile[{1, 2, 3, 4}, 1/2]
Out[6]= 2
```

While Median averages the two central values, so the two deliberately differ

```mathematica
In[7]:= Median[{1, 2, 3, 4}]
Out[7]= 5/2
```

A list of probabilities gives a list of quantiles

```mathematica
In[8]:= Quantile[{1, 2, 3, 4}, {1/4, 1/2, 3/4}]
Out[8]= {1, 2, 3}
```

The Quartiles parameterization interpolates — here it is the median

```mathematica
In[9]:= Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]
Out[9]= 5/2
```

Data need not be pre-sorted; a machine-real sample gives a machine-real quantile

```mathematica
In[10]:= Quantile[{5., 1., 4., 2., 3.}, 0.9]
Out[10]= 5.0
```

The list form on real data

```mathematica
In[11]:= Quantile[{1., 2., 3., 4., 5.}, {0.25, 0.5, 0.75}]
Out[11]= {2.0, 3.0, 4.0}
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

**Algorithm.** `builtin_quantile` takes `Quantile[data, q]`, `Quantile[data,
{q1, …}]`, or `Quantile[data, q, {{a,b},{c,d}}]`. A visible `NDArray` / packed
argument is materialised to the exact `List` path with `pack_unpack`
(correctness-first — no buffer kernel yet; the audit baselines carry the declared
reason). `data` must be a non-empty `List`.

1. **Matrix** (first element is itself a `List`) — transpose and recurse
   columnwise, carrying `q` and the parameter matrix, so each column is quantiled
   independently.
2. **Validation** — every element must be real-numeric (`Quantile::rectn`
   otherwise, including an already-evaluated `Complex[re, im]` with nonzero `im`).
3. **Parameters** — `parse_param_matrix` validates a `{{a,b},{c,d}}` spec (head and
   shape checked, all four entries real-numeric); the default is `{{0,0},{1,0}}`,
   Hyndman–Fan **Type 1**, the left-continuous inverse CDF.
4. **Sort once** with `pack_eval_plain` (so a large machine-number sort's packed
   result is read through `.args`).
5. Each `q` goes through `quantile_one` → `stats_quantile_point` (shared with
   `Quartiles`, in `stats_common.c`). An exact-irrational `q` (e.g. `1/Sqrt[2]`) is
   read numerically through `N[]`; `q ∉ [0, 1]` raises `Quantile::q100`.

`stats_quantile_point` computes `h = a + (n+b)q`, edge-clamps (`h ≤ 1` → first
element, `h ≥ n` → last), takes `j = Floor[h]` clamped to `[1, n−1]`, `g = h − j`,
and the weight `w = c + d·g`. The upper index is `j` when `g = 0` (the two
neighbours coincide at integer `h`) else `j+1`. At `w = 0` or `w = 1` the element
is **selected** outright. For `w ∈ [0, 1]` it returns the convex combination
`(1−w)·A[j] + w·A[j+1]`; outside `[0, 1]` (which the parameterization permits,
though no standard type uses it) it returns the difference form `A[j] + w·(A[j+1] −
A[j])`. Each form is used only where it is numerically safe — the convex form
overflows on neighbours of opposite huge magnitude, the difference form on `w`
outside the unit interval. All arithmetic runs in the evaluator (exact in, exact
out); `double`s are read only for the clamp, index, and weight-regime decisions.

**Data structures.** `Expr` trees through `eval_and_free`; one sorted copy of the
data; an `Expr**` of per-`q` results for the list form. The four parameters are
kept as fresh `Expr*` copies freed on every exit path.

**Complexity / limits.** `O(n log n)` to sort, then `O(1)` per quantile point.
`Median` deliberately differs from `Quantile[…, 1/2]` on even-length data (`Median`
averages the two central order statistics; the default `Quantile` type selects one).
No `Compile[]` lowering and no packed kernel — a visible `NDArray` is unpacked
first.

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

- R. J. Hyndman and Y. Fan, *Sample quantiles in statistical packages*, The American Statistician **50** (1996) 361–365.
- Source: [`src/stats/quantile.c`](https://github.com/stblake/mathilda/blob/main/src/stats/quantile.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)

## Notes & additional examples

### Notes

The default parameters are `{{0, 0}, {1, 0}}` — Hyndman–Fan **Type 1**, the
left-continuous inverse CDF — so `Quantile[{1, 2, 3, 4}, 1/2]` selects the order
statistic `x₍⌈nq⌉₎ = 2`, whereas `Median` averages the two central values to get
`5/2`. The two are meant to differ on even-length data. `Quartiles[data]` is
`Quantile[data, {1/4, 1/2, 3/4}, {{1/2, 0}, {0, 1}}]`, which is why the fourth
example reproduces the median.

With `h = a + (n + b)q` and weight `w = c + d (h − ⌊h⌋)`, the result is the
edge-clamped interpolation between `x₍⌊h⌋₎` and `x₍⌈h⌉₎`; for `w ∈ [0, 1]` it is
evaluated as the convex combination `(1 − w) x₍⌊h⌋₎ + w x₍⌈h⌉₎`, which (unlike the
algebraically equal difference form) does not overflow on `Real` neighbours that
straddle zero near the machine range. `Quantile` requires real numeric data and
`q ∈ [0, 1]` (`Quantile::q100` otherwise), computes per column for a matrix, and
materialises a visible `NDArray` to the exact `List` path.
