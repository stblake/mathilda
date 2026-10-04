# InterquartileRange

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InterquartileRange[data]`**

gives the difference between the upper and lower quartiles of the elements in data.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= 4

In[2]:= InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]
Out[2]= {6, 6, 6}
```

### Applications (3)

The exact interquartile range of an integer sample

```mathematica
In[3]:= InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[3]= 4
```

A machine-real sample gives a machine-real range

```mathematica
In[4]:= InterquartileRange[{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}]
Out[4]= 5.0
```

On a matrix the IQR is computed per column

```mathematica
In[5]:= InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]
Out[5]= {6, 6, 6}
```

## Algorithm

interquartilerange.c -- InterquartileRange[].

InterquartileRange[data] = q3 - q1 of Quartiles[data] (Wolfram's IQR uses the Quartiles parameterization {{1/2,0},{0,1}}, so composing over the Quartiles builtin is the definition, not a shortcut). Matrix input recurses per-column BEFORE the vector path: a k-column matrix's Quartiles result is a k-list of triples, and for k == 3 a bare "is it a 3-list" guard would confuse it with a vector's {q1,q2,q3} and compute Quartiles(col3)-Quartiles(col1) -- the silent wrong answer the STATS-1 plan review flagged. The vector-path guard also requires all three quartiles to be SCALARS for the same reason. NDArray input is materialised via pack_unpack (correctness-first, no kernel). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

**Algorithm.** `InterquartileRange[data]` is `q̂₃ − q̂₁`, the difference between the
upper and lower quartiles. It is defined *by composition* over the `Quartiles`
builtin — Wolfram's IQR uses the `Quartiles` parameterization `{{1/2,0},{0,1}}`, so
calling `Quartiles` is the definition, not a shortcut. `builtin_interquartilerange`:

1. materialises a visible `NDArray` / packed argument via `pack_unpack`
   (correctness-first, no kernel yet);
2. requires `data` to be a non-empty `List`;
3. **handles the matrix case first** — if the first element is itself a `List`, it
   recurses columnwise via `stats_apply_columnwise("InterquartileRange", data)`.
   This ordering is deliberate: a 3-column matrix's `Quartiles` result is a 3-list
   of triples, and a bare "is it a 3-list?" test would mistake it for a vector's
   `{q1, q2, q3}` and compute `Quartiles(col₃) − Quartiles(col₁)` — a silent wrong
   answer the plan review flagged;
4. checks every element is real-numeric (`InterquartileRange::rectn` otherwise);
5. evaluates `Quartiles[data]`, requires it to be a 3-list of real-numeric
   **scalars** (the same guard against the matrix/vector confusion), and returns
   `q₃ − q₁`.

**Data structures.** `Expr` trees through `eval_and_free`; one `Quartiles` result
list. The matrix recursion reuses `stats_apply_columnwise`
(`Map[InterquartileRange, Transpose[matrix]]`).

**Complexity / limits.** Dominated by the single `Quartiles` call (`O(n log n)` to
sort). A robust scale estimator — insensitive to the tails beyond the quartiles.
No packed buffer kernel and no `Compile[]` lowering: a visible `NDArray` is
unpacked to the exact `List` path first.

- `Protected`.
- A robust scale estimator: insensitive to outliers beyond the quartiles.
- For `MatrixQ` data the IQR is computed per column.
- Requires numeric data (`InterquartileRange::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [Quartiles](../../statistics/Quartiles/), [MatrixQ](../../expression-information/MatrixQ/)

- R. J. Hyndman and Y. Fan, *Sample quantiles in statistical packages*, The American Statistician **50** (1996) 361–365.
- Source: [`src/stats/interquartilerange.c`](https://github.com/stblake/mathilda/blob/main/src/stats/interquartilerange.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)

## Notes & additional examples

### Notes

`InterquartileRange[data]` is `q̂₃ − q̂₁`, the spread between the upper and lower
quartiles, computed by composition over `Quartiles` (and so using the same
`{{1/2, 0}, {0, 1}}` parameterization). It is a robust scale estimator — insensitive
to the data beyond the quartiles — and so a good companion to `MedianDeviation`.

The matrix case is handled before the vector case on purpose: a three-column
matrix's per-column quartiles are a 3-list of triples, which a naive shape test
would confuse with a vector's `{q₁, q₂, q₃}`. The head requires numeric data
(`InterquartileRange::rectn` otherwise) and materialises a visible `NDArray` to the
exact `List` path.
