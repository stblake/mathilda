---
references:
  - "R. J. Hyndman and Y. Fan, *Sample quantiles in statistical packages*, The American Statistician **50** (1996) 361–365."
source: src/stats/interquartilerange.c
---
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
