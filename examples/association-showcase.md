# Association showcase: a split-apply-combine pipeline

This walks a realistic data-analysis pipeline over a synthetic transaction log
(100,000 `{category, amount}` records) using **only** Mathilda's Association
toolchain. Every heavy step is hash-backed and amortised `O(n)`.

Reproduce with:

```bash
./Mathilda < examples/association_showcase.m
```

## The pipeline

```mathematica
SeedRandom[42];
categories = {"food", "rent", "travel", "fun", "misc"};
txns = Table[{categories[[RandomInteger[{1, 5}]]], RandomInteger[{1, 100}]}, {100000}];
```

**1. Transactions per category** — a hash tally:

```mathematica
In[]:= Counts[txns[[All, 1]]]
Out[]= <|"rent" -> 20003, "misc" -> 19832, "fun" -> 20244, "travel" -> 19844, "food" -> 20077|>
```

**2. Total spend per category** — split-apply-combine: group by category, reduce
each group by summing its amounts:

```mathematica
In[]:= spend = GroupBy[txns, First, Total[#[[All, 2]]] &]
Out[]= <|"rent" -> 1011632, "misc" -> 996142, "fun" -> 1020869, "travel" -> 1000921, "food" -> 1013468|>
```

**3. Rank categories by spend** (descending), and take the top one:

```mathematica
In[]:= ReverseSort[spend]
Out[]= <|"fun" -> 1020869, "food" -> 1013468, "rent" -> 1011632, "travel" -> 1000921, "misc" -> 996142|>

In[]:= TakeLargest[spend, 1]
Out[]= <|"fun" -> 1020869|>
```

**4. Average transaction size per category** — a different reducer over the same
groups:

```mathematica
In[]:= N[GroupBy[txns, First, Mean[#[[All, 2]]] &]]
Out[]= <|"rent" -> 50.574, "misc" -> 50.229, "fun" -> 50.4282, "travel" -> 50.4395, "food" -> 50.4791|>
```

## Timings (100,000 records)

Measured on Apple Silicon (arm64), `-O3` build:

| Step                                   | Time    |
|----------------------------------------|---------|
| `Counts[txns[[All,1]]]`                | ~30 ms  |
| `GroupBy[txns, First, Total[...] &]`   | ~64 ms  |
| `ReverseSort[spend]`                   | ~0.004 ms |

`Counts` and `GroupBy` are the `O(n)` hash passes over all 100k records;
`GroupBy`'s extra cost is running the classifier `First` and the reducer once
per group. `ReverseSort` acts on the 5-entry summary, so it is effectively free.
The whole pipeline stays linear in the number of records — see
[`association-benchmarks.md`](association-benchmarks.md) for the scaling study.

## Toolchain used

A single, coherent data pipeline touched: `Table`/`Part` (column extraction via
`[[All, k]]`), `Counts`, `GroupBy` (with reducer), `Total`, `Mean`,
`ReverseSort`, `TakeLargest` — all operating on associations with Wolfram-faithful
semantics.

## Loadable notebook

A guided, four-notebook tour of the association toolchain is available as a
loadable library: [`associations.lb`](associations.lb) — open it in the Mathilda
app to explore basics, aggregation, transform/rank, and pattern matching interactively.
