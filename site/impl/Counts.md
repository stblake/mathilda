---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_counts` returns `<|element -> count, ...|>` in order of
first appearance. For a boxed `List` it makes a single hash pass: each element is
looked up in a `KeyIndex`, a new distinct value is recorded with count 1 and
every repeat increments the stored counter, after which the distinct keys and
their counts are emitted as rules. An `Association` argument is counted over its
values (`Counts[Values[assoc]]`).

**Data structures.** The distinct-value table is a `KeyIndex` open-addressing
hash set over borrowed `Expr*` keys paired with an `int64_t` count array; the
result is a canonical `Association` built directly from the rule array.

**Complexity / limits.** `O(n)` for a list of `n` elements. A packed/`NDArray`
argument takes the fast path: `counts_from_ndarray` runs `ndred_tally`
(`ndreduce.c`) over the raw `int64`/`float64` words — direct-indexed when the
value range allows, hashed otherwise — and relabels its `{key, count}` pairs as
`key -> count` rules, so the per-element work happens on machine words rather
than boxed `Expr`s. Tally declines the dtypes and ranks it cannot key faithfully
(complex, rank > 1, non-finite floats) and returns the List-path answer, which is
rewritten the same way.
