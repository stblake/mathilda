---
source: src/partitions.c
---
**Algorithm.** `builtin_integerpartitions` collapses all of its call forms —
`IntegerPartitions[n]`, `[n, k]`, `[n, {k}]`, `[n, {kmin, kmax}]`, `[n, {kmin, kmax, dk}]`,
an optional parts set `sspec`, and a final count `m` — onto one recursive count-vector
enumerator (`ip_recurse`) over an ordered, reversed set of allowed parts. The default
`sspec` is `Range[n]` reversed (`floor(n), …, 1`); the recursion assigns a count to each
part in descending order, pruned by a length budget from a finite `kmax` and a tight numeric
bound (`ip_floor_div`) that is valid only when the remaining parts are single-signed.
Complete partitions — `remaining == 0` with a length inside `[kmin, kmax]` stepped by `dk` —
are emitted in reverse-lexicographic order, and a final `m` slices the first `m` (or, when
negative, the last `|m|`) of them.

**Data structures.** Everything runs in exact GMP rationals (`mpq_t`), so `n` and each
`s_i` may be an integer, big integer, rational or negative value handled by one code path.
The enumeration context `PartCtx` holds the reversed parts array, precomputed suffix-sign
flags (which license the numeric bound), the current count vector, and a growable `Expr**`
of emitted `List` partitions. The builtin only reads `res`; a `NULL` return (bad arguments,
symbolic input, infinite result) leaves the call unevaluated. There is no ND/packed/`Compile`
path — it is a structural list generator, and `IntegerPartitions` is `Protected` but not
`Listable`.

**Complexity / limits.** Output-sensitive: proportional to the number of partitions
enumerated, times `O(nparts)` per emission. An infinite result (a part of `0`, or mixed
signs with an unbounded part count) is detected up front and raises `IntegerPartitions::undef`;
a `take` count exceeding the number found warns `IntegerPartitions::take`; `0` or more than
`4` arguments emit `IntegerPartitions::argb`. Symbolic or real `n` / `s_i` leave the call
unevaluated.
