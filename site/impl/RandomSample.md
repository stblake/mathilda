---
source: src/random.c
---
**Algorithm.** `builtin_randomsample` (in `src/random.c`) selects elements *without replacement*. The uniform form uses `fisher_yates_sample(total, n)`: a partial Fisher–Yates shuffle of an index array `[0..total)` that performs only the first `n` swaps (each swap picks `j` uniformly from the remaining suffix via `random_index`), returning the first `n` shuffled indices. With no count it returns a full random permutation. The size argument may be an integer or `UpTo[n]` (clamped to the list length via `is_upto`).

The weighted form `RandomSample[{w1,...}->{e1,...}, n]` uses `weighted_sample_without_replacement`, which repeatedly draws by inverse-CDF over the live weights (`u = U(0,1)*total`, linear scan accumulating cumulative weight) and zeroes the chosen weight so it cannot be picked again — i.e. sequential weighted sampling without replacement, O(n·count). Both paths draw from the shared Mersenne Twister state and deep-copy the selected elements into a `List[...]`. Requesting more than the available count (without `UpTo`) returns `NULL` (unevaluated).
