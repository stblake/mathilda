---
source: src/list/take_drop.c
---
**Algorithm.** `builtin_take` extracts a subsequence per level. `get_seq_spec_indices`
turns one sequence spec into an explicit index list: `n` (first `n`), `-n` (last `n`),
`{m}` (just element `m`), `{m, n}` (from `m` to `n`), `{m, n, s}` (stepped),
`UpTo[k]`, `All` and `None`, with negative indices resolved from the end. Multiple
specs descend level by level — `apply_take_drop` recurses into each retained element
with the remaining specs — so `Take[matrix, r, c]` slices rows then columns. An
`NDArray` first argument takes the buffer fast path (`ndstruct_take`, a contiguous
leading-axis slice). `Drop` is the same machinery with the complementary index set.

**Data structures.** A small `int64_t` index array per level from the spec; a fresh
node with the original head holding the copied elements.

**Complexity / limits.** O(output size) per level. An out-of-range spec (e.g. taking
more than `Length`) is not silently clamped — `get_seq_spec_indices` returns false,
`builtin_take` returns `NULL`, and the call is left unevaluated (except `UpTo`, which
caps at the available length by design).
