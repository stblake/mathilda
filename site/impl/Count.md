---
source: src/patterns.c
---
**Algorithm.** `builtin_count` (`src/patterns.c`) tallies the subexpressions matching a pattern. Level-spec (default `{1,1}` — immediate elements only), the `Heads -> True|False` option, and the argument shapes mirror `Cases`. The worker `do_count_at_level` recurses depth-first into the head (when `heads`) and every argument, and at each in-range node calls `match(e, pattern, env)` from `src/match.c`, incrementing a `size_t` counter on success; level membership for negative specs is resolved via `get_expr_depth_patterns`. Unlike `Cases`/`Position` it stores nothing and has no result limit — it just returns the final count as an `EXPR_INTEGER`. `Count[pat]` with one argument returns the operator form `Function[Count[#1, pat]]`.

**Data structures.** A single `size_t` accumulator; one `MatchEnv` per node tested. No results buffer.
