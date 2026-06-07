---
source: src/list.c
---
`builtin_differences` computes successive differences, keeping the input head. One pass (`diff_once`) emits `elem[i+s] - elem[i]` for step `s` (reversed for negative `s`), each subtraction built as `Subtract` and reduced via `eval_and_free` so integers, rationals, doubles, symbolics, and matrix rows all combine. `Differences[list, n, s]` applies `diff_once` `n` times with step `s` (`diff_n_step`); `Differences[list, {n1, n2, ...}]` applies per-level first differences recursively into each element (`diff_levels`), e.g. for multidimensional arrays. A list no longer than `|s|` yields the empty list. Non-integer or negative `n`, or step `0`, return `NULL`. This is the additive analog of `Ratios` in the same file.
