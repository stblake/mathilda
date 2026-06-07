---
source: src/match.c
---
`BlankSequence` is a pattern object head, not an evaluator builtin: `__` parses to `BlankSequence[]` and `__h` to `BlankSequence[h]` (and `x__` to `Pattern[x, BlankSequence[...]]`). It matches a sequence of **one or more** consecutive arguments. The matcher recognises it in `src/match.c` via `is_sequence_blank`, which sets `min_len = 1` and reports the optional head constraint. Sequence matching against an argument list is handled by `match_args_internal`, which backtracks over the possible partitions of the argument run; an optional head `h` requires every element of the matched run to have head `h` (same atomic-type-to-head mapping as `Blank`). Contrast `BlankNullSequence` (`min_len = 0`).
