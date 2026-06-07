---
source: src/match.c
---
`Repeated[p]` (`p..`) is a pattern object handled entirely inside the matcher, not by a builtin. `is_repeated` in `src/match.c` recognises the `Repeated` head, sets the matched-length range to `[1, ∞)` by default, and parses an optional count spec: `Repeated[p, max]` gives `[1, max]`, `Repeated[p, {n}]` gives exactly `n`, `Repeated[p, {min, max}]` gives `[min, max]` (with `Infinity` allowed as an open upper bound). The argument-sequence matcher then matches a run of consecutive arguments each satisfying `p`, using the standard backtracking that explores valid run lengths within the range.
