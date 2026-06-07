---
source: src/match.c
---
`RepeatedNull[p]` (`p...`) is the zero-or-more variant of `Repeated`, handled in the matcher. `is_repeated` in `src/match.c` recognises the head and sets the default length range to `[0, ∞)` (the only difference from `Repeated`'s `[1, ∞)`); the same optional count specs apply (`max`, `{n}`, `{min, max}` with `Infinity` permitted). The argument matcher backtracks over run lengths down to and including the empty run.
