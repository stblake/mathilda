---
source: src/match.c
---
`Optional` is a pattern wrapper, not a function. The matcher in `match.c` peels `Optional[p]` / `Optional[p, default]` off a pattern (in the same wrapper-stripping loop that handles `Pattern`, `Shortest`, `Longest`), sets `is_optional`, and records the default. When the optional argument is absent at that position, the bound symbol is filled with the explicit `default` when given (`opt_container->args[1]`), otherwise with `get_default_value(pat_head, pos, total)` — which supplies the head's identity (0 for `Plus`, 1 for `Times`, etc., the head's `Default[]` value). When the argument *is* present it matches `p` normally. This is the mechanism behind the `x_.` / `x_:def` surface syntax. `Optional` is in the set of pattern heads `eval.c` leaves unevaluated.
