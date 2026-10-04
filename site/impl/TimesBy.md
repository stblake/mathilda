---
source: src/core.c
---
**Algorithm.** `builtin_timesby` (`src/core.c`) implements `x *= dx` through the
shared `increment_core` helper with `mode = IC_MUL` and `pre = true`. The mode is
kept rather than a `negate` flag precisely because `TimesBy`/`DivideBy` differ
from `AddTo`/`SubtractFrom` in the combining **head**, not the sign: every other
step of the six increment operators is identical. `increment_core` resolves the
lvalue to its backing symbol (`lvalue_symbol_name` handles a plain symbol or a
`Part[sym, ...]` target) and, if that symbol has no OwnValue, emits
`TimesBy::rvalue` and returns `NULL` (the call stays unevaluated).

Otherwise it evaluates the current value, builds and evaluates `Times[old, dx]`,
then writes the result back via an evaluated `Set` (whose `HoldFirst` preserves a
`Part` lvalue so the mutation lands in place). The "pre" flag returns the **new**
value.

**Attributes & limits.** `TimesBy` is `HoldFirst | Protected` so the target is
not pre-evaluated; it requires exactly two arguments.
