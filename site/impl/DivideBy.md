---
source: src/core.c
---
**Algorithm.** `builtin_divideby` (`src/core.c`) implements `x /= dx` through the
shared `increment_core` helper with `mode = IC_DIV` and `pre = true` — the same
worker behind `AddTo`/`SubtractFrom`/`TimesBy`, which differ only in how the
current value and `dx` are combined. `increment_core` first resolves the lvalue
to the symbol that actually holds the value (`lvalue_symbol_name` accepts a plain
symbol or a `Part[sym, ...]` target); if that symbol has no existing OwnValue it
emits `DivideBy::rvalue` and returns `NULL`, leaving the expression unevaluated.

Otherwise it evaluates the lvalue to the current value, builds and evaluates
`Times[old, dx^-1]` (so list threading and symbolic simplification happen exactly
as for the written-out form), and writes the result back through an evaluated
`Set`. `Set`'s `HoldFirst` preserves a compound lvalue shape such as
`Part[list, i]` for the assignment to update in place. The "pre" flag means the
**new** value is returned.

**Attributes & limits.** `DivideBy` is `HoldFirst | Protected` so the target is
not pre-evaluated; it requires exactly two arguments.
