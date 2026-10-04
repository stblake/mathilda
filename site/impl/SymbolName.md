---
source: src/core.c
---
**Algorithm.** `builtin_symbolname` takes one argument and, when it is a symbol,
returns the symbol's short name as a string with any context prefix stripped. The
implementation finds the last backtick in the interned name with `strrchr(n, '`')`
and returns everything after it (or the whole name if there is no backtick), so
`SymbolName[Global`x]` and `SymbolName[P`Private`x]` both give `"x"`. Anything that
is not a symbol — a number, a compound expression — declines (returns `NULL`) and
the call is left unevaluated, which is the same observable a caller testing the
result's head would get from Wolfram's `SymbolName::sym` path.

**Data structures.** None beyond the argument. The returned value is a freshly
allocated `EXPR_STRING` copied from the tail of the interned name; the name buffer
itself is not modified.

**Complexity / limits.** `O(k)` in the length of the name (one backtick scan). The
argument is evaluated first, so `SymbolName[s]` where `s` has a value reports on
that value, not on `s`; apply it to the bare symbol to name the symbol itself.
Wrong arity declines.
