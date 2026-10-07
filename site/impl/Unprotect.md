---
source: src/core.c
---
**Algorithm.** `builtin_unprotect` (`src/core.c`) is the mirror of `Protect`: it
calls the shared driver `core_protect_unprotect(res, protecting = false)`, which
walks the argument list (symbols, strings, or a flat `List` of them) and applies
`core_unprotect_one` to each name. `core_unprotect_one` does nothing if the
symbol is not currently `Protected`; otherwise it
clears the `ATTR_PROTECTED` bit and bumps the rule epoch, returning `true` only
when the bit was actually cleared.

The driver returns a `List` of the names (as strings) whose protection state
changed, so `Unprotect` of a symbol that was never protected yields `{}`. Once
the bit is cleared, `Set`/`SetDelayed` and the clearing heads will again accept
the symbol, which is the usual prelude to redefining or extending a built-in's
behaviour.

**Attributes & limits.** `Unprotect` carries `HoldAll | Protected`, so its
arguments arrive as unevaluated names.
