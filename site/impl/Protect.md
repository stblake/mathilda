---
source: src/core.c
---
**Algorithm.** `builtin_protect` (`src/core.c`) calls the shared driver
`core_protect_unprotect(res, protecting = true)`. That driver walks the argument
list — each spec being a symbol, a string, or a flat `List` of them — and applies
`core_protect_one` to every name. `core_protect_one` leaves a `Locked` symbol
untouched and does nothing if the symbol is already `Protected`; otherwise it
sets the `ATTR_PROTECTED` bit and bumps the rule epoch (invalidating the
evaluation cache), returning `true` only when the bit was *newly* set.

The driver collects the names whose state actually changed into a growable
`Expr**` buffer and returns them as a `List` of strings — matching the Wolfram
convention where `Protect` reports exactly what it altered, so re-protecting an
already-protected symbol yields `{}`. The `Protected` attribute is what the
evaluator and `Set` consult to refuse redefinition of a symbol.

**Attributes & limits.** `Protect` carries `HoldAll | Protected`, so its
arguments reach the handler as unevaluated symbol names.
