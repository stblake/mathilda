---
source: src/core.c
---
`builtin_addto` (`src/core.c`) implements `x += dx` via the shared `increment_core` helper (negate=false, pre=true). `increment_core` requires the lvalue to be a symbol with an existing OwnValue (else `AddTo::rvalue`), evaluates the current value, builds and evaluates `Plus[old, dx]`, then writes the new value back through an evaluated `Set` call (Set's `HoldFirst` preserves complex lvalue shapes like `Part[list, i]`). The "pre" flag means it returns the new value. `AddTo` itself is `ATTR_HOLDFIRST` so the target is not pre-evaluated.
