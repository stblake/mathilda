---
source: src/core.c
---
`builtin_decrement` (`src/core.c`) implements `x--` via the shared `increment_core` helper with a delta of `1`, negate=true, pre=false. `increment_core` requires the target be a symbol with an existing OwnValue (else `Decrement::rvalue`), evaluates the old value, forms and evaluates `Plus[old, Times[-1, 1]]`, and writes it back through an evaluated `Set`. Because pre=false it returns the *old* value (post-decrement). `Decrement` is `ATTR_HOLDFIRST`. The pre-form `--x` is the separate `builtin_predecrement` (pre=true).
