---
source: src/core.c
---
`builtin_predecrement` (`--x`) calls the shared `increment_core(lhs, 1, negate=true, pre=true, "PreDecrement")`. It writes back `Plus[old, Times[-1, 1]]` via an evaluated `Set` and, because `pre=true`, returns the new (decremented) value. Requires an existing OwnValue on the target, else emits `PreDecrement::rvalue` and returns `NULL`. Carries `ATTR_HOLDFIRST | ATTR_PROTECTED`.
