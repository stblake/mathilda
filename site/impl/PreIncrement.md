---
source: src/core.c
---
`builtin_preincrement` (`++x`) calls the shared `increment_core(lhs, 1, negate=false, pre=true, "PreIncrement")`. The `pre=true` flag makes it return the *new* value: it reads the old OwnValue via `evaluate`, computes and writes back `Plus[old, 1]` through an evaluated `Set` (which preserves lvalue shape via Set's `HoldFirst`), then returns the new value and frees the old. Requires the target to already have an OwnValue, else emits `PreIncrement::rvalue` and returns `NULL`. Carries `ATTR_HOLDFIRST | ATTR_PROTECTED`.
