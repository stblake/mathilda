---
source: src/core.c
---
`builtin_subtractfrom` (`x -= dx`) delegates to the shared `increment_core(lhs, dx, negate=true, pre=true, "SubtractFrom")`. That helper requires `lhs` to be an lvalue with an existing OwnValue (else it emits `SubtractFrom::rvalue` and returns `NULL`), reads the old value via `evaluate`, builds and evaluates `Plus[old, Times[-1, dx]]`, writes the result back through an evaluated `Set[lhs, new]` (preserving lvalue shape such as `Part[...]` via Set's `HoldFirst`), and returns the new value. Has `ATTR_HOLDFIRST` so the target symbol is not evaluated before mutation.
