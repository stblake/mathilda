---
source: src/match.c
---
`OneIdentity` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_ONEIDENTITY` bitflag. The bit is consumed by the pattern matcher: in `match.c`, when matching a single argument against a sequence/optional pattern under a head carrying `ATTR_ONEIDENTITY`, `f[x]` is treated as equivalent to `x` (and an optional/defaulted argument lets `x` match `f[x, default]`). This is what lets `a + b_.` match `a` (with `b` bound to the identity 0) for `Plus`, or `a x_.` match `a` for `Times`. The matcher checks `def->attributes & ATTR_ONEIDENTITY` at the relevant decision points. `purefunc.c` maps `SYM_OneIdentity` → `ATTR_ONEIDENTITY` for pure-function attribute specs.
