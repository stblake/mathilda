---
source: src/core.c
---
**Algorithm.** `$TimeUnit` is not a function but a read-only system constant,
registered in `system_constants_init` (`src/core.c`) by
`register_system_constant("$TimeUnit", expr_new_real(1.0 / (double)CLOCKS_PER_SEC))`.
Its value is the reciprocal of the C library's `CLOCKS_PER_SEC`, i.e. the
resolution of the `clock()`-based CPU timers that `Timing` and `TimeUsed` read,
and the granularity `Pause` documents itself against. On a standard POSIX host
`CLOCKS_PER_SEC` is `1000000`, so `$TimeUnit` is `1.*10^-6`.

**Data structures.** A single `EXPR_REAL` bound as the symbol's `OwnValue`;
nothing is recomputed per reference. `register_system_constant` also sets
`ATTR_PROTECTED`, so the binding cannot be reassigned.

**Complexity / limits.** O(1) lookup. The value is a machine `double` fixed at
start-up from a compile-time constant; it is a scalar system parameter, not a
numeric kernel, so there is no packed/`Compile[]` surface to provide.
