---
source: src/datetime.c
---
**Algorithm.** `builtin_time_used` returns `(double)clock() / CLOCKS_PER_SEC` as
an `EXPR_REAL`: the total CPU seconds consumed so far in the current session. It
reads the same processor clock that `Timing` brackets, so — like `Timing` — it
does **not** advance during `Pause` or other idle waits (where no CPU is burned),
and it counts CPU time summed over threads.

**Data structures / limits.** A single `clock()` call and a divide; takes no
arguments (returns `NULL` for any other arity). `ATTR_PROTECTED`. Its resolution
is `$TimeUnit` (`1/CLOCKS_PER_SEC`). The value is non-deterministic (it grows
with work done); it is a session CPU clock, not a numeric kernel, so there is no
packed/`Compile[]` surface.
