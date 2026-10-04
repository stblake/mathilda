---
source: src/datetime.c
---
**Algorithm.** `builtin_session_time` returns `dt_wall_seconds() -
g_session_start` as an `EXPR_REAL`: the elapsed wall-clock seconds since the
kernel started. `g_session_start` is a file-static captured once in
`datetime_init` (at start-up) from the same monotonic source
(`clock_gettime(CLOCK_MONOTONIC)`) that `AbsoluteTiming` uses. Because the zero
point and the reading share that monotonic clock, time spent inside `Pause` is
counted here, and no NTP step or manual clock change can make the interval go
backwards.

**Data structures / limits.** One `double` subtraction; `struct timespec` under
`dt_wall_seconds`. Takes no arguments (returns `NULL`, leaving it unevaluated,
for any other arity). `ATTR_PROTECTED`. The value is non-deterministic by
definition (it grows with elapsed time); it is a session clock, not a numeric
kernel, so there is no packed/`Compile[]` surface.
