---
source: src/datetime.c
---
**Algorithm.** `builtin_absolute_timing` evaluates its single (held) argument
bracketed by two readings of `dt_wall_seconds()` and returns the two-element
`List` `{elapsed, result}`, where `elapsed` is an `EXPR_REAL`. `dt_wall_seconds`
reads `clock_gettime(CLOCK_MONOTONIC)` and returns `tv_sec + tv_nsec*1e-9`.

**Why wall-clock, not `clock()`.** This is the deliberate difference from
`Timing`. `clock()` reports CPU time summed across threads, so any operation
using `nd_parallel_for`/`nd_parallel_reduce` or the platform BLAS reads roughly
*cores* × its true duration — making `Timing` unusable for the threaded NDArray
paths. `CLOCK_MONOTONIC` (rather than `CLOCK_REALTIME`) also means an NTP step
or manual clock change during a long evaluation cannot produce a negative
interval. If no monotonic clock is available the code falls back to
`clock()/CLOCKS_PER_SEC`, which at least is a duration.

**Data structures & limits.** `struct timespec` plus a two-element `Expr*` List
(`SYM_List`); O(1) overhead around the evaluation it measures.
`ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD` (`src/attr.c`), so the
argument is timed, not pre-evaluated. The elapsed component is non-deterministic
by nature; it is a measurement wrapper, not a numeric kernel, so it exposes no
packed/`Compile[]` surface.
