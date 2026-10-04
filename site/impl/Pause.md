---
source: src/datetime.c
---
**Algorithm.** `builtin_pause` blocks for at least `n` seconds of wall-clock
time and returns `Null`. `pause_seconds` coerces the argument to a machine
`double`, accepting integers, reals, bignums, rationals, MPFR reals, and any
`NumericQ` symbolic form (`Pi`, `Sqrt[2]`, …) via `numericalize` — so `Pause[1/4]`
and `Pause[Pi]` behave like Mathematica; a genuinely non-numeric argument leaves
`Pause[x]` unevaluated. The double is split into whole seconds and nanoseconds
and handed to `nanosleep`, re-armed with the reported remainder whenever it
returns early on a signal (`EINTR`), so the full duration is always observed.
`n <= 0` (or non-finite) waits not at all, matching `Pause[0]`.

**CPU vs wall clock.** `nanosleep` consumes no CPU, so the elapsed time is
counted by the wall clocks (`AbsoluteTiming`, `SessionTime`) but is invisible to
the CPU clocks (`Timing`, `TimeUsed`) — the two clock sources distinguish the
cases for free, with no special-casing. Accuracy is bounded by the OS timer
granularity (`$TimeUnit`-scale).

**Data structures / limits.** `struct timespec` and a scalar; returns
`expr_new_symbol(SYM_Null)`. `ATTR_PROTECTED`. A blocking wall-clock wait, not a
numeric kernel — no packed/`Compile[]` surface.
