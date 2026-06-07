---
source: src/core.c
---
**Algorithm.** `builtin_time_constrained` implements `TimeConstrained[expr, t]` / `[expr, t, failexpr]` under `HoldAll`: only `t` is evaluated up front (via `tc_parse_time`, which accepts integers, reals, bigints, `Rational`, and `Infinity`/`DirectedInfinity[1]`). `Infinity` evaluates the body unconstrained; a non-positive budget aborts immediately (returning `failexpr` or `$Aborted`).

For a finite positive budget it installs a `SIGPROF` handler (`tc_sigprof_handler`) and arms an `ITIMER_PROF` interval timer for `t` seconds (saving and later restoring any pre-existing handler/timer so nested calls compose). The body is run inside `tc_run_guarded`, which holds the `sigsetjmp` landing pad; on timeout the signal handler `siglongjmp`s back out. As a portability backstop (e.g. WSL where `ITIMER_PROF` is unreliable) it also arms a cooperative wall-clock deadline `tc_deadline` checked by `tc_check_deadline` at each rewrite step using `CLOCK_MONOTONIC`. After a timeout it restores the saved recursion depth (the unwound `evaluate` never decremented it) and returns `failexpr` (evaluated afterward with a fresh budget) or `$Aborted`.

**Limits.** The timer measures CPU (profiling) time, not pure wall-clock; budgets are capped near 2147483 s. On the timeout path the partially-built body tree is leaked — an inherent consequence of unwinding with `siglongjmp`.
