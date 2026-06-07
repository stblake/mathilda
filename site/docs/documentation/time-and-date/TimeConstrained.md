# TimeConstrained

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TimeConstrained[expr, t]
    evaluates expr, stopping after t seconds.
TimeConstrained[expr, t, failexpr]
    returns failexpr if the time constraint is not met.

TimeConstrained generates an interrupt to abort the evaluation of
expr if the evaluation is not completed within the specified time.
TimeConstrained evaluates failexpr only if the evaluation is
aborted. TimeConstrained returns $Aborted if the evaluation is
aborted and no failexpr is specified.

TimeConstrained[expr, Infinity] imposes no time constraint.
TimeConstrained may give different results on different occasions
within a single session, for example as a result of different
conditions of internal system caches.
TimeConstrained takes account only of CPU time spent inside the
main Mathilda kernel process; it does not include additional
threads or processes.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_time_constrained` implements `TimeConstrained[expr, t]` / `[expr, t, failexpr]` under `HoldAll`: only `t` is evaluated up front (via `tc_parse_time`, which accepts integers, reals, bigints, `Rational`, and `Infinity`/`DirectedInfinity[1]`). `Infinity` evaluates the body unconstrained; a non-positive budget aborts immediately (returning `failexpr` or `$Aborted`).

For a finite positive budget it installs a `SIGPROF` handler (`tc_sigprof_handler`) and arms an `ITIMER_PROF` interval timer for `t` seconds (saving and later restoring any pre-existing handler/timer so nested calls compose). The body is run inside `tc_run_guarded`, which holds the `sigsetjmp` landing pad; on timeout the signal handler `siglongjmp`s back out. As a portability backstop (e.g. WSL where `ITIMER_PROF` is unreliable) it also arms a cooperative wall-clock deadline `tc_deadline` checked by `tc_check_deadline` at each rewrite step using `CLOCK_MONOTONIC`. After a timeout it restores the saved recursion depth (the unwound `evaluate` never decremented it) and returns `failexpr` (evaluated afterward with a fresh budget) or `$Aborted`.

**Limits.** The timer measures CPU (profiling) time, not pure wall-clock; budgets are capped near 2147483 s. On the timeout path the partially-built body tree is leaked — an inherent consequence of unwinding with `siglongjmp`.

- `HoldAll`, `Protected`.
- Returns the value of `expr` if it completes within `t` seconds of CPU time.
- Returns `$Aborted` if the budget is exhausted and no `failexpr` is provided.
- Returns the (then-evaluated) value of `failexpr` if the budget is exhausted and the three-argument form is used. `failexpr` is **not** evaluated when the body completes in time.
- `TimeConstrained[expr, Infinity]` imposes no time constraint.
- Zero, negative, or non-numeric (NaN) time budgets abort immediately without evaluating `expr`; the abort path still produces `$Aborted` (or `failexpr` for the 3-argument form).
- The time budget is measured with `ITIMER_PROF` and counts only CPU time (user + kernel) consumed by the Mathilda kernel process. Wall-clock time spent in I/O, sleep, or other processes is not charged.
- A cooperative wall-clock deadline (`CLOCK_MONOTONIC` + budget) is also armed and checked once per evaluator rewrite step, as a portability backstop on hosts where `ITIMER_PROF` is unreliable (notably WSL 1, whose syscall-translation layer under-counts CPU time and delivers `SIGPROF` many seconds late). On real Linux / macOS the `SIGPROF` normally fires first and the cooperative check is a cheap no-op; on broken hosts the cooperative check enforces the deadline between rewrite steps. The only case that escapes both layers is a single long-running C builtin (e.g. `FactorInteger` on a huge composite) on a broken host -- it must wait for the late `SIGPROF`.
- May give different results on different occasions within a single session, for example as a result of different conditions of internal system caches.
- Nested `TimeConstrained` calls compose: each call saves and restores the previous `SIGPROF` handler, `ITIMER_PROF` state, and the cooperative-deadline state, so an inner abort does not disturb an outer time budget.
- The abort is implemented by `siglongjmp`-ing out of the in-flight evaluator. Expression nodes allocated by the aborted computation are not reclaimed; this is the documented Mathematica behaviour.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TimeConstrained[1 + 1, 5]
Out[1]= 2

In[2]:= TimeConstrained[2^10, 5]
Out[2]= 1024
```

### Notes

`TimeConstrained[expr, t]` returns the result of `expr` if it finishes within `t` seconds; otherwise evaluation is aborted (returning `$Aborted` by default).
