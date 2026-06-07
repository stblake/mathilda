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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
