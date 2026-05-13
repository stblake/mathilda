# Time and Date

## Timing
Evaluates `expr` and returns a list of the time in seconds used, together with the result obtained.
- `Timing[expr]`

**Features**:
- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{timing, result}`.
- Includes only CPU time spent evaluating the expression.

## RepeatedTiming
Evaluates `expr` repeatedly and returns a list of the average time in seconds used, together with the result obtained.
- `RepeatedTiming[expr]`
- `RepeatedTiming[expr, t]`

**Features**:
- `HoldFirst`, `Protected`, `SequenceHold`.
- Returns `{average_timing, result}`.
- Does repeated evaluation for at least `t` seconds. Default is 1 second.
- Gives a trimmed mean of the timings obtained, discarding lower and upper quartiles.
- Always evaluates `expr` at least four times.

## TimeConstrained
Evaluates `expr`, generating an interrupt to abort the evaluation if it has not completed within the time budget.
- `TimeConstrained[expr, t]`
- `TimeConstrained[expr, t, failexpr]`

**Features**:
- `HoldAll`, `Protected`.
- Returns the value of `expr` if it completes within `t` seconds of CPU time.
- Returns `$Aborted` if the budget is exhausted and no `failexpr` is provided.
- Returns the (then-evaluated) value of `failexpr` if the budget is exhausted and the three-argument form is used. `failexpr` is **not** evaluated when the body completes in time.
- `TimeConstrained[expr, Infinity]` imposes no time constraint.
- Zero, negative, or non-numeric (NaN) time budgets abort immediately without evaluating `expr`; the abort path still produces `$Aborted` (or `failexpr` for the 3-argument form).
- The time budget is measured with `ITIMER_PROF` and counts only CPU time (user + kernel) consumed by the PicoCAS kernel process. Wall-clock time spent in I/O, sleep, or other processes is not charged.
- May give different results on different occasions within a single session, for example as a result of different conditions of internal system caches.
- Nested `TimeConstrained` calls compose: each call saves and restores the previous `SIGPROF` handler and `ITIMER_PROF` state so an inner abort does not disturb an outer time budget.
- The abort is implemented by `siglongjmp`-ing out of the in-flight evaluator. Expression nodes allocated by the aborted computation are not reclaimed; this is the documented Mathematica behaviour.

## AbsoluteTime
Gives the total number of seconds since the beginning of January 1, 1900.
- `AbsoluteTime[]` -- current wall-clock time, in the local time zone.
- `AbsoluteTime[date]` -- absolute time corresponding to the given date specification.

**Supported date specifications**:
- `{y, m, d, h, m, s}` -- `DateList`-style specification. Trailing entries may be elided; missing fields default to `{_, 1, 1, 0, 0, 0}`.
- `time` -- a number (`AbsoluteTime` specification); returned unchanged.

**Features**:
- `Protected`.
- Year and month must be integer-valued; day, hour, minute, and second may be noninteger.
- Out-of-range date components are converted to standard normalized form, e.g. `AbsoluteTime[{2022, 2, 31}] == AbsoluteTime[{2022, 3, 3}] == 3855254400`.
- Performs no corrections for time zones, daylight saving time, or leap seconds.
- Returns an integer when every component is integer-valued and the total is exact; otherwise returns a real.

