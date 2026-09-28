# Time and Date

## Timing
Evaluates `expr` and returns a list of the time in seconds used, together with the result obtained.
- `Timing[expr]`

**Features**:
- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{timing, result}`.
- Includes only CPU time spent evaluating the expression, **summed over threads**.
  A threaded NDArray path or a BLAS call is therefore over-reported by roughly the
  core count — use `AbsoluteTiming` to measure how long something actually took.

## AbsoluteTiming
Evaluates `expr` and returns a list of the absolute number of seconds of elapsed
wall-clock time, together with the result obtained.
- `AbsoluteTiming[expr]`

**Features**:
- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{seconds, result}`.
- Elapsed real time from a monotonic clock, so a clock adjustment during a long
  evaluation cannot produce a negative interval.
- This, not `Timing`, is the right measurement for anything threaded: the
  multithreaded reductions and elementwise kernels, `Dot` and the LAPACK-backed
  decompositions all run on several cores at once.

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
- The time budget is measured with `ITIMER_PROF` and counts only CPU time (user + kernel) consumed by the Mathilda kernel process. Wall-clock time spent in I/O, sleep, or other processes is not charged.
- A cooperative wall-clock deadline (`CLOCK_MONOTONIC` + budget) is also armed and checked once per evaluator rewrite step, as a portability backstop on hosts where `ITIMER_PROF` is unreliable (notably WSL 1, whose syscall-translation layer under-counts CPU time and delivers `SIGPROF` many seconds late). On real Linux / macOS the `SIGPROF` normally fires first and the cooperative check is a cheap no-op; on broken hosts the cooperative check enforces the deadline between rewrite steps. The only case that escapes both layers is a single long-running C builtin (e.g. `FactorInteger` on a huge composite) on a broken host -- it must wait for the late `SIGPROF`.
- May give different results on different occasions within a single session, for example as a result of different conditions of internal system caches.
- Nested `TimeConstrained` calls compose: each call saves and restores the previous `SIGPROF` handler, `ITIMER_PROF` state, and the cooperative-deadline state, so an inner abort does not disturb an outer time budget.
- The abort is implemented by `siglongjmp`-ing out of the in-flight evaluator. Expression nodes allocated by the aborted computation are not reclaimed; this is the documented behaviour.

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

## DateList
Gives a broken-down date/time as `{year, month, day, hour, minute, second}`.
- `DateList[]` -- current local date and time.
- `DateList[date]` -- date list for the given date specification.

**Supported date specifications**:
- `{y, m, d, h, m, s}` -- a `DateList` specification. Trailing entries may be elided; missing fields default to `{_, 1, 1, 0, 0, 0}`.
- `time` -- a number, taken as an `AbsoluteTime` (seconds since 1900) and inverted.
- `"string"` -- a `DateString` specification, parsed heuristically.
- `{"string", {e1, ...}}` -- a date string parsed with explicit format elements.

**Features**:
- `Protected`. Returns five integers and a real second, so a fractional or current-time result carries sub-second precision (e.g. `{2026, 9, 28, 12, 24, 10.196}`).
- The inverse of `AbsoluteTime`: `DateList[AbsoluteTime[spec]]` and `DateList[spec]` agree, and out-of-range/fractional fields normalize identically (`DateList[{2022,0}] == {2021,12,1,0,0,0.}`, `DateList[{2022,1,0}] == {2021,12,31,0,0,0.}`, `DateList[{2022,3,15.5}] == {2022,3,15,12,0,0.}`).
- Year and month must be integer-valued; a non-integer month gives `DateList::arg` and leaves the input unevaluated. Day, hour, minute, and second may be noninteger.
- **Strings**: tokens are separated by any non-alphanumeric characters. `DateList["28 Sep, 2026"] == {2026,9,28,0,0,0.}`. A purely numeric string with no four-digit year is order-ambiguous: the US `Month/Day/Year` order is used and `DateList::ambig` is issued (`DateList["05/10/1"] == {2001,5,10,0,0,0.}`, two-digit years mapping to `2000+yy`); a leading four-digit token reads as ISO `Year/Month/Day`.
- **Format elements**: `"Year"`, `"YearShort"`, `"Quarter"`, `"Month"`, `"MonthName"`, `"Day"`, `"DayName"`, `"Hour"`, `"Hour12"`, `"AMPM"`, `"Minute"`, `"Second"`, `"Millisecond"`. They are read in the order given; any non-element string between them is treated as a separator. Unfilled fields default to `{current year, 1, 1, 0, 0, 0}`. E.g. `DateList[{"09/28/26",{"Day","Month","YearShort"}}] == {2028,4,9,0,0,0.}` (month `28` reduces to April 2028), `DateList[{"9/28/2026",{"Month","/","Day","/","Year"}}] == {2026,9,28,0,0,0.}`, `DateList[{"2/15",{"Month","Day"}}]` fills the current year.
- Performs no corrections for time zones, daylight saving time, or leap seconds.

## Pause
Pauses for at least `n` seconds, then returns `Null`.
- `Pause[n]`

**Features**:
- `Protected`.
- Sleeps on a wall-clock timer (`nanosleep`) that consumes no CPU, so the elapsed
  time is counted by `AbsoluteTiming` and `SessionTime` but **not** by `Timing`
  or `TimeUsed`. `Timing[Pause[1]]` reports ≈ 0; `AbsoluteTiming[Pause[1]]`
  reports ≈ 1.
- Accurate down to a granularity of at least `$TimeUnit` seconds.
- The sleep resumes across signal interruptions, so the full duration is always
  observed ("at least `n` seconds").
- `n` may be any non-negative number: an integer, real, rational (`Pause[1/4]`),
  or a `NumericQ` symbolic form (`Pause[Pi]`). Zero or negative `n` returns
  immediately. A non-numeric argument leaves `Pause[x]` unevaluated.

## SessionTime
Gives the total number of seconds of wall-clock time elapsed since the beginning
of the current Mathilda session.
- `SessionTime[]`

**Features**:
- `Protected`.
- Measured from a monotonic clock captured at kernel start-up; includes time
  spent in `Pause`.

## TimeUsed
Gives the total number of seconds of CPU time used so far in the current Mathilda
session.
- `TimeUsed[]`

**Features**:
- `Protected`.
- CPU time via `clock()`; does not advance during `Pause` or other idle waits.

## $TimeUnit
Gives the minimum time interval in seconds recorded on the computer system.

**Features**:
- `Protected` (read-only system constant).
- A real equal to the resolution of the `clock()`-based timers
  (`1 / CLOCKS_PER_SEC`), the granularity `Pause` documents itself against.

