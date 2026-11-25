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

```
In[1]:= AbsoluteTiming[Sum[i, {i, 1, 1000000}]][[2]]
Out[1]= 500000500000
In[2]:= AbsoluteTiming[Integrate[x^2, x]][[2]]
Out[2]= 1/3 x^3
In[3]:= Length[AbsoluteTiming[1 + 1]]
Out[3]= 2
```

The first element of the pair is the elapsed wall-clock time, which varies
between runs; extract the reproducible result with `[[2]]`.


## EchoTiming
Evaluates `expr`, prints the elapsed wall-clock seconds as `>> t`, and returns
the result.
- `EchoTiming[expr]`
- `EchoTiming[expr, label]` prints `>> label t`.

**Features**:
- `HoldFirst`, `Protected`, `SequenceHold`.
- Uses the same monotonic clock as `AbsoluteTiming`. The line is printed
  through `Print` (see `Echo`).

```
In[1]:= EchoTiming[Pause[0.1]; 7, "slept:"]
>> slept: 0.100147
Out[1]= 7
```

## RepeatedTiming
- An exact time budget (`1/10`, `Pi/100`) is read through `N` (audit RepeatedTiming#1).
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
- A time budget that is not a number or `Infinity` (`x`, `-Infinity`, `1 + I`) prints `TimeConstrained::timc` and the call stays unevaluated (audit TimeConstrained#2).
Evaluates `expr`, generating an interrupt to abort the evaluation if it has not completed within the time budget.
- `TimeConstrained[expr, t]`
- `TimeConstrained[expr, t, failexpr]`

**Features**:
- `HoldAll`, `Protected`.
- Returns the value of `expr` if it completes within `t` seconds of elapsed (wall-clock) time, as in the reference.
- Returns `$Aborted` if the budget is exhausted and no `failexpr` is provided.
- Returns the (then-evaluated) value of `failexpr` if the budget is exhausted and the three-argument form is used. `failexpr` is **not** evaluated when the body completes in time.
- `TimeConstrained[expr, Infinity]` imposes no time constraint.
- **The timeout is the ordinary `Abort[]`** (v0.6623). It does not jump. When the deadline passes, the
  TimeConstrained's frame on the constraint stack (`src/memacct.h`, shared with `MemoryConstrained`) is
  *tripped*. The evaluator's poll then turns the trip into the `Abort[]` sentinel, every frame between it and
  the TimeConstrained leaves through its own cleanup, and the TimeConstrained consumes the abort and returns
  `failexpr`. This matches the reference:
  - nested budgets: the outermost expired one answers;
  - an `AbortProtect` *inside* the body runs to completion first, while one enclosing the TimeConstrained
    does not stop it;
  - `CheckAbort` runs its `failexpr`, then (under the default `PropagateAborts -> Automatic`) the timeout
    goes on. With `PropagateAborts -> False` it is absorbed, and that frame never fires again;
  - `Catch` and `Check` do not stop it;
  - a timeout during a lazily loaded package's `Get` puts `$Context` / `$ContextPath` back.
- **What changed.** Before, the cooperative check `siglongjmp`'d out of the evaluation, skipping every C
  frame in between. That leaked whatever those frames held, and it bypassed `CheckAbort` and `AbortProtect`.
  Now:
  - `TimeConstrained[While[True], 0.3]` leaks nothing under valgrind.
  - The jump is only the backstop. It is taken when a computation that reaches no evaluator step (a pure C
    loop), or that ignores the `Abort[]` it gets back, still runs `TC_ASYNC_GRACE_US` (50 ms CPU) after the
    deadline, and never while an `AbortProtect` inside the body runs. That path alone can still leak.
  - A kernel that declines on the deadline must not fall back to a slower path: the MPFR `Det` /
    `LUDecomposition` / `LinearSolve` callers now stop there.
  - Messages from a computation that is being aborted are dropped by the funnel, as they already were for
    Ctrl-C, so a C loop that took the `Abort[]` as data cannot make an enclosing `Check` fire.
- **A FLINT exception aborts like a timeout** (v0.6090).
  - **What happens.** FLINT reports an internal error (for example `qqbar_pow_fmpq: division by zero`) through
    `flint_throw`, whose default is `abort()`. Mathilda installs a handler (`tc_flint_throw`, `src/core.c`) that
    unwinds to the innermost recovery point instead and reports `FLINT::exc`, honouring the `Quiet` in force where
    FLINT threw. `Check` sees the message.
  - **Inside `TimeConstrained`:** the body is abandoned and the call returns `failexpr` or `$Aborted`.
  - **Outside it:** every top-level evaluation runs under a fault-only recovery point (`tc_evaluate_fault_guarded`;
    no signal-mask save, no measurable cost), so the input evaluates to `$Aborted` and the kernel carries on.
  - **In heavy methods:** `Integrate`'s Trager and ParallelMixedTower method stages evaluate under their own recovery
    point, so a fault there is a decline and the cascade continues.
  - **Where it cannot apply:** off the main thread there is no recovery point, and FLINT's default applies.
- Zero, negative, or non-numeric (NaN) time budgets abort immediately without evaluating `expr`; the abort path still produces `$Aborted` (or `failexpr` for the 3-argument form).
- **The budget is elapsed time** (v0.463): the deadline is one `CLOCK_MONOTONIC` instant, `now + t`. Two layers enforce it. (1) A cooperative check once per evaluator rewrite step, at heavy-kernel loop heads (`tc_deadline_passed`), and in `Pause`, which sleeps no further than the deadline. (2) On POSIX, an `ITIMER_PROF`/`SIGPROF` timer that preempts a long C builtin which never returns to the evaluator. That timer counts process CPU time summed over threads, so it serves only as a *wake-up*: its handler compares the monotonic clock with the deadline and, when woken early (N busy threads accrue CPU N times faster than wall time), re-arms for the remainder. Threaded `FactorInteger` regions (QS, ECM, p±1, BlakeRationalBaseDescent) poll the deadline from their workers — GMP-ECM through its `stop_asap` hook — and the main thread takes the timeout as the region ends.
- Before v0.463 the CPU timer *was* the budget. Threaded ECM/QS runs therefore aborted 10–16× early in elapsed time: `TimeConstrained[FactorInteger`ECM[NextPrime[3 10^29] NextPrime[7 10^29]], 3]` gave `$Aborted` after 0.3 s. Time spent idle was not charged either: `TimeConstrained[Pause[5], 1]` slept the full 5 s and returned `Null`. Both now stop at the budget.
- **The timeout is taken at a safe point** (v0.592). At the deadline the `SIGPROF` handler only records the timeout and gives the computation 50 ms of CPU to reach the evaluator's next rewrite-step check or a guarded GMP allocation, and the jump happens there. When that grace runs out (a long C loop), the handler reads the interrupted PC and frame-pointer chain and unwinds only once no frame lies in the allocator, the C library or pthread, retrying every millisecond until then (`src/tc_unwind_safe.c`; macOS and Linux on x86-64 and arm64). An asynchronous jump out of `malloc` leaves the allocator's lock held, and the next allocation then traps the process.
- On hosts whose `ITIMER_PROF` is unreliable (notably WSL 1, which under-counts CPU time and delivers `SIGPROF` late), and always on native Windows, which has no `SIGPROF`, the cooperative layer alone enforces the deadline between rewrite steps. A single long-running C builtin that neither returns nor polls is the one case it cannot interrupt.
- May give different results on different occasions within a single session, for example as a result of different conditions of internal system caches.
- Nested `TimeConstrained` calls compose: each call saves and restores the previous `SIGPROF` handler, the timer, and the deadline state, so an inner abort does not disturb an outer time budget.
- A timeout restores the **context state**: the `Begin`/`BeginPackage` stack and the values of `$Context` and `$ContextPath` go back to what they were when the body started (v0.6501).
  - Before this, a timeout that landed while a package was being read (Integrate's lazily loaded Roach engine, for example) left the session in `Roach`Private`` with `$ContextPath = {"System`"}`. Every later bare symbol then resolved to a new private symbol, so a script's own variables read back unevaluated.
  - A body that *completes* keeps its context changes, so `TimeConstrained[Begin["X`"], 1]` leaves `X`` current.
- A timeout also restores every C nesting counter a subsystem has registered with `tc_register_unwind_counter` (v0.451) — the unwind skips their decrements, and a counter left raised changes behaviour for the rest of the session (one aborted nested `DSolve` used to leave `DSolve`'s recursion depth at 2, disabling its top-level-only methods and its failure memo).
- **Nesting CLAMPS, it does not refund** (v0.253). An inner budget can only ever *shorten* the enclosing one: the inner deadline is `min(now + inner, outer deadline)`, and since the outer's deadline is an absolute instant, the time the inner call spent is charged to it automatically (the outer timer is re-armed for `outer deadline − now` on exit). So `TimeConstrained[expr, t]` is an upper bound on the whole call, however many nested budgets `expr` enters. A shorter inner budget still wins, and a computation that fits inside its budget is unaffected. Before this, both layers lifted the caller's deadline merely by entering an inner scope — the `ITIMER_PROF` layer reinstalled the outer timer at the value it held when the inner call *started*, and the cooperative layer set an inner absolute deadline that could be later than the outer's — so `TimeConstrained[TimeConstrained[<long loop>, 30], 3]` ran the loop to completion in 24.5 s. Since `Integrate`, `Simplify`, `NIntegrate` and `DSolve` all bound their own sub-steps this way, no budget was an upper bound in practice.
- The abort is implemented by `siglongjmp`-ing out of the in-flight evaluator. Expression nodes allocated by the aborted computation are not reclaimed; this is the documented behaviour.

```mathematica
In[1]:= TimeConstrained[2 + 2, 1, Print["never printed"]; -1]
Out[1]= 4

In[2]:= TimeConstrained[While[True], 0.2, "timed out"]
Out[2]= "timed out"

In[3]:= TimeConstrained[Print["not run"], 0]
Out[3]= $Aborted

In[4]:= TimeConstrained[{TimeConstrained[While[True], 0.1, "inner"], 7}, 5]
Out[4]= {"inner", 7}

In[5]:= AbsoluteTiming[TimeConstrained[Pause[3], 0.3]]
Out[5]= {0.300412, $Aborted}
```

## Internal`TimeConstrainedRemaining
The seconds left before the innermost active `TimeConstrained` deadline.
- `Internal`TimeConstrainedRemaining[]`

**Features**:
- `Protected`. Returns a Real: `0.` once the deadline has passed, and `Infinity` outside any `TimeConstrained`.
- A nested budget is already clamped to the outer one, so the value is the time the
  *enclosing* computation has left: `TimeConstrained[TimeConstrained[Internal`TimeConstrainedRemaining[], 50], 2]`
  is about `2.`.
- For package code that holds a verified answer and has an optional step to spend time on.
  When a nested `TimeConstrained` expires because the outer budget is exhausted, the
  *outer* call aborts, and the answer is lost with it. Bronstein's real-form presentation
  therefore takes at most a third of what is left (v0.6526).

## AbsoluteTime
- An all-integer `List` spec whose total passes 2^53 is summed exactly in GMP (audit AbsoluteTime#4).
Gives the total number of seconds since the beginning of January 1, 1900.
- `AbsoluteTime[]` -- current wall-clock time, in the local time zone, as a machine real with sub-second resolution (so `AbsoluteTime[] - t0` times a computation).
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

Accepts date strings and `{"string", {elements}}` specs (the DateList parsers), and
exact rational times and fields.

## DateList
- The fraction of the seconds field is carried apart from the absolute total, so `DateList[{2022, 5, 6, 7, 8, 9.123456789}]` keeps `9.123456789` (audit DateList#8).
Gives a broken-down date/time as `{year, month, day, hour, minute, second}`.
- `DateList[]` -- current local date and time.
- `DateList[date]` -- date list for the given date specification.
- The proleptic Gregorian calendar holds for every year: the day-number
  conversions use floor division, so `DateList[{-5000, 3, 15}]` round-trips
  (audit DateList#1).

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

```
In[1]:= DateList[{2026, 10, 4}]
Out[1]= {2026, 10, 4, 0, 0, 0.0}
In[2]:= DateList[0]
Out[2]= {1900, 1, 1, 0, 0, 0.0}
In[3]:= DateList[{2022, 2, 31}]
Out[3]= {2022, 3, 3, 0, 0, 0.0}
In[4]:= DateList["28 Sep, 2026"]
Out[4]= {2026, 9, 28, 0, 0, 0.0}
```

Free-form strings: numbers joined by `:` are the clock time, `AM`/`PM` applies a
12-hour clock (also to a lone trailing hour, `"Sep 28 2026 3 PM"`), an ISO `T`
separator and ordinal suffixes are skipped, a time alone is today at that time, and an
all-numeric date whose first field exceeds 12 is day-first (no ambiguity warning).
Rational AbsoluteTimes and fields are accepted.

## UnixTime

- An exact Integer AbsoluteTime is shifted exactly: `UnixTime[10^30]` is `10^30 - 2208988800` (it went through a double; audit UnixTime#4).
Gives the total number of seconds since the Unix epoch, 1970-01-01 00:00:00 GMT.
- `UnixTime[]` -- current time, the true POSIX epoch second.
- `UnixTime[date]` -- Unix time corresponding to the given date specification.
- A whole number of seconds outside the int64 range is an exact big integer
  (`UnixTime[2.^63 + 2208988800.]` is `2^63`; audit UnixTime#1).

**Supported date specifications**:
- `{y, m, d, h, m, s}` -- a `DateList` specification. Trailing entries may be elided; missing fields default to `{_, 1, 1, 0, 0, 0}`.
- `time` -- a number, taken as an `AbsoluteTime` (seconds since 1900).
- `"string"` -- a `DateString` specification, parsed heuristically.
- `{"string", {e1, ...}}` -- a date string parsed with explicit format elements.

**Features**:
- `Protected`.
- `AbsoluteTime` shifted by the fixed 1900→1970 epoch offset `days_since_1900[1970,1,1] * 86400 == 2208988800`: `UnixTime[spec] == AbsoluteTime[spec] - 2208988800`, so `UnixTime[{1970,1,1,0,0,0}] == 0` and `UnixTime[{2022,1,1,0,0,0}] == 1640995200`.
- **Always returns an integer**, the nearest whole second — a fractional-second spec is rounded rather than kept as a real (a deliberate difference from `AbsoluteTime`).
- `UnixTime[]` reports the true POSIX second (GMT), matching the operating system clock (`date +%s`); the date-list and string forms are timezone-free calendar arithmetic and interpret a spec exactly as `DateList` does.
- Year and month must be integer-valued; a non-integer month gives `UnixTime::arg` and leaves the input unevaluated. Day, hour, minute, and second may be noninteger.
- Out-of-range date components are reduced to standard normalized form before rounding, e.g. `UnixTime[{2022, 2, 31}] == UnixTime[{2022, 3, 3}] == 1646265600`.
- Performs no corrections for time zones, daylight saving time, or leap seconds.

```
In[1]:= UnixTime[{1970, 1, 1, 0, 0, 0}]
Out[1]= 0
In[2]:= UnixTime[{2022, 1, 1, 0, 0, 0}]
Out[2]= 1640995200
In[3]:= UnixTime[{2022, 2, 31}]
Out[3]= 1646265600
In[4]:= AbsoluteTime[{2000, 1, 1}] - UnixTime[{2000, 1, 1}]
Out[4]= 2208988800
```

Rational times and fields are accepted; the nearest second rounds a half to even.

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
- The sleep resumes across signal interruptions, so the full duration is
  observed ("at least `n` seconds"), except an interactive interrupt (Ctrl-C),
  which ends it within 50 ms (not inside `AbortProtect`).
- `n` may be any non-negative number: an integer, real, rational (`Pause[1/4]`),
  or a `NumericQ` symbolic form (`Pause[Pi]`). Zero or negative `n` returns
  immediately. A non-numeric argument leaves `Pause[x]` unevaluated.

```
In[1]:= Pause[0] === Null
Out[1]= True
In[2]:= Pause[1/100] === Null
Out[2]= True
```

## SessionTime
Gives the total number of seconds of wall-clock time elapsed since the beginning
of the current Mathilda session.
- `SessionTime[]`

**Features**:
- `Protected`.
- Measured from a monotonic clock captured at kernel start-up; includes time
  spent in `Pause`.

```
In[1]:= SessionTime[] >= 0
Out[1]= True
In[2]:= Head[SessionTime[]]
Out[2]= Real
```

The elapsed-time value itself grows with every call, so these examples report
its sign and head rather than a fixed number.

## TimeUsed
Gives the total number of seconds of CPU time used so far in the current Mathilda
session.
- `TimeUsed[]`

**Features**:
- `Protected`.
- CPU time via `clock()`; does not advance during `Pause` or other idle waits.

```
In[1]:= TimeUsed[] >= 0
Out[1]= True
In[2]:= Head[TimeUsed[]]
Out[2]= Real
```

The CPU-time value grows with the work done, so these examples report its sign
and head rather than a fixed number.

## $TimeUnit
Gives the minimum time interval in seconds recorded on the computer system.

**Features**:
- `Protected` (read-only system constant).
- A real equal to the resolution of the `clock()`-based timers
  (`1 / CLOCKS_PER_SEC`), the granularity `Pause` documents itself against.

```
In[1]:= $TimeUnit
Out[1]= 1e-06
In[2]:= Head[$TimeUnit]
Out[2]= Real
```

The value is machine-dependent, but on a standard POSIX host `CLOCKS_PER_SEC` is
`1000000`, so `$TimeUnit` is `1.*10^-6`.

