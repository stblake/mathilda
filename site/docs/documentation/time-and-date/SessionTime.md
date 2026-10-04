# SessionTime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SessionTime[] gives the total number of seconds of wall-clock time elapsed since the beginning of the current Mathilda session.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SessionTime[] >= 0
Out[1]= True

In[2]:= Head[SessionTime[]]
Out[2]= Real
```

### Applications (2)

Always a Real

```mathematica
In[3]:= Head[SessionTime[]]
Out[3]= Real
```

Grows monotonically from kernel start-up

```mathematica
In[4]:= TrueQ[Positive[SessionTime[]]]
Out[4]= True
```

## Options & behaviour

The elapsed-time value itself grows with every call, so these examples report
its sign and head rather than a fixed number.

## Implementation notes

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

- `Protected`.
- Measured from a monotonic clock captured at kernel start-up; includes time
  spent in `Pause`.

**Attributes:** `Protected`.

## References

**See also:** [Pause](../../time-and-date/Pause/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`SessionTime[]` gives the wall-clock seconds elapsed since the current Mathilda
session began. It is measured from a monotonic clock captured at kernel
start-up — the same clock source as `AbsoluteTiming` — so it increases on every
call and includes time spent sleeping in `Pause`.

The exact value is non-deterministic (it depends on how long the session has
been running), so useful examples report its structure — its head, or its sign —
rather than a fixed number. Contrast `TimeUsed`, which counts only CPU time and
does not advance while the session is idle.
