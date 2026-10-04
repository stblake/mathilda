# Pause

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pause[n] pauses for at least n seconds.`**

<details>
<summary>Notes</summary>

Pause is accurate only down to a granularity of at least $TimeUnit seconds. The time elapsed during the execution of Pause is counted in SessionTime and AbsoluteTiming, but not in TimeUsed or Timing.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Pause[0] === Null
Out[1]= True

In[2]:= Pause[1/100] === Null
Out[2]= True
```

### Applications (2)

Pause returns Null; zero waits not at all

```mathematica
In[3]:= Pause[0] === Null
Out[3]= True
```

A rational (or any NumericQ) duration is accepted

```mathematica
In[4]:= Pause[1/100] === Null
Out[4]= True
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [AbsoluteTiming](../../time-and-date/AbsoluteTiming/), [SessionTime](../../time-and-date/SessionTime/), [Timing](../../time-and-date/Timing/), [TimeUsed](../../time-and-date/TimeUsed/), [$TimeUnit](../../time-and-date/$TimeUnit/), [NumericQ](../../expression-information/NumericQ/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`Pause[n]` blocks for at least `n` seconds of wall-clock time, then returns
`Null`. The argument may be any non-negative number — integer, real, rational
(`Pause[1/4]`), or a `NumericQ` symbolic form (`Pause[Pi]`); zero or negative
returns immediately, and a non-numeric argument leaves `Pause[x]` unevaluated.

The wait is a `nanosleep`, which consumes no CPU. So the elapsed time is counted
by `AbsoluteTiming` and `SessionTime` but **not** by `Timing` or `TimeUsed`:
`Timing[Pause[1]]` reports about `0`, while `AbsoluteTiming[Pause[1]]` reports
about `1`. The sleep resumes across signal interruptions, so the full duration is
always observed. Accuracy is bounded by the timer granularity `$TimeUnit`.
