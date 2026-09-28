# Pause

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pause[n] pauses for at least n seconds.`**

<details>
<summary>Notes</summary>

Pause is accurate only down to a granularity of at least $TimeUnit seconds. The time elapsed during the execution of Pause is counted in SessionTime and AbsoluteTiming, but not in TimeUsed or Timing.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)
