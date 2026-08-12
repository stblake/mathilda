# AbsoluteTiming

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`AbsoluteTiming[expr] evaluates expr, and returns a list of the absolute number of seconds of elapsed wall-clock time, together with the result obtained.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{seconds, result}`.
- Elapsed real time from a monotonic clock, so a clock adjustment during a long
  evaluation cannot produce a negative interval.
- This, not `Timing`, is the right measurement for anything threaded: the
  multithreaded reductions and elementwise kernels, `Dot` and the LAPACK-backed
  decompositions all run on several cores at once.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## See also

[HoldAll](../../expression-information/HoldAll/), [SequenceHold](../../expression-information/SequenceHold/), [Timing](../../time-and-date/Timing/), [Dot](../../linear-algebra/Dot/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
