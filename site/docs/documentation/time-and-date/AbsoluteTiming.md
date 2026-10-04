# AbsoluteTiming

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AbsoluteTiming[expr] evaluates expr, and returns a list of the absolute number of seconds of elapsed wall-clock time, together with the result obtained.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AbsoluteTiming[Sum[i, {i, 1, 1000000}]][[2]]
Out[1]= 500000500000

In[2]:= AbsoluteTiming[Integrate[x^2, x]][[2]]
Out[2]= 1/3 x^3

In[3]:= Length[AbsoluteTiming[1 + 1]]
Out[3]= 2
```

### Applications (3)

The result is the reproducible second element

```mathematica
In[4]:= AbsoluteTiming[Factorial[20]][[2]]
Out[4]= 2432902008176640000
```

Always a two-element {seconds, result} list

```mathematica
In[5]:= Length[AbsoluteTiming[1 + 1]]
Out[5]= 2
```

The elapsed time is never negative

```mathematica
In[6]:= First[AbsoluteTiming[Pause[0]]] >= 0
Out[6]= True
```

## Options & behaviour

The first element of the pair is the elapsed wall-clock time, which varies
between runs; extract the reproducible result with `[[2]]`.

## Implementation notes

**Algorithm.** `builtin_absolute_timing` evaluates its single (held) argument
bracketed by two readings of `dt_wall_seconds()` and returns the two-element
`List` `{elapsed, result}`, where `elapsed` is an `EXPR_REAL`. `dt_wall_seconds`
reads `clock_gettime(CLOCK_MONOTONIC)` and returns `tv_sec + tv_nsec*1e-9`.

**Why wall-clock, not `clock()`.** This is the deliberate difference from
`Timing`. `clock()` reports CPU time summed across threads, so any operation
using `nd_parallel_for`/`nd_parallel_reduce` or the platform BLAS reads roughly
*cores* × its true duration — making `Timing` unusable for the threaded NDArray
paths. `CLOCK_MONOTONIC` (rather than `CLOCK_REALTIME`) also means an NTP step
or manual clock change during a long evaluation cannot produce a negative
interval. If no monotonic clock is available the code falls back to
`clock()/CLOCKS_PER_SEC`, which at least is a duration.

**Data structures & limits.** `struct timespec` plus a two-element `Expr*` List
(`SYM_List`); O(1) overhead around the evaluation it measures.
`ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD` (`src/attr.c`), so the
argument is timed, not pre-evaluated. The elapsed component is non-deterministic
by nature; it is a measurement wrapper, not a numeric kernel, so it exposes no
packed/`Compile[]` surface.

- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{seconds, result}`.
- Elapsed real time from a monotonic clock, so a clock adjustment during a long
  evaluation cannot produce a negative interval.
- This, not `Timing`, is the right measurement for anything threaded: the
  multithreaded reductions and elementwise kernels, `Dot` and the LAPACK-backed
  decompositions all run on several cores at once.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [SequenceHold](../../expression-information/SequenceHold/), [Timing](../../time-and-date/Timing/), [Dot](../../linear-algebra/Dot/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`AbsoluteTiming[expr]` returns `{seconds, result}`, where `seconds` is the
elapsed **wall-clock** time measured from a monotonic clock. Only the second
element is reproducible — the timing varies between runs — so extract it with
`[[2]]`.

This, not `Timing`, is the right measurement for anything threaded. `Timing`
reports CPU time summed over threads, so the multithreaded reductions and
elementwise kernels, `Dot`, and the LAPACK-backed decompositions all read
roughly *core-count* times their true duration there; `AbsoluteTiming` reads the
time you actually waited. Because the clock is monotonic, a clock adjustment
mid-evaluation cannot produce a negative interval, and time spent sleeping in
`Pause` is included.
