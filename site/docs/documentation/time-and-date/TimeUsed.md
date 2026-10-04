# TimeUsed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TimeUsed[] gives the total number of seconds of CPU time used so far in the current Mathilda session.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= TimeUsed[] >= 0
Out[1]= True

In[2]:= Head[TimeUsed[]]
Out[2]= Real
```

### Applications (2)

Always a Real

```mathematica
In[3]:= Head[TimeUsed[]]
Out[3]= Real
```

CPU seconds consumed so far this session

```mathematica
In[4]:= TimeUsed[] >= 0
Out[4]= True
```

## Options & behaviour

The CPU-time value grows with the work done, so these examples report its sign
and head rather than a fixed number.

## Implementation notes

**Algorithm.** `builtin_time_used` returns `(double)clock() / CLOCKS_PER_SEC` as
an `EXPR_REAL`: the total CPU seconds consumed so far in the current session. It
reads the same processor clock that `Timing` brackets, so — like `Timing` — it
does **not** advance during `Pause` or other idle waits (where no CPU is burned),
and it counts CPU time summed over threads.

**Data structures / limits.** A single `clock()` call and a divide; takes no
arguments (returns `NULL` for any other arity). `ATTR_PROTECTED`. Its resolution
is `$TimeUnit` (`1/CLOCKS_PER_SEC`). The value is non-deterministic (it grows
with work done); it is a session CPU clock, not a numeric kernel, so there is no
packed/`Compile[]` surface.

- `Protected`.
- CPU time via `clock()`; does not advance during `Pause` or other idle waits.

**Attributes:** `Protected`.

## References

**See also:** [Pause](../../time-and-date/Pause/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`TimeUsed[]` gives the total CPU seconds consumed so far in the current session,
read from the same processor clock `Timing` brackets. Like `Timing`, it does
**not** advance during `Pause` or other idle waits — only work that actually
burns CPU moves it — and it counts CPU time summed over threads.

The value is non-deterministic (it grows with the work the session has done), so
examples show its structure rather than a fixed number. For elapsed real time
since the session began, including idle and `Pause` time, use `SessionTime`
instead; its resolution is `$TimeUnit`.
