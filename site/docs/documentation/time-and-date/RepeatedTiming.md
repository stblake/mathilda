# RepeatedTiming

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RepeatedTiming[expr] evaluates expr repeatedly and returns a list of the average time in seconds used, together with the result obtained.
RepeatedTiming[expr, t] does repeated evaluation for at least t seconds.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_repeated_timing` takes `(expr)` or `(expr, t)` where `t` (default 1.0 s) is the target measurement budget. It repeatedly `evaluate`s the expression, recording each per-run CPU time via `clock()`/`CLOCKS_PER_SEC` into a growable array, looping until at least 4 runs have completed *and* the accumulated wall budget `t` is reached. The first run's result is kept; later results are freed. It then `qsort`s the timings and returns a **trimmed mean**: it discards the lowest and highest quartile (`count/4` from each end) and averages the middle, returning `{averageSeconds, firstResult}` as a `List` of `EXPR_REAL` and the evaluated expression. Measures CPU time per run.

- `HoldFirst`, `Protected`, `SequenceHold`.
- Returns `{average_timing, result}`.
- Does repeated evaluation for at least `t` seconds. Default is 1 second.
- Gives a trimmed mean of the timings obtained, discarding lower and upper quartiles.
- Always evaluates `expr` at least four times.

**Attributes:** `HoldFirst`, `Protected`, `SequenceHold`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= RepeatedTiming[Sum[i, {i, 100}]]
Out[1]= {4.53806e-05, 5050}
```

### Notes

`RepeatedTiming` evaluates the expression many times and returns `{averageSeconds, result}`, giving a steadier estimate than `Timing` for fast operations. The timing value varies between runs; only the result element is reproducible.
