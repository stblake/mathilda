# Timing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Timing[expr] evaluates expr, and returns a list of the time in seconds used, together with the result obtained.`**

<details>
<summary>Notes</summary>

Timing reports CPU time summed over threads; use AbsoluteTiming to measure elapsed time.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= Timing[Sum[i, {i, 1000}]]
Out[1]= {0.000244, 500500}
```

The timing is non-deterministic, so extract the reproducible result with `Part`:

```mathematica
In[1]:= Timing[Sum[i, {i, 1, 1000000}]][[2]]
Out[1]= 500000500000
```

```mathematica
In[1]:= Timing[D[Tan[x]^10, x]][[2]]
Out[1]= 10 Sec[x]^2 Tan[x]^9
```

## Implementation notes

`builtin_timing` brackets a single `evaluate(arg)` call with `clock()` (CPU time, `CLOCKS_PER_SEC`) and returns `{seconds, result}` as a two-element `List`, where `seconds` is `(end - start)/CLOCKS_PER_SEC` as an `EXPR_REAL`. It measures processor time, not wall-clock, and times a single evaluation only. (Note: the argument is evaluated explicitly inside the builtin; `Timing` is not given Hold attributes here.)

- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{timing, result}`.
- Includes only CPU time spent evaluating the expression, **summed over threads**.
  A threaded NDArray path or a BLAS call is therefore over-reported by roughly the
  core count — use `AbsoluteTiming` to measure how long something actually took.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## See also

[HoldAll](../../expression-information/HoldAll/), [SequenceHold](../../expression-information/SequenceHold/), [AbsoluteTiming](../../time-and-date/AbsoluteTiming/)

## References

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)
- Tests: [`tests/test_ludecomposition_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ludecomposition_machine.c)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)

## Notes & additional examples

### Notes

`Timing` returns `{seconds, result}`. The first element is the CPU time spent and varies between runs and machines; only the second element (the computed result) is reproducible.
