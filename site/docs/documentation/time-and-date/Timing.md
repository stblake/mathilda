# Timing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Timing[expr] evaluates expr, and returns a list of the time in seconds used, together with the result obtained.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_timing` brackets a single `evaluate(arg)` call with `clock()` (CPU time, `CLOCKS_PER_SEC`) and returns `{seconds, result}` as a two-element `List`, where `seconds` is `(end - start)/CLOCKS_PER_SEC` as an `EXPR_REAL`. It measures processor time, not wall-clock, and times a single evaluation only. (Note: the argument is evaluated explicitly inside the builtin; `Timing` is not given Hold attributes here.)

- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{timing, result}`.
- Includes only CPU time spent evaluating the expression.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)

## Notes & additional examples

### Worked examples

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

### Notes

`Timing` returns `{seconds, result}`. The first element is the CPU time spent and varies between runs and machines; only the second element (the computed result) is reproducible.
