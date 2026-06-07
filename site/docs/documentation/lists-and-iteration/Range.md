# Range

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Range[n]
    generates the list {1, 2, 3, ..., n}.
Range[n, m]
    generates the list {n, n + 1, ..., m - 1, m}.
Range[n, m, d]
    uses step d.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_range` (in `src/list.c`) generates the arithmetic sequence for `Range[imax]` (origin 1, step 1), `Range[imin, imax]`, and `Range[imin, imax, di]`. Bounds may be integers, reals, or rationals (parsed via `is_rational`); a `double` view of each is used only for the loop-termination test (`val <= max_val + 1e-14`, or the reversed test for negative step, with a 1,000,000-element cap). The element values themselves are built exactly: a running `curr_e` starts at `imin` and is advanced each step by `evaluate(Plus[curr_e, di])`, so integer and rational ranges stay exact while any real bound promotes the elements to `EXPR_REAL`. A zero step, or an empty oriented range, yields `{}`; the result is wrapped as `List[...]`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Range[5]
Out[1]= {1, 2, 3, 4, 5}
```

```mathematica
In[1]:= Range[2, 10, 2]
Out[1]= {2, 4, 6, 8, 10}
```

```mathematica
In[1]:= Range[0, 1, 1/4]
Out[1]= {0, 1/4, 1/2, 3/4, 1}
```

### Notes

`Range[n]` produces the integers `1` through `n`. `Range[n, m]` runs from `n` to
`m` in unit steps, and `Range[n, m, d]` uses step `d`. The step may be an exact
rational, in which case the result stays exact rather than being converted to
floating point. `Range` is the most direct way to build the index list that
`Map`, `Select`, or `Fold` then consume.
