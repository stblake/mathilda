# Range

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Range[n]`**

generates the list {1, 2, 3, ..., n}.

**`Range[n, m]`**

generates the list {n, n + 1, ..., m - 1, m}.

**`Range[n, m, d]`**

uses step d.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Range[5]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= Range[2, 10, 2]
Out[2]= {2, 4, 6, 8, 10}

In[3]:= Range[0, 1, 1/4]
Out[3]= {0, 1/4, 1/2, 3/4, 1}

In[4]:= Range[10, 1, -1]
Out[4]= {10, 9, 8, 7, 6, 5, 4, 3, 2, 1}

In[5]:= Map[#^2 &, Range[5]]
Out[5]= {1, 4, 9, 16, 25}

In[6]:= Total[Range[100]]
Out[6]= 5050
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Clip to [0.25, 0.75] over 4x10^6 | 575 s | 1.95 s | 0.953 s |
| return {real, int, mask}, then Total | 61.2 s | 0.344 s | 0.983 s |
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| return {real, int}, then Total | 42 s | 0.219 s | 0.41 s |
| return {real, int, mask}, discarded | 41.5 s | 0.244 s | 0.972 s |
| return ragged {n, 1000, 100}, then Total | 20.2 s | 0.033 s | 0.051 s |

## Implementation notes

`builtin_range` (in `src/list.c`) generates the arithmetic sequence for `Range[imax]` (origin 1, step 1), `Range[imin, imax]`, and `Range[imin, imax, di]`. Bounds may be integers, reals, or rationals (parsed via `is_rational`); a `double` view of each is used only for the loop-termination test (`val <= max_val + 1e-14`, or the reversed test for negative step, with a 1,000,000-element cap). The element values themselves are built exactly: a running `curr_e` starts at `imin` and is advanced each step by `evaluate(Plus[curr_e, di])`, so integer and rational ranges stay exact while any real bound promotes the elements to `EXPR_REAL`. A zero step, or an empty oriented range, yields `{}`; the result is wrapped as `List[...]`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_backtrack.c`](https://github.com/stblake/mathilda/blob/main/tests/test_backtrack.c)
- Tests: [`tests/test_clip.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clip.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Range[n]` produces the integers `1` through `n`. `Range[n, m]` runs from `n` to
`m` in unit steps, and `Range[n, m, d]` uses step `d`. The step may be an exact
rational, in which case the result stays exact rather than being converted to
floating point. `Range` is the most direct way to build the index list that
`Map`, `Select`, or `Fold` then consume.
