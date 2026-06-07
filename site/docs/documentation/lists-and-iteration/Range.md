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

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
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
