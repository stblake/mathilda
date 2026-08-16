# UpTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UpTo[n]`**

is a symbolic specification that represents up to n objects or positions. If n objects or positions are available, all are used. If fewer are available, only those available are used.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= Take[{a, b, c, d}, UpTo[2]]
Out[1]= {a, b}

In[2]:= Take[{a, b}, UpTo[5]]
Out[2]= {a, b}

In[3]:= Take[Range[10], UpTo[100]]
Out[3]= {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
```

## Implementation notes

`UpTo[n]` is an inert specification object with no builtin handler — it is interpreted by the consumers that accept a count or position. List extractors (`Take`/`Drop` in `src/list.c`, `Part` ranges in `src/part.c`) detect `UpTo[n]` with an integer `n` and clamp the request to whatever is available: if at least `n` elements/positions exist all `n` are used, otherwise only those present, without raising the out-of-range error a bare `n` would. (Note: `SVD`'s `"UpTo"` target clamps to matrix rank, a separate use of the name.)

**Attributes:** none registered.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_findclusters_ndim.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findclusters_ndim.c)
- Tests: [`tests/test_findclusters_scalar_pin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findclusters_scalar_pin.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`UpTo[n]` is a count specification meaning "as many as `n`, but no error if fewer are available." With `Take`, requesting more elements than exist (Out[2], and the saturating example above) returns all of them rather than failing.
