# UpTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
UpTo[n]
    is a symbolic specification that represents up to n objects or positions. If n objects or positions are available, all are used. If fewer are available, only those available are used.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`UpTo[n]` is an inert specification object with no builtin handler — it is interpreted by the consumers that accept a count or position. List extractors (`Take`/`Drop` in `src/list.c`, `Part` ranges in `src/part.c`) detect `UpTo[n]` with an integer `n` and clamp the request to whatever is available: if at least `n` elements/positions exist all `n` are used, otherwise only those present, without raising the out-of-range error a bare `n` would. (Note: `SVD`'s `"UpTo"` target clamps to matrix rank, a separate use of the name.)

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Take[{a, b, c, d}, UpTo[2]]
Out[1]= {a, b}

In[2]:= Take[{a, b}, UpTo[5]]
Out[2]= {a, b}
```

### Notes

`UpTo[n]` is a count specification meaning "as many as `n`, but no error if fewer are available." With `Take`, requesting more elements than exist (Out[2]) returns all of them rather than failing.
