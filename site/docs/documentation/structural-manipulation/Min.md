# Min

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Min[x1, x2, ...]
    yields the numerically smallest of the xi.
Min[{x1, x2, ...}, {y1, ...}, ...]
    yields the smallest element of any of the lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Min[9, 2]
Out[1]= 2

In[2]:= Min[{4, 1, 7, 2}]
Out[2]= 1

In[3]:= Max[Infinity, 5]
Out[3]= Infinity
```

## Implementation notes

**Algorithm.** `builtin_min` mirrors `Max`: it flattens `List` arguments, scans real-numeric
atoms for the minimum (via `expr_compare`), collects distinct symbolic terms, and treats
`Infinity`/`-Infinity`/`Overflow[]` as identity/absorbing elements. All-numeric input returns
the single smallest value; mixed input returns `Min[...]` over the numeric minimum and the
remaining symbolic terms, or `NULL` if nothing simplified. Empty `Min[]` is `Infinity`.

- `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.
- Flattens `List` arguments.
- `Min[]` returns `Infinity`.
- `Max[]` returns `-Infinity`.
- Handles `Infinity` and `-Infinity`.
- Simplifies numeric arguments to a single value.

**Attributes:** `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
