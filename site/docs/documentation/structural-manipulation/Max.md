# Max

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Max[x1, x2, ...]
    yields the numerically largest of the xi.
Max[{x1, x2, ...}, {y1, ...}, ...]
    yields the largest element of any of the lists.
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

**Algorithm.** `builtin_max` flattens any `List` arguments into a flat argument sequence, then
scans for the maximum among real-numeric atoms (compared with `expr_compare`) while collecting
distinct non-numeric/symbolic terms. `Infinity`/`-Infinity` and `Overflow[]` are handled as
absorbing/identity elements. If everything reduces to numbers it returns the single largest
value; otherwise it returns `Max[...]` over the numeric maximum plus the remaining symbolic
terms (returning `NULL` to stay unevaluated when nothing simplified). Empty `Max[]` is
`-Infinity`.

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
