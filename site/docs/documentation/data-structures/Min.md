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
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

## Implementation notes

**Algorithm.** `builtin_min` mirrors `Max`: it flattens `List` arguments, scans real-numeric
atoms for the minimum (via `expr_compare`), collects distinct symbolic terms, and treats
`Infinity`/`-Infinity`/`Overflow[]` as identity/absorbing elements. All-numeric input returns
the single smallest value; mixed input returns `Min[...]` over the numeric minimum and the
remaining symbolic terms, or `NULL` if nothing simplified. Empty `Min[]` is `Infinity`.

**Attributes:** `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Min[3, 7, 2]
Out[1]= 2
```

```mathematica
In[1]:= Min[1/3, 2/7, 5/11]
Out[1]= 2/7
```

```mathematica
In[1]:= Min[x, 0, Infinity]
Out[1]= Min[0, x]
```

```mathematica
In[1]:= Min[{}]
Out[1]= Infinity
```

### Notes

`Min[x1, x2, ...]` returns the numerically smallest argument, and `Min` of
several lists returns the smallest element across all of them. Comparisons are
exact, so rationals are ordered without converting to floating point —
`Min[1/3, 2/7, 5/11]` correctly picks `2/7`. With symbolic arguments `Min` stays
unevaluated but still prunes operands it can decide: `Min[x, 0, Infinity]` drops
`Infinity` (which can never be a minimum) and returns `Min[0, x]`. The empty case
`Min[{}]` returns `Infinity`, the identity element of minimisation — the value
that leaves any subsequent `Min` unchanged.
