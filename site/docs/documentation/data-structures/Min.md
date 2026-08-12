# Min

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Min[x1, x2, ...]`**

yields the numerically smallest of the xi.

**`Min[{x1, x2, ...}, {y1, ...}, ...]`**

yields the smallest element of any of the lists.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

### Applications (4)

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

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Clip to [0.25, 0.75] over 4x10^6 | 575 s | 1.95 s | 0.953 s |
| MapThread[Max] over 4x10^6 | 14.8 s | 692 s | 0.772 s |
| MapThread[Min] over 4x10^6 | 14.7 s | 687 s | 0.769 s |
| integer Mod over 4x10^6 | 3.88 s | 0.504 s | 3.28 s |
| a b + a over 4x10^6 | 0.754 s | 1.07 s | 1.41 s |
| a + b over 4x10^6 | 0.383 s | 0.516 s | 0.74 s |

## Implementation notes

**Algorithm.** `builtin_min` mirrors `Max`: it flattens `List` arguments, scans real-numeric
atoms for the minimum (via `expr_compare`), collects distinct symbolic terms, and treats
`Infinity`/`-Infinity`/`Overflow[]` as identity/absorbing elements. All-numeric input returns
the single smallest value; mixed input returns `Min[...]` over the numeric minimum and the
remaining symbolic terms, or `NULL` if nothing simplified. Empty `Min[]` is `Infinity`.

**Attributes:** `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## See also

[Max](../../data-structures/Max/), [MinMax](../../data-structures/MinMax/)

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Min[x1, x2, ...]` returns the numerically smallest argument, and `Min` of
several lists returns the smallest element across all of them. Comparisons are
exact, so rationals are ordered without converting to floating point —
`Min[1/3, 2/7, 5/11]` correctly picks `2/7`. With symbolic arguments `Min` stays
unevaluated but still prunes operands it can decide: `Min[x, 0, Infinity]` drops
`Infinity` (which can never be a minimum) and returns `Min[0, x]`. The empty case
`Min[{}]` returns `Infinity`, the identity element of minimisation — the value
that leaves any subsequent `Min` unchanged.
