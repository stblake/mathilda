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
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

## Implementation notes

**Algorithm.** `builtin_max` flattens any `List` arguments into a flat argument sequence, then
scans for the maximum among real-numeric atoms (compared with `expr_compare`) while collecting
distinct non-numeric/symbolic terms. `Infinity`/`-Infinity` and `Overflow[]` are handled as
absorbing/identity elements. If everything reduces to numbers it returns the single largest
value; otherwise it returns `Max[...]` over the numeric maximum plus the remaining symbolic
terms (returning `NULL` to stay unevaluated when nothing simplified). Empty `Max[]` is
`-Infinity`.

**Attributes:** `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Max[3, 7, 2]
Out[1]= 7
```

```mathematica
In[1]:= Max[{1, 5}, {9, 2}]
Out[1]= 9
```

```mathematica
In[1]:= Max[2^100, 3^60, 5^40]
Out[1]= 1267650600228229401496703205376
```

```mathematica
In[1]:= Max[Abs[Eigenvalues[{{2, 1}, {1, 2}}]]]
Out[1]= 3
```

```mathematica
In[1]:= Max[x, 3, x]
Out[1]= Max[3, x]
```

### Notes

`Max[x1, x2, ...]` returns the numerically largest argument, and `Max` of a
collection of lists returns the largest element across all of them. Comparisons
are exact: `Max[2^100, 3^60, 5^40]` resolves a contest between three large
bignums (with automatic GMP promotion) and returns `2^100`, the actual winner.
Because `Max` is variadic and flattens lists, it composes naturally with other
operations — `Max[Abs[Eigenvalues[...]]]` computes the spectral radius of a
matrix in one line. When some arguments are non-numeric symbols, `Max` keeps the
call symbolic but still simplifies what it can: it discards duplicate operands
and drops any argument provably smaller than another, so `Max[x, 3, x]` collapses
to `Max[3, x]`.
