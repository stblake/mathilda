# Max

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Max[x1, x2, ...]`**

yields the numerically largest of the xi.

**`Max[{x1, x2, ...}, {y1, ...}, ...]`**

yields the largest element of any of the lists.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

### Applications (5)

```mathematica
In[2]:= Max[3, 7, 2]
Out[2]= 7

In[3]:= Max[{1, 5}, {9, 2}]
Out[3]= 9

In[4]:= Max[2^100, 3^60, 5^40]
Out[4]= 1267650600228229401496703205376

In[5]:= Max[Abs[Eigenvalues[{{2, 1}, {1, 2}}]]]
Out[5]= 3

In[6]:= Max[x, 3, x]
Out[6]= Max[3, x]
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

**Algorithm.** `builtin_max` flattens any `List` arguments into a flat argument sequence, then
scans for the maximum among real-numeric atoms (compared with `expr_compare`) while collecting
distinct non-numeric/symbolic terms. `Infinity`/`-Infinity` and `Overflow[]` are handled as
absorbing/identity elements. If everything reduces to numbers it returns the single largest
value; otherwise it returns `Max[...]` over the numeric maximum plus the remaining symbolic
terms (returning `NULL` to stay unevaluated when nothing simplified). Empty `Max[]` is
`-Infinity`.

**Attributes:** `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Min](../../data-structures/Min/), [MinMax](../../data-structures/MinMax/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)

## Notes & additional examples

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
