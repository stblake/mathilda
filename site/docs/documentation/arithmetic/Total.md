# Total

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Total[list]`**

gives the total of the elements in list.

**`Total[list, n]`**

totals all elements down to level n.

**`Total[list, {n}]`**

totals elements at level n.

**`Total[list, {n1, n2}]`**

totals elements at levels n1 through n2.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Total[{a, b, c, d}]
Out[1]= a + b + c + d

In[2]:= Total[{{1, 2}, {3, 4}}]
Out[2]= {4, 6}

In[3]:= Total[{{1, 2}, {3, 4}}, 2]
Out[3]= 10

In[4]:= Total[{{1, 2}, {3}}, 2]
Out[4]= 6
```

### Applications (5)

```mathematica
In[1]:= Total[{1, 2, 3, 4}]
Out[1]= 10
```

On a matrix, the one-argument form sums the rows (a column total):

```mathematica
In[1]:= Total[{{1, 2}, {3, 4}, {5, 6}}]
Out[1]= {9, 12}
```

A level specification controls the depth of summation: `{2}` sums each column
instead, giving the per-column totals.

```mathematica
In[1]:= Total[{{1, 2}, {3, 4}, {5, 6}}, {2}]
Out[1]= {3, 7}
```

The tenth row of Pascal's triangle sums to a power of two:

```mathematica
In[1]:= Total[Table[Binomial[10, k], {k, 0, 10}]]
Out[1]= 1024
```

Summing exact rationals stays exact — a partial sum of the Basel series as a
single reduced fraction:

```mathematica
In[1]:= Total[Table[1/k^2, {k, 1, 100}]]
Out[1]= 1589508694133037873112297928517553859702383498543709859889432834803818131090369901/972186144434381030589657976672623144161975583995746241782720354705517986165248000
```

## Performance

Measured on arm64 Darwin at commit `2dea9cc05`.

| case | n | time |
|---|---:|---:|
| list of machine reals | 1,000 | 5 us |
| list of machine reals | 10,000 | 7 us |
| list of machine reals | 100,000 | 20 us |

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Clip to [0.25, 0.75] over 4x10^6 | 575 s | 1.95 s | 0.953 s |
| Dot 6x6 x 6x6 x 10000 | 338 s | 6.36 s | 4.01 s |
| Eigenvalues 600x600 (general) | 270 s | 163 s | 79.2 s |
| return {real, int, mask}, then Total | 61.2 s | 0.344 s | 0.983 s |
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| return {real, int}, then Total | 42 s | 0.219 s | 0.41 s |

## Implementation notes

`builtin_total` sums the elements of a list, optionally restricted to a level range. `Total[list]` sums the top level; `Total[list, n]` sums levels 1..n; `Total[list, {n}]` sums exactly level n; `Total[list, {n1, n2}]` sums a range; `Total[list, Infinity]` sums all levels. Negative level indices count from the bottom using the list's depth (`get_depth_for_total`). The chosen levels' elements are gathered and combined with `Plus` (so the usual numeric/symbolic Plus folding applies — `Total` is just structural element collection feeding `Plus`). `Total` carries `ATTR_PROTECTED`.

- `Protected`.
- `Total[list]` is equivalent to `Apply[Plus, list]`. The list elements are
  already evaluated, so the level-1 sum collects like terms directly (skipping a
  redundant re-evaluation and sort of every term); summing a large symbolic list
  is fast even when the result collapses to a few distinct terms.
- `Total[list, n]` totals all elements down to level `n`.
- `Total[list, {n}]` totals elements at level `n` only.
- Supports negative levels to count from the bottom (`-1` is the last dimension).
- Handles ragged arrays correctly by summing from the inside out when multiple levels are specified.
- `Total[list, Infinity]` totals all atoms in the expression.
- `Total[{}]` is `0` — the additive identity, at any level spec, matching
  `Plus @@ {}`. It answered `{}` until 2026-07-31, which then propagated as a
  non-number through anything consuming it (a `Select` that matched nothing is
  the usual way to reach it).
- Reads a packed buffer directly (`ndred_total`), exactly on an `int64` buffer:
  the accumulation abandons the whole result on overflow so the List path
  re-runs it and GMP answers.

**Attributes:** `Protected`.

## See also

[Select](../../data-structures/Select/)

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
