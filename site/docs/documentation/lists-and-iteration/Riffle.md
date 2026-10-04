# Riffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Riffle[list, x]`**

Interleaves x into the gaps between successive elements of list, giving {e1, x, e2, x, ..., x, en}. Nothing is placed before the first or after the last element, so a list of length 0 or 1 comes back unchanged.

**`Riffle[list, {x1, x2, ...}]`**

Uses the xi cyclically, filling the n - 1 gaps left to right; separators beyond the last gap are unused. The head of list is preserved.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Riffle[{1, 2, 3}, 0]
Out[1]= {1, 0, 2, 0, 3}

In[2]:= Riffle[{a, b, c, d}, {x, y}]
Out[2]= {a, x, b, y, c, x, d}

In[3]:= Riffle[{a, b, c}, {x, y, z}]
Out[3]= {a, x, b, y, c}

In[4]:= Riffle[{a}, {x, y}]
Out[4]= {a}

In[5]:= Riffle[{a, b}, {}]
Out[5]= {a, b}

In[6]:= Riffle[f[a, b], x]
Out[6]= f[a, x, b]
```

### Applications (4)

One separator in every gap

```mathematica
In[7]:= Riffle[{1, 2, 3, 4}, 0]
Out[7]= {1, 0, 2, 0, 3, 0, 4}
```

Interleave a symbol between the elements

```mathematica
In[8]:= Riffle[Range[5], x]
Out[8]= {1, x, 2, x, 3, x, 4, x, 5}
```

A separator list cycles through the gaps

```mathematica
In[9]:= Riffle[{a, b, c, d, e}, {1, 2, 3}]
Out[9]= {a, 1, b, 2, c, 3, d, 1, e}
```

The comma-join idiom

```mathematica
In[10]:= StringJoin[Riffle[{"a", "b", "c"}, ", "]]
Out[10]= "a, b, c"
```

## Algorithm

Riffle — interleave separators into the gaps of a list.

Mathematica semantics:

```text
  Riffle[list, x]                 x is placed in every gap:
                                  Riffle[{1,2,3}, 0] -> {1, 0, 2, 0, 3}
  Riffle[list, {x1, ..., xk}]     the xi are consumed in order and cycle
                                  back to x1 after xk, filling the gaps
                                  left to right:
                                  Riffle[{a,b,c,d}, {x,y}] ->
                                    {a, x, b, y, c, x, d}
```

### The Gap Invariant

Separators go only BETWEEN consecutive elements — never before the first and never after the last. A list of n elements therefore has exactly n - 1 gaps, and the output has 2n - 1 slots. Two consequences drive the code below:

```text
  - n <= 1 means there are no gaps at all, so the result is the input
    unchanged whatever the separator is. This case is checked BEFORE the
    2n - 1 output sizing, because with n == 0 that expression underflows
    size_t to SIZE_MAX and the allocation would be nonsense.
  - separators past the last gap are never indexed, so
    Riffle[{a,b,c}, {x,y,z}] -> {a, x, b, y, c} simply never reaches z.
```

An empty separator list has nothing to interleave, so it also passes the list through unchanged; that check doubles as the guard that keeps the cycling index from dividing by zero.

The head of the first argument is preserved rather than forced to List, so Riffle[f[a,b], x] gives f[a, x, b]. That is also what makes Riffle[{}, 0] come back as {} with no special case.

### Performance

One pass, O(n) element copies, and a single exactly-sized allocation — the output length is known up front from n, so no growable buffer is needed.

NOT HANDLED: packed arrays (EXPR_NDARRAY) are a distinct representation from List and are left unevaluated here; see md-2aa.

## Implementation notes

**Algorithm.** `builtin_riffle` interleaves separators into the gaps of a list.
`Riffle[list, x]` places `x` in every gap; `Riffle[list, {x1, ..., xk}]` consumes
the `xi` in order and cycles back to `x1`, filling gaps left to right. A list of
`n` elements has exactly `n - 1` gaps, so the output has `2n - 1` slots: the gap
following element `i` takes separator index `i mod k`, and separators past the
last gap are simply never indexed.

**Edge invariants.** `n <= 1` (no gaps) or an empty separator list copies the
input through unchanged — checked *before* the `2n - 1` sizing, since with
`n == 0` that expression underflows `size_t`. The head of the first argument is
preserved rather than forced to `List`, so `Riffle[f[a, b], x]` gives
`f[a, x, b]` and `Riffle[{}, 0]` is `{}` with no special case.

**Data structures / limits.** One pass, O(n) element copies, a single
exactly-sized allocation (the output length is known up front). A packed/`NDArray`
first argument takes the `ndstruct_riffle` buffer fast path, falling back to
`ndstruct_delist_repack`; `Riffle` is on `pack.c`'s `AWARE` list. An atom first
argument stays unevaluated. No special attributes beyond `ATTR_PROTECTED`.

- `Protected`.
- Separators go only **between** consecutive elements — never before the first
  and never after the last. A list of $n$ elements has exactly $n - 1$ gaps, so
  the result has $2n - 1$ elements.
- A list of length 0 or 1 has no gaps, so it comes back **unchanged** whatever
  the separator is.
- With a `List` separator of length $k$, gap $i$ (1-based) receives
  `x[((i - 1) mod k) + 1]`. Separators beyond the last gap are simply unused, so
  `Riffle[{a, b, c}, {x, y, z}]` never places `z`.
- An **empty** separator list supplies nothing, so the list passes through
  unchanged.
- Only a `List` second argument cycles. Any other head is a single separator, so
  `Riffle[{a, b, c}, f[x, y]]` puts the whole `f[x, y]` in each gap.
- The object `list` need not have head `List`; its head is preserved on the
  result.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/)

- Source: [`src/list/riffle.c`](https://github.com/stblake/mathilda/blob/main/src/list/riffle.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`Riffle[list, x]` places `x` in each of the gaps between consecutive elements,
and `Riffle[list, {x1, ..., xk}]` cycles through the separators left to right. A
list of `n` elements has `n - 1` gaps, so separators go only *between* elements —
never before the first or after the last — and the output has `2n - 1` slots.

A single element (or an empty list) has no gaps, so the list is returned
unchanged. The head of the first argument is preserved, so
`Riffle[f[a, b], x]` gives `f[a, x, b]`. Interleaving a separator and then
`StringJoin`-ing is the usual way to build a delimited string.
