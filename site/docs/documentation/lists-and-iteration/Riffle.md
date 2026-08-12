# Riffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Riffle[list, x]`**

Interleaves x into the gaps between successive elements of list, giving {e1, x, e2, x, ..., x, en}. Nothing is placed before the first or after the last element, so a list of length 0 or 1 comes back unchanged.

**`Riffle[list, {x1, x2, ...}]`**

Uses the xi cyclically, filling the n - 1 gaps left to right; separators beyond the last gap are unused. The head of list is preserved.

## Examples (6)

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

## See also

[List](../../other-advanced/List/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
