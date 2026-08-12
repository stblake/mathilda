# Riffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Riffle[list, x]
    Interleaves x into the gaps between successive
    elements of list, giving {e1, x, e2, x, ..., x, en}. Nothing is
    placed before the first or after the last element, so a list of
    length 0 or 1 comes back unchanged.
Riffle[list, {x1, x2, ...}]
    Uses the xi cyclically, filling the
    n - 1 gaps left to right; separators beyond the last gap are
    unused. The head of list is preserved.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

- `Protected`.
- Separators go only **between** consecutive elements — never before the first

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
