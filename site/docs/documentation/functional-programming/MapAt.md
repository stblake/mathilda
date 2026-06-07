# MapAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MapAt[f, expr, n]
    applies f to the element at position n in expr. Negative n counts from the end.
MapAt[f, expr, {i, j, ...}]
    applies f to the part of expr at position {i, j, ...}.
MapAt[f, expr, {{i1, j1, ...}, {i2, j2, ...}, ...}]
    applies f to the parts of expr at each of the listed positions.

Positions may contain All or Span specifications. MapAt[f, expr, 0]
applies f to the head of expr. Repeated positions apply f repeatedly.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, 2]
Out[1]= {a, f[b], c, d}

In[2]:= MapAt[f, {a, b, c, d}, {{1}, {4}}]
Out[2]= {f[a], b, c, f[d]}

In[3]:= MapAt[f, {{a, b, c}, {d, e}}, {2, 1}]
Out[3]= {{a, b, c}, {f[d], e}}

In[4]:= MapAt[f, {{a, b, c}, {d, e}}, {All, 2}]
Out[4]= {{a, f[b], c}, {d, f[e]}}

In[5]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, -3]
Out[5]= {{a, b, c}, h[{d, e}], f, g}

In[6]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {2, 1}]
Out[6]= {{a, b, c}, {h[d], e}, f, g}

In[7]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{2}, {1}}]
Out[7]= {h[{a, b, c}], h[{d, e}], f, g}

In[8]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{1, 1}, {2, 2}, {3}}]
Out[8]= {{h[a], b, c}, {d, h[e]}, h[f], g}
```

## Implementation notes

- `Protected`.
- Path components may be integers (positive or negative), `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- `MapAt[f, expr, {}]` applies `f` to the entire expression itself.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
