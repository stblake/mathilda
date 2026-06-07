# Span

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
i;;j represents a span of elements i through j. i;;j;;k represents a span in steps of k.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= {a, b, c, d, e, f, g, h}[[2;;5]]
Out[1]= {b, c, d, e}

In[2]:= {a, b, c, d, e, f, g, h}[[1;;-1;;3]]
Out[2]= {a, d, g}

In[3]:= t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = x; t
Out[3]= {a, x, x, x, x, f, g, h}

In[4]:= t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = {p, q, r, s}; t
Out[4]= {a, p, q, r, s, f, g, h}

In[5]:= Range[10][[3;;All]]
Out[5]= {3, 4, 5, 6, 7, 8, 9, 10}
```

## Implementation notes

**Algorithm.** `Span` is not a standalone builtin — it is the head the parser produces for the
`i;;j` and `i;;j;;k` syntax (`OP_SPAN`, precedence 290, in `parse.c`). A `Span[start, end,
step]` expression is interpreted only inside `Part`/`Extract`: `expr_part` (in `part.c`)
resolves it against a concrete list, mapping negative endpoints to `len + k + 1`, clamping
`UpTo`/`All` endpoints, computing the element count from the step, and emitting the selected
elements. On its own a bare `Span` stays unevaluated.

- `m[[i;;j;;k]]` is equivalent to `Take[m, {i, j, k}]` but evaluated natively within `Part`.
- `m[[i;;j]] = v` can be used to assign `v` iteratively over a span of elements. If `v` is a list, elements are assigned sequentially. If `v` is a non-list expression, it is assigned uniformly to all elements in the span.
- When used in `Part`, negative `i` and `j` count from the end.
- `i` and `j` can be of the form `UpTo[n]` to restrict endpoints to the actual length of the list.
- Any argument of `Span[...]` can be `All`.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
