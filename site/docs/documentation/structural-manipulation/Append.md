# Append

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Append[expr, elem] adds elem to the end of expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Flatten[{{a, b}, {c, {d, e}}}]
Out[1]= {a, b, c, d, e}

In[2]:= Flatten[{{a, b}, {c, {d, e}}}, 1]
Out[2]= {a, b, c, {d, e}}

In[3]:= Flatten[f[f[a], b], -1, f]
Out[3]= f[a, b]
```

## Implementation notes

`builtin_append` (in `src/core.c`) requires a 2-arg call `Append[expr, elem]` whose first argument is an `EXPR_FUNCTION`. It allocates a fresh argument array one slot longer than `expr`, deep-copies every existing argument plus `elem` into it, and rebuilds a new `EXPR_FUNCTION` with the same head. Works on any head, not just `List`. Returns `NULL` (unevaluated) when the first argument is atomic.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
