# Through

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Through[p[f1, f2, ...][x1, x2, ...]]
    distributes the trailing argument list across the inner functions,
    giving p[f1[x1, x2, ...], f2[x1, x2, ...], ...].
Through[expr, h]
    distributes only when the outer head equals h.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Through[{f, g, h}[x]]
Out[1]= {f[x], g[x], h[x]}

In[2]:= Through[(f + g)[x, y]]
Out[2]= f[x, y] + g[x, y]
```

## Implementation notes

- `Protected`.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
