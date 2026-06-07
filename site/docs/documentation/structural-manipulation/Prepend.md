# Prepend

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Prepend[expr, elem] adds elem to the beginning of expr.
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

**Algorithm.** `builtin_prepend` builds a new function node with the same head whose first
element is a copy of the new element followed by copies of all original arguments. Two-argument
form only; returns `NULL` if the first argument is an atom. (`PrependTo` is the mutating
variant that writes the result back to a symbol's OwnValue.)

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
