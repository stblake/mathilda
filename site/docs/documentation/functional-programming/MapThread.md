# MapThread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MapThread[f, {{a1, a2, ...}, {b1, b2, ...}, ...}]
    gives {f[a1, b1, ...], f[a2, b2, ...], ...}, applying f to
    corresponding elements of the lists.
MapThread[f, {e1, e2, ...}, n]
    applies f to the parts of the ei at level n.

The ei must all have the same shape down through level n. MapThread is
a generalization of Map to functions of several variables; it takes
the function and its argument lists separately, unlike Thread. Lists
of associations with identical keys thread over their values.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MapThread[f, {{a, b, c}, {x, y, z}}]
Out[1]= {f[a, x], f[b, y], f[c, z]}

In[2]:= MapThread[f, {{{a, b}, {c, d}}, {{u, v}, {s, t}}}, 2]
Out[2]= {{f[a, u], f[b, v]}, {f[c, s], f[d, t]}}

In[3]:= MapThread[Plus, {{a, b, c}, {u, v, w}, {x, y, z}}]
Out[3]= {a + u + x, b + v + y, c + w + z}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
