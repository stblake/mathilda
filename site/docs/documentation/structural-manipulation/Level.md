# Level

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Level[expr, levelspec]
    gives a list of all subexpressions of expr on levels specified by levelspec.
Level[expr, levelspec, f]
    applies f to the sequence of subexpressions.

Level uses standard level specifications:
  n            levels 1 through n
  Infinity     levels 1 through Infinity
  {n}          level n only
  {n1, n2}     levels n1 through n2

Level[expr, {-1}] gives a list of all "atomic" objects in expr.
A positive level n consists of all parts of expr specified by n indices.
A negative level -n consists of all parts of expr with depth n.
Level 0 corresponds to the whole expression.
With the option setting Heads->True, Level includes heads of expressions and their parts.
Level traverses expressions in depth-first order, so that the subexpressions in the final list are ordered lexicographically by their indices.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Level[a + f[x, y^n], {-1}]
Out[1]= {a, x, y, n}

In[2]:= Level[a + f[x, y^n], 2]
Out[2]= {a, x, y^n, f[x, y^n]}

In[3]:= Level[x^2 + y^3, 3, Heads -> True]
Out[3]= {Plus, Power, x, 2, x^2, Power, y, 3, y^3}
```

## Implementation notes

- `Protected`.
- Default option: `Heads -> False`.
- Standard level specifications:
  - `n`: levels 1 through `n`.
  - `Infinity`: levels 1 through `Infinity`.
  - `{n}`: level `n` only.
  - `{n1, n2}`: levels `n1` through `n2`.
- Positive level `n` refers to distance from the top (level 0 is the whole expression).
- Negative level `-n` refers to distance from the bottom (depth `n`).
- Level `-1` corresponds to atomic objects.
- Lists subexpressions in post-order (depth-first), resulting in lexicographic ordering of indices.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
