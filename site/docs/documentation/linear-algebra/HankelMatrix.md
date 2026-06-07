# HankelMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HankelMatrix[n] gives the n x n Hankel matrix with first row and column the integers 1..n.
HankelMatrix[{c1, ..., cm}] gives the m x m Hankel matrix with first column the given list.
HankelMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Hankel matrix with first column the first list and last row the second.
A Hankel matrix is constant along its antidiagonals; entries are copied verbatim.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HankelMatrix[4]
Out[1]= {{1, 2, 3, 4}, {2, 3, 4, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}}

In[2]:= HankelMatrix[{a, b, c, d}]
Out[2]= {{a, b, c, d}, {b, c, d, 0}, {c, d, 0, 0}, {d, 0, 0, 0}}

In[3]:= HankelMatrix[{x, y, z}, {z, a, b, c, d}]
Out[3]= {{x, y, z, a, b}, {y, z, a, b, c}, {z, a, b, c, d}}

In[4]:= HankelMatrix[{1, 1 + 2 I, 3 + 4 I}]
Out[4]= {{1, 1 + 2*I, 3 + 4*I}, {1 + 2*I, 3 + 4*I, 0}, {3 + 4*I, 0, 0}}

In[5]:= N[HankelMatrix[3]]
Out[5]= {{1.0, 2.0, 3.0}, {2.0, 3.0, 0.0}, {3.0, 0.0, 0.0}}
```

## Implementation notes

- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
