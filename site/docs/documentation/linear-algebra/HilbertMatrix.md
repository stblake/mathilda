# HilbertMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HilbertMatrix[n] gives the n x n Hilbert matrix with entries 1/(i + j - 1).
HilbertMatrix[{m, n}] gives the m x n Hilbert matrix.
Entries are exact Rationals unless the WorkingPrecision option requests MachinePrecision or a digit count.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HilbertMatrix[3]
Out[1]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}

In[2]:= HilbertMatrix[{3, 5}]
Out[2]= {{1, 1/2, 1/3, 1/4, 1/5}, {1/2, 1/3, 1/4, 1/5, 1/6}, {1/3, 1/4, 1/5, 1/6, 1/7}}

In[3]:= HilbertMatrix[3, WorkingPrecision -> MachinePrecision]
Out[3]= {{1.0, 0.5, 0.333333}, {0.5, 0.333333, 0.25}, {0.333333, 0.25, 0.2}}

In[4]:= Det[HilbertMatrix[3]]
Out[4]= 1/2160

In[5]:= Inverse[HilbertMatrix[3]]
Out[5]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}
```

## Implementation notes

- `Protected`.
- Entries are exact `Rational`s by default. The matrix is symmetric and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
