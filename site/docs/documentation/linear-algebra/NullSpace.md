# NullSpace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NullSpace[m]
    gives a list of vectors that forms a basis for the null
    space of the matrix m (i.e. vectors v such that m . v == 0).
NullSpace[m, Method -> "<name>"]
    runs a specific elimination algorithm.  Accepted method
    names are the same as RowReduce / LinearSolve / Inverse:
      "Automatic"                 — alias for "DivisionFreeRowReduction" (default)
      "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan
      "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot
      "CofactorExpansion"         — identity-if-invertible (falls back to
                                     DivisionFreeRowReduction on singular /
                                     rectangular m)

NullSpace works on both numerical and symbolic matrices.  The
matrix m may be square or rectangular.  When m has full column
rank the result is the empty list `{}`.

Basis vectors are returned with the rightmost free column
first.  For exact integer / rational input each basis vector
is scaled to clear integer denominators, so the result is
integer-valued whenever the input is integer-valued.  For
symbolic input the basis vectors are left in their natural
rational form.

An unknown method name emits NullSpace::method and leaves the
call unevaluated.  A non-rank-2 / empty matrix emits
NullSpace::matrix and the call is left unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, -2, 1}}

In[2]:= NullSpace[{{a, b}, {2 a, 2 b}}]
Out[2]= {{-b/a, 1}}

In[3]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[3]= {}

In[4]:= NullSpace[{{a, b, c}, {c, b, a}, {0, 0, 0}}]
Out[4]= {{1, -(a/b + c/b), 1}}

In[5]:= NullSpace[{{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}}]
Out[5]= {{12, -23, -5, 5}}

In[6]:= NullSpace[IdentityMatrix[5]]
Out[6]= {}

In[7]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[7]= {0, 0, 0}
```

## Implementation notes

- `Protected`.
- Returns a list of linearly-independent vectors whose span equals

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
