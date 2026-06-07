# VandermondeMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
VandermondeMatrix[{x1, ..., xn}] gives the n x n Vandermonde matrix with entry (i, j) equal to xi^(j-1).
VandermondeMatrix[{x1, ..., xn}, k] gives the n x k Vandermonde matrix.
The nodes need not be numerical or distinct; columns are successive powers, so the first column is all ones.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= VandermondeMatrix[{x1, x2, x3, x4}]
Out[1]= {{1, x1, x1^2, x1^3}, {1, x2, x2^2, x2^3}, {1, x3, x3^2, x3^3}, {1, x4, x4^2, x4^3}}

In[2]:= VandermondeMatrix[{2, 3, 5}]
Out[2]= {{1, 2, 4}, {1, 3, 9}, {1, 5, 25}}

In[3]:= VandermondeMatrix[{a, b, c}, 2]
Out[3]= {{1, a}, {1, b}, {1, c}}

In[4]:= Factor[Det[VandermondeMatrix[{a, b, c}]]]
Out[4]= (-a + b) (-a + c) (-b + c)

In[5]:= LinearSolve[VandermondeMatrix[{1, 2, 3}], {6, 11, 18}]
Out[5]= {3, 2, 1}
```

## Implementation notes

- `Protected`.
- The nodes need not be numerical and need not be distinct. Symbolic nodes stay

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
