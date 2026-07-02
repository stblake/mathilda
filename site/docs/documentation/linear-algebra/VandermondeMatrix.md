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
Out[4]= -(a - b) (a - c) (b - c)

In[5]:= LinearSolve[VandermondeMatrix[{1, 2, 3}], {6, 11, 18}]
Out[5]= {3, 2, 1}
```

## Implementation notes

**Algorithm.** `builtin_vandermondematrix` constructs the matrix of successive powers of a node list: entry `(i,j) = x_i^(j−1)`. `VandermondeMatrix[{x1,...,xn}]` is `n × n`; `VandermondeMatrix[{x}, k]` is `n × k`. The builder `vm_build` emits each entry via `vm_entry`: the exponent-0 column is the literal Integer `1` (so `0^0` reads as 1, matching interpolation semantics), and every other entry is a `Power[x_i, j]` node which the evaluator later folds (numeric powers to their value, `Power[x,1]` to `x`), leaving symbolic nodes as clean `Power` expressions. Nodes are deep-copied and need not be numeric or distinct.

**Limits.** Zero arguments emit `VandermondeMatrix::argt`. The single-matrix structured-array conversion form `VandermondeMatrix[vmat]` is unsupported (Mathilda has no structured-array representation), so a list-of-lists argument (`vm_is_matrix`) is left unevaluated.

- `Protected`.
- The nodes need not be numerical and need not be distinct. Symbolic nodes stay

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/vandermondemat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/vandermondemat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= VandermondeMatrix[{1, 2, 3, 4}]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}, {1, 4, 16, 64}}
```

Symbolic nodes give the classic generic Vandermonde matrix, with successive powers along each row:

```mathematica
In[1]:= VandermondeMatrix[{a, b, c}]
Out[1]= {{1, a, a^2}, {1, b, b^2}, {1, c, c^2}}
```

Its determinant factors into the product of all pairwise node differences — the celebrated closed form `Product[xj - xi, {i < j}]`:

```mathematica
In[1]:= Factor[Det[VandermondeMatrix[{a, b, c}]]]
Out[1]= (-a + b) (-a + c) (-b + c)
```

The identity scales to higher dimensions, here giving all six pairwise differences for four nodes:

```mathematica
In[1]:= Factor[Det[VandermondeMatrix[{a, b, c, d}]]]
Out[1]= (-a + b) (-a + c) (-b + c) (-a + d) (-b + d) (-c + d)
```

### Notes

`VandermondeMatrix[{x1, ..., xn}]` gives the `n x n` matrix with entry `(i, j)` equal to `xi^(j-1)`; the two-argument form `VandermondeMatrix[nodes, k]` produces an `n x k` rectangular block. Nodes need not be numeric or distinct. The determinant vanishes exactly when two nodes coincide, which is why the Vandermonde system is invertible precisely for distinct interpolation points.
