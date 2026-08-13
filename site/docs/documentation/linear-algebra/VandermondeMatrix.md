# VandermondeMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VandermondeMatrix[{x1, ..., xn}] gives the n x n Vandermonde matrix with entry (i, j) equal to xi^(j-1).`**

**`VandermondeMatrix[{x1, ..., xn}, k] gives the n x k Vandermonde matrix.`**

<details>
<summary>Notes</summary>

The nodes need not be numerical or distinct; columns are successive powers, so the first column is all ones.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (4)

```mathematica
In[6]:= VandermondeMatrix[{1, 2, 3, 4}]
Out[6]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}, {1, 4, 16, 64}}

In[7]:= VandermondeMatrix[{a, b, c}]
Out[7]= {{1, a, a^2}, {1, b, b^2}, {1, c, c^2}}

In[8]:= Factor[Det[VandermondeMatrix[{a, b, c}]]]
Out[8]= (-a + b) (-a + c) (-b + c)

In[9]:= Factor[Det[VandermondeMatrix[{a, b, c, d}]]]
Out[9]= (-a + b) (-a + c) (-b + c) (-a + d) (-b + d) (-c + d)
```

## Options & behaviour

### Diagnostics

An empty node list, a non-list first argument, a non-positive integer `k`, the
single-argument structured-array conversion form `VandermondeMatrix[vmat]`
(unsupported — Mathilda has no structured-array representation), and any
over-arity call are returned unevaluated.

## Algorithm

VandermondeMatrix — a matrix whose rows are the successive powers of a sequence of nodes.

```text
  VandermondeMatrix[{x1, ..., xn}]      n x n Vandermonde matrix on the
                                        nodes x_i.
  VandermondeMatrix[{x1, ..., xn}, k]   n x k Vandermonde matrix.
```

The (1-based) entry (i, j) is x_i^(j-1), so the first column is all ones, the second column is the nodes themselves, the third their squares, and so

```text
on.  The nodes need not be numerical and need not be distinct: the entries
```

are built as Power[x_i, j-1] nodes (the first column emitted as the literal integer 1, so 0^0 reads as 1 to match the interpolation semantics) and the evaluator then simplifies them — numeric powers fold to their value, Power[x, 1] folds to x, leaving symbolic nodes as clean Power expressions.

Vandermonde matrices arise in polynomial interpolation and in computing moments in the monomial basis: LinearSolve[V, b] recovers the coefficients of the polynomial through the points {x_i, b_i}.

The single-argument structured-array conversion form, VandermondeMatrix[vmat], is not supported (Mathilda has no structured-array representation); a single matrix argument (a list of lists) is therefore left unevaluated.

Diagnostics mirror Wolfram's surface text:

```text
  - zero arguments  ->  VandermondeMatrix::argt
```

## Implementation notes

**Algorithm.** `builtin_vandermondematrix` constructs the matrix of successive powers of a node list: entry `(i,j) = x_i^(j−1)`. `VandermondeMatrix[{x1,...,xn}]` is `n × n`; `VandermondeMatrix[{x}, k]` is `n × k`. The builder `vm_build` emits each entry via `vm_entry`: the exponent-0 column is the literal Integer `1` (so `0^0` reads as 1, matching interpolation semantics), and every other entry is a `Power[x_i, j]` node which the evaluator later folds (numeric powers to their value, `Power[x,1]` to `x`), leaving symbolic nodes as clean `Power` expressions. Nodes are deep-copied and need not be numeric or distinct.

**Limits.** Zero arguments emit `VandermondeMatrix::argt`. The single-matrix structured-array conversion form `VandermondeMatrix[vmat]` is unsupported (Mathilda has no structured-array representation), so a list-of-lists argument (`vm_is_matrix`) is left unevaluated.

- `Protected`.
- The nodes need not be numerical and need not be distinct. Symbolic nodes stay
  as `Power` expressions; numeric powers fold to their value. Precision comes
  from the nodes themselves (e.g. ``2`20``) or from wrapping the result in `N`.
- The first column is the literal integer `1` (`xi^0`), so a zero node reads as
  `1` there rather than `Indeterminate` — matching the interpolation semantics.
- For distinct nodes the matrix is non-confluent and `LinearSolve[V, b]`
  recovers the coefficients `{a0, a1, ...}` of the polynomial
  `p(x) = a0 + a1 x + ...` through the points `{xi, bi}`. The determinant is the
  product of node differences `prod_{i<j} (xj - xi)`.
- An all-integer or all-machine-real node list writes the rank-2 result buffer
  directly. The integer path uses checked powers and falls back to the exact
  bignum path on `int64` overflow, so `VandermondeMatrix[{2, 100}, 20]` is still
  exact. (A *mixed* int+real list keeps its unpacked, two-head result, since the
  `Power` cells of an integer node stay exact beside a real one.) The
  single-vector form **compiles** (`Compile[{{v, _Real, 1}}, VandermondeMatrix[v]]`,
  rank 1 → rank 2). See [`packed-arrays.md`](../packed-arrays/index.md).

**Attributes:** `Protected`.

## References

**See also:** [Power](../../arithmetic/Power/), [N](../../arithmetic/N/)

- Source: [`src/linalg/vandermondemat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/vandermondemat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_vandermondematrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vandermondematrix.c)

## Notes & additional examples

### Notes

`VandermondeMatrix[{x1, ..., xn}]` gives the `n x n` matrix with entry `(i, j)` equal to `xi^(j-1)`; the two-argument form `VandermondeMatrix[nodes, k]` produces an `n x k` rectangular block. Nodes need not be numeric or distinct. The determinant vanishes exactly when two nodes coincide, which is why the Vandermonde system is invertible precisely for distinct interpolation points.
