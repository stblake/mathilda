# Cross

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cross[a, b]
    gives the vector cross product of two length-3 vectors.
Cross[a1, a2, ..., a(n-1)]
    gives the generalized (n-1)-fold cross product in n dimensions,
    i.e. the unique vector orthogonal to all inputs whose components
    are the signed cofactor minors of the matrix [a1; a2; ...; en].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Cross[{1, 2, -1}, {-1, 1, 0}]
Out[1]= {1, 1, 3}

In[2]:= Cross[{1, Sqrt[3]}]
Out[2]= {-Sqrt[3], 1}

In[3]:= Cross[{3.2, 4.2, 5.2}, {0.75, 0.09, 0.06}]
Out[3]= {-0.216, 3.708, -2.862}

In[4]:= Cross[{1.3 + I, 2, 3 - 2 I}, {6. + I, 4, 5 - 7 I}]
Out[4]= {-2 - 6*I, 6.5 - 4.9*I, -6.8 + 2.0*I}

In[5]:= Cross[{1, 2, 3}, {4, 5}]
Out[5]= Cross[{1, 2, 3}, {4, 5}]
```

## Implementation notes

**Algorithm.** `builtin_cross` implements the generalised cross product. Given `m = n-1` input vectors of length `n` (each must be a `List` of equal length, exactly one less than the number of arguments — otherwise it emits `Cross::nonn1` and returns `NULL`), the `i`-th component of the result is the signed minor obtained by stacking the input rows and deleting column `i`. Each minor is evaluated by `laplace_det` (shared with `Det`), and a sign `(-1)^(m+i)` is applied by wrapping the determinant in `Times[-1, …]` and calling `eval_and_free`. The classic 3-vector case `Cross[u, v]` is the `m = 2`, `n = 3` instance.

**Data structures.** Each minor is assembled into a flat `Expr**` array of `m*m` element pointers (borrowed from the input vectors, not copied) and passed to `laplace_det` with an explicit column-index list. Entries flow through symbolically, so the result is exact/symbolic whenever the inputs are. The output is a `List` of `n` components.

- `Protected`.
- Returns the cross product or totally antisymmetric product of $n-1$ vectors of length $n$.
- Works for symbolic and numerical inputs.
- Outputs `Cross::nonn1` error message if inputs are not equal-length vectors or if the number of arguments is not one less than their length.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/cross.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/cross.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
