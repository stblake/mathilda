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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cross[{1,0,0},{0,1,0}]
Out[1]= {0, 0, 1}
```

```mathematica
In[1]:= Cross[{a1,a2,a3},{b1,b2,b3}]
Out[1]= {-a3 b2 + a2 b3, -(-a3 b1 + a1 b3), -a2 b1 + a1 b2}
```

```mathematica
In[1]:= Cross[{2,1,-1},{1,-1,2}]
Out[1]= {1, -5, -3}
```

```mathematica
In[1]:= Cross[{1,2}]
Out[1]= {-2, 1}
```

```mathematica
In[1]:= Cross[{1,2,3,4},{5,6,7,8},{9,10,11,13}]
Out[1]= {4, -8, 4, 0}
```

### Notes

`Cross[a, b]` is the usual 3-vector cross product, returned in fully symbolic cofactor form so it works for indeterminate components. The general `Cross[a1, ..., a(n-1)]` form gives the unique vector in *n* dimensions orthogonal to all *n*-1 inputs, computed as the signed cofactor minors of the matrix whose remaining row is the basis vectors — so `Cross[{1,2}]` in the plane returns the 90-degree rotation `{-2, 1}`, and three 4-vectors yield a fourth vector perpendicular to all of them.
