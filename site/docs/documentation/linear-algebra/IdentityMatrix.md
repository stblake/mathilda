# IdentityMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IdentityMatrix[n] gives the n x n identity matrix.
IdentityMatrix[{m, n}] gives the m x n identity matrix.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[2]:= IdentityMatrix[{2, 3}]
Out[2]= {{1, 0, 0}, {0, 1, 0}}
```

## Implementation notes

`builtin_identitymatrix` accepts either an integer `n` (square `n×n`) or a pair `{m, n}` of integers, and constructs a `List` of `List`s with `Integer` `1` on the main diagonal (`i == j`) and `0` elsewhere. Non-integer or malformed dimension specs are returned unevaluated (`expr_copy(res)`). The output is exact integer entries; no numeric or symbolic processing is involved.

- `Protected`.
- Generates exact integer outputs (`1` on main diagonal, `0` elsewhere).
- Will remain unevaluated if arguments are symbolic or negative.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/construct.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/construct.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

A two-element argument gives a rectangular identity (1s on the main diagonal,
0s elsewhere):

```mathematica
In[1]:= IdentityMatrix[{2, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}}
```

It is the multiplicative identity for matrix products — multiplying any matrix
by a conformant identity leaves it unchanged:

```mathematica
In[1]:= IdentityMatrix[4] . HilbertMatrix[4] == HilbertMatrix[4]
Out[1]= True
```

### Notes

`IdentityMatrix[n]` gives the `n x n` identity; `IdentityMatrix[{m, n}]` gives the `m x n` rectangular identity. Use it as a seed for matrix algebra (e.g. `MatrixPower`, characteristic-matrix constructions) and as the neutral element of `Dot`.
