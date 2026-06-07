# PseudoInverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PseudoInverse[m]
    finds the Moore-Penrose pseudoinverse of a rectangular matrix m.
PseudoInverse[m, Tolerance -> t]
    specifies that singular values smaller than t times the maximum
    singular value should be dropped.  With the default setting
    Tolerance -> Automatic, the rationalisation precision of the
    input is used (Real -> 53 bits, MPFR -> input precision).

For non-singular square matrices m, the pseudoinverse coincides
with the standard inverse: PseudoInverse[m] == Inverse[m].

PseudoInverse works on exact (Integer / Rational / Complex)
matrices and on approximate (Real / MPFR) matrices.  For exact
input the result is exact; for inexact input the input is
rationalised, the pseudoinverse is computed in exact arithmetic
via a full-rank decomposition, and the result is numericalised
back to the input precision.

Algorithm: row-reduce m to identify rank r and a full-rank
decomposition m = B . C with B m x r and C r x n.  Then
    PseudoInverse[m] = ConjugateTranspose[C] . Inverse[C . ConjugateTranspose[C]]
                                            . Inverse[ConjugateTranspose[B] . B]
                                            . ConjugateTranspose[B].
When m is the zero matrix the pseudoinverse is the corresponding
zero matrix of transposed shape.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PseudoInverse[{{1,2},{3,4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}

In[2]:= PseudoInverse[{{1,2},{3,4}}] == Inverse[{{1,2},{3,4}}]
Out[2]= True

In[3]:= PseudoInverse[{{1,2,3},{4,5,6},{7,8,9}}]
Out[3]= {{-23/36, -1/6, 11/36}, {-1/18, 0, 1/18}, {19/36, 1/6, -7/36}}

In[4]:= PseudoInverse[{{0,0,0},{0,0,0}}]
Out[4]= {{0, 0}, {0, 0}, {0, 0}}

In[5]:= PseudoInverse[{{2,3},{2,2},{3,1},{4,3}}]
Out[5]= {{-29/134, -2/67, 22/67, 17/134}, {49/134, 8/67, -21/67, -1/134}}

In[6]:= PseudoInverse[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {0, 0, 0}}]
Out[6]= {{-0.0385185, 0.0966633, 0.0}, {0.210183, -0.0659894, 0.0}, {0.117363, 0.0282303, 0.0}}
```

## Implementation notes

- `Protected`.
- Works on rectangular and rank-deficient matrices over the integers,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
