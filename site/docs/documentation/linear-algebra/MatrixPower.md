# MatrixPower

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MatrixPower[m, n]
    gives the n-th matrix power of the square matrix m.
    MatrixPower[m, n, v] gives the n-th matrix power of the matrix m applied to the vector v.
    When n is negative, MatrixPower finds powers of the inverse of the matrix m.
    MatrixPower[m, 0] gives IdentityMatrix[Length[m]].
    Fractional matrix powers are not currently supported.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MatrixPower[{{a, b}, {c, d}}, 2]
Out[1]= {{a^2 + b c, a b + b d}, {a c + c d, b c + d^2}}

In[2]:= MatrixPower[{{a, b}, {c, d}}, 2] == {{a, b}, {c, d}} . {{a, b}, {c, d}}
Out[2]= True

In[3]:= MatrixPower[{{a, b}, {c, d}}, -2] == Inverse[{{a, b}, {c, d}}] . Inverse[{{a, b}, {c, d}}]
Out[3]= True

In[4]:= MatrixPower[{{1, 1}, {1, 2}}, 10]
Out[4]= {{4181, 6765}, {6765, 10946}}

In[5]:= MatrixPower[{{1.2, 2.5, -3.2}, {0.7, -9.4, 5.8}, {-0.2, 0.3, 6.4}}, 5]
Out[5]= {{-1208.61, 19598.2, -12658.4}, {5784.51, -83315.1, 35420.6}, {-559.11, 1960.12, 11511.9}}

In[6]:= MatrixPower[{{2, 3, 0}, {4, 9, 0}, {0, 0, 4}}, 14]
Out[6]= {{25881337259836, 54508871401413, 0}, {72678495201884, 153068703863133, 0}, {0, 0, 268435456}}

In[7]:= MatrixPower[{{a, b}, {0, c}}, 4]
Out[7]= {{a^4, a^2 (a b + b c) + c^2 (a b + b c)}, {0, c^4}}

In[8]:= MatrixPower[{{1, 1}, {1, 2}}, 3, {1, 0}]
Out[8]= {5, 8}
```

## Implementation notes

**Algorithm.** `builtin_matrixpower` computes `MatrixPower[m, e]` (and `MatrixPower[m, e, v]`) for integer exponents by binary exponentiation (square-and-multiply). After validating that `m` is a non-empty square rank-2 tensor (`MatrixPower::matsq`) and that the exponent is an Integer or BigInt that fits in `int64_t`, it loops over the bits of `|e|`, accumulating the product and repeatedly squaring the running matrix; each matrix product goes through `dot2` (the `Dot` kernel, via the local `mat_dot` wrapper) followed by a full `evaluate` pass so entries simplify. `e == 0` returns `IdentityMatrix[n]`. Negative exponents first compute `Inverse[m]` (in `inv.c`) as the base; if `Inverse` returns unevaluated (singular), the call is abandoned. If a third vector argument is supplied, the final matrix is dotted with it.

**Data structures / limits.** Operates directly on `Expr*` `List`-of-`List` matrices; no dense numeric buffer. Fractional exponents (`Rational` or `Real`) are explicitly unsupported — they emit `MatrixPower::fract` and return unevaluated rather than going through an eigendecomposition; symbolic and BigInt-too-large exponents also return unevaluated. Cost is O(log|e|) matrix multiplies.

- `Protected`.
- `MatrixPower[m, n]` effectively evaluates the product of a matrix with itself `n` times.
- When `n` is negative, `MatrixPower` finds powers of the inverse of the matrix `m`.
- `MatrixPower[m, 0]` returns `IdentityMatrix[Length[m]]`.
- `MatrixPower` works only on square matrices.
- Uses binary exponentiation (repeated squaring) for efficient computation.
- Fractional matrix powers are not currently supported and generate a warning.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — powers of matrices.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — repeated multiplication and inversion.
- Source: [`src/linalg/matpow.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/matpow.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatrixPower[{{1, 1}, {0, 1}}, 3]
Out[1]= {{1, 3}, {0, 1}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 3}}, 2]
Out[1]= {{4, 0}, {0, 9}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 4}}, -1]
Out[1]= {{1/2, 0}, {0, 1/4}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 3}}, 2, {1, 1}]
Out[1]= {4, 9}
```

### Notes

`MatrixPower[m, n]` computes the `n`-th power by repeated matrix multiplication. A negative exponent (third example) raises the inverse to the corresponding positive power, so it requires a non-singular matrix; `MatrixPower[m, 0]` returns the identity matrix of the right size. The three-argument form `MatrixPower[m, n, v]` applies the matrix power to a vector — equivalent to `MatrixPower[m, n] . v` but without forming the full power matrix. Fractional exponents are not currently supported.
