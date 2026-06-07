# SingularValueDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SingularValueDecomposition[m]
    gives the singular value decomposition of a matrix m as a
    list {u, sigma, v}, where sigma is a diagonal matrix and
    m == u . sigma . ConjugateTranspose[v].  u and v have
    orthonormal columns.

SingularValueDecomposition[m, k]
    gives the SVD associated with the k largest singular values
    (or |k| smallest when k is negative).
SingularValueDecomposition[m, UpTo[k]]
    gives the SVD for as many of the k largest singular values
    as are available (up to MatrixRank[m]).

SingularValueDecomposition[{m, a}]
    gives the generalized singular value decomposition of m
    with respect to a as {{u, ua}, {sigma, sigma_a}, v} such
    that m == u . sigma . ConjugateTranspose[v] and
    a == ua . sigma_a . ConjugateTranspose[v].  Uses LAPACK
    dggsvd3 / zggsvd3 for real / complex machine-precision
    inputs; exact-numeric input is numericalised to 53 bits and
    routed through the same path.  High-precision MPFR input is
    currently downgraded to machine precision and emits the
    ::gmpdwn warning.  Free-symbolic input emits ::nogsymb and
    leaves the call unevaluated.

SingularValueDecomposition works on every input family supported
by the rest of the linear-algebra builtins:
    - exact integer / rational matrices (output stays exact;
      singular values are Sqrt[...] forms when irrational)
    - complex matrices (u and v are unitary in the Hermitian
      inner product)
    - machine-precision Real matrices (LAPACK dgesdd / zgesdd,
      or dggsvd3 / zggsvd3 for the generalized form)
    - arbitrary-precision MPFR matrices (one-sided Jacobi SVD
      at the input precision, real and complex)
    - free-symbolic matrices (eigendecomposition of m^H . m;
      the call is left unevaluated when no closed form exists)

Options:
    Tolerance -> t       :  zero out singular values below t
    TargetStructure ->
      "Dense"            :  u, sigma, v all dense (default)
      "Structured"       :  sigma returned as DiagonalMatrix[{..}]

A non-rank-2 or empty matrix emits
SingularValueDecomposition::matrix and the call is left
unevaluated.  Generalized SVD with mismatched column counts
emits ::matdims.  An out-of-range k or UpTo[k] emits ::sval.
Unknown option keys / values emit ::opts and leave the call
unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SingularValueDecomposition[{{1, 2}, {1, 2}}]
Out[1]= {{{1/Sqrt[2], 1/Sqrt[2]}, {1/Sqrt[2], -1/Sqrt[2]}}, {{Sqrt[10], 0}, {0, 0}}, {{1/Sqrt[5], -2/Sqrt[5]}, {2/Sqrt[5], 1/Sqrt[5]}}}

In[2]:= SingularValueDecomposition[{{1.1, 2, 5}, {3, -11, 4.2}}, 1]
Out[2]= {{{0.0195749}, {0.999808}}, {{12.1526}}, {{0.248586}, {-0.901763}, {0.353593}}}

In[3]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}, 30]]
Out[3]= {{{0.4045535848337569316424487226274, -0.9145142956773044526791769738123}, {0.9145142956773044526791769738107, 0.4045535848337569316424487226246}}, {{5.464985704219042650451188493292, 0.0}, {0.0, 0.3659661906262578204229643842617}}, {{0.5760484367663207913310985819436, 0.8174155604703632730886523884647}, {0.8174155604703632730886523884647, -0.5760484367663207913310985819436}}}

In[4]:= SingularValueDecomposition[{{1.0, 0}, {1.0, 10^-14}}, Tolerance -> 10^-15]
Out[4]= {{{-0.707107, -0.707107}, {-0.707107, 0.707107}}, {{1.41421, 0.0}, {0.0, 7.07107e-15}}, {{-1.0, -5e-15}, {-5e-15, 1.0}}}
```

## Implementation notes

- `Protected`.
- Options: `Tolerance -> t` zeroes out singular values with

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
