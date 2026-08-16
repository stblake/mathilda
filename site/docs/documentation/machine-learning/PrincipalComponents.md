# PrincipalComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrincipalComponents[matrix] gives the rows of matrix in principal-component coordinates, components ordered by decreasing variance. Rows are observations and columns are variables. Method -> "Correlation" standardises each variable to unit variance first, which is what you want when the columns have different units; the default "Covariance" does not.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
Out[1]= {{-2.82843, 0.0}, {-1.41421, 0.0}, {0.0, 0.0}, {1.41421, 0.0}, {2.82843, 0.0}}

In[2]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}]; {Variance[Map[First, p]], Variance[Map[Last, p]]}
Out[2]= {13.4817, 0.518295}
```

## Implementation notes

- Rows are observations, columns are variables. A flat list declines: one variable has
  no components to rotate.
- **The transform is an orthogonal rotation, so total variance is preserved** — it is
  redistributed into the leading components, not created or destroyed. Rank-deficient
  input therefore puts exactly zero variance in the trailing components: five points
  on a line give a second coordinate of `0`.
- `Method -> "Correlation"` standardises each variable to unit variance first, which
  is what you want when the columns have incommensurable units — otherwise the
  variable with the largest raw scale dominates for no statistical reason. The default
  `"Covariance"` does not. Note that "variance explained" means a different thing
  under each, so eigenvalues are not comparable across the two.
- An unrecognised `Method` returns unevaluated rather than silently choosing, so a
  typo is visible instead of quietly changing the statistics.
- **Eigenvector signs are canonical.** An eigenvector is defined only up to sign, and
  LAPACK's `dsyev` and the in-house Jacobi fallback do not agree on which they return;
  each component is flipped so its largest-magnitude loading is positive, so the
  output does not depend on whether the binary was linked against LAPACK.

**Attributes:** `Protected`.

## References

- Source: [`src/ml/pca.c`](https://github.com/stblake/mathilda/blob/main/src/ml/pca.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)
