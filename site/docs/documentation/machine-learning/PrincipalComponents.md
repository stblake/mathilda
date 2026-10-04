# PrincipalComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrincipalComponents[matrix] gives the rows of matrix in principal-component coordinates, components ordered by decreasing variance. Rows are observations and columns are variables. Method -> "Correlation" standardises each variable to unit variance first, which is what you want when the columns have different units; the default "Covariance" does not.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
Out[1]= {{-2.82843, 0.0}, {-1.41421, 0.0}, {0.0, 0.0}, {1.41421, 0.0}, {2.82843, 0.0}}

In[2]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}]; {Variance[Map[First, p]], Variance[Map[Last, p]]}
Out[2]= {13.4817, 0.518295}
```

### Applications (4)

Collinear points leave the second component at exactly 0

```mathematica
In[3]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
Out[3]= {{-2.82843, 0.0}, {-1.41421, 0.0}, {0.0, 0.0}, {1.41421, 0.0}, {2.82843, 0.0}}
```

Rows are observations, columns variables

```mathematica
In[4]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}]
Out[4]= {{-4.81235, -0.203256}, {-1.22355, -0.550395}, {-1.36609, 0.856616}, {3.77227, -0.754988}, {3.62972, 0.652023}}
```

Variance concentrates in the leading component

```mathematica
In[5]:= {Variance[Map[First, p]], Variance[Map[Last, p]]}
Out[5]= {13.4817, 0.518295}
```

Standardise to unit variance first when the units differ

```mathematica
In[6]:= PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}, Method -> "Correlation"]
Out[6]= {{-1.83134, -0.064173}, {-0.502674, -0.208145}, {-0.451939, 0.33347}, {1.36761, -0.301383}, {1.41835, 0.240232}}
```

## Implementation notes

**Algorithm.** `builtin_principal_components` reads an `n × dim` matrix (a flat list
declines — one variable has nothing to rotate) and calls `ml_pca`. That routine centres
the columns with `ml_standardize`, and for `Method -> "Correlation"` *also* rescales
each column to unit variance — the single flag that is the whole difference between a
covariance and a correlation analysis. It forms the symmetric `dim × dim` covariance
(`n − 1` divisor) and diagonalises it with `ml_sym_eigen_desc`, which calls LAPACK's
`dsyev` when linked and otherwise falls back to an in-house cyclic Jacobi sweep
(`ml_jacobi`, written here because the linalg symmetric solvers are static to their
units). Eigenpairs are then selection-sorted into **descending** eigenvalue order (so
the leading component carries the most variance regardless of which backend ran), and
each eigenvector is **sign-canonicalised** — flipped so its largest-magnitude loading is
positive — because an eigenvector is defined only up to sign and `dsyev` and Jacobi
disagree on which they return. Finally each centred row is projected onto the ordered
eigenvectors. An unrecognised `Method` returns unevaluated rather than silently choosing.

**Data structures.** Plain `double` buffers throughout: a centred copy of the data, the
`dim × dim` covariance, an eigenvalue vector and a row-major eigenvector matrix. The
transformed rows are returned as a plain `List` (`ml_list_matrix`).

**Complexity / limits.** `O(n · dim²)` to form the covariance and `O(dim³)` to
diagonalise (the Jacobi fallback runs up to 100 sweeps). Because the transform is an
orthogonal rotation, total variance is preserved — rank-deficient input puts exactly
zero variance in the trailing components. Eigenvalues are *not* comparable between the
Covariance and Correlation methods, since "variance explained" means a different thing
under each.

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

- I. T. Jolliffe, *Principal Component Analysis*, 2nd ed. (Springer, 2002).
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013), §8.5 (the cyclic Jacobi eigenvalue method).
- Source: [`src/ml/pca.c`](https://github.com/stblake/mathilda/blob/main/src/ml/pca.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)

## Notes & additional examples

### Notes

The transform is an orthogonal rotation, so total variance is preserved — it is
redistributed into the leading components, not created or destroyed. Rank-deficient
input therefore puts exactly zero variance in the trailing components, as the collinear
first example shows.

`Method -> "Correlation"` standardises each variable to unit variance first, which is
what you want when the columns have incommensurable units, otherwise the variable with
the largest raw scale dominates for no statistical reason; the default `"Covariance"`
does not. Eigenvalues are not comparable across the two.

Eigenvector signs are canonical — each component is flipped so its largest-magnitude
loading is positive — so the output does not depend on whether the binary was linked
against LAPACK or fell back to the in-house Jacobi solver.
