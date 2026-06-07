# Fit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Fit[data, {f1, ..., fn}, vars] finds a least-squares fit a1 f1 + ... + an fn to data for functions of the variables vars (a symbol or list of symbols).
Fit[{m, v}] gives the coefficient vector minimizing ||m.a - v|| for design matrix m and response vector v.
Data may be a list of values {v1, ...} (coordinates 1, 2, ...), univariate pairs {{x, v}, ...}, or multivariate rows {{x, y, ..., v}, ...}.
Options: WorkingPrecision (Automatic | n | Infinity), FitRegularization ({"Tikhonov"|"L2"|"RidgeRegression"|"LASSO"|"L1", lambda}), NormFunction.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
