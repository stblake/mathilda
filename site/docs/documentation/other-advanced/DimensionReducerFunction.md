# DimensionReducerFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReducerFunction[method, {means, loadings...}, featureCount, reducedDimension] is the reusable reducer DimensionReduction returns. Apply it to a feature vector, or to a matrix of them, to project into the reduced space.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)
