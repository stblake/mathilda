# PredictorFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PredictorFunction[method, coefficients, featureCount] is the fitted object Predict returns. Apply it to a feature vector to get a prediction, or to "Method", "Coefficients" or "FeatureCount" to read it. A one-feature model also accepts a bare scalar.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)
