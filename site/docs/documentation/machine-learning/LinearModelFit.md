# LinearModelFit

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`LinearModelFit[data] fits a linear model with an intercept and returns a PredictorFunction carrying its coefficients. Wolfram's version returns a FittedModel with regression diagnostics (RSquared, standard errors, ANOVA); those are not implemented and are not approximated. The coefficients are the same ones Predict finds.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [PredictorFunction](../../other-advanced/PredictorFunction/), [Predict](../../machine-learning/Predict/)

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)
