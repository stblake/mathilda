# NormalDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NormalDistribution[mu, sigma] represents a normal distribution; NormalDistribution[] is the standard normal. Unlike a fitted model it prints its parameters in full, because they are what the user specified rather than an implementation detail.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
