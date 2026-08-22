# LearnedDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LearnedDistribution[method, parameters, dimension, extra] is the fitted distribution LearnDistribution returns. Use it with PDF. Unlike a SPECIFIED distribution such as NormalDistribution[mu, sigma], it prints elided, because its parameters are derived rather than user-supplied.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
