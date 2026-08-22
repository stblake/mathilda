# RandomVariate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomVariate[dist] draws one value from dist; RandomVariate[dist, n] draws a list of n. Supports NormalDistribution[mu, sigma] and UniformDistribution[{lo, hi}], each also usable with no arguments for the standard case. Draws come from the same stream as RandomReal, so SeedRandom makes them reproducible. A non-positive standard deviation, or an inverted range, returns unevaluated rather than producing NaNs.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- **Draws come from the same stream as `RandomReal`**, so `SeedRandom` makes them
  reproducible. A sampler with its own generator would silently ignore `SeedRandom`
  while `RandomReal` honoured it — reproducibility half-working is worse than not
  working.
- Normal deviates use Box–Muller in its polar form, which needs no `sin`/`cos`. There
  was no Gaussian deviate anywhere in the tree before this.
- **A non-positive standard deviation, or an inverted range, returns unevaluated** —
  not `NaN`, which would propagate silently through a whole sample and surface much
  later as a strange plot. `RandomVariate[dist, 0]` is a valid empty request.

**Attributes:** `Protected`.

## References

**See also:** [RandomReal](../../random-number-generation/RandomReal/), [SeedRandom](../../random-number-generation/SeedRandom/)

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
