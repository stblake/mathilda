---
source: src/ml/dist.c
---
**Definition.** `NormalDistribution[mu, sigma]` represents a normal (Gaussian)
distribution with mean `mu` and standard deviation `sigma`; `NormalDistribution[]`
is the standard normal (`mu = 0`, `sigma = 1`). It is an inert, `Protected`
distribution head — no builtin of its own — consumed by `PDF` and `RandomVariate`.
Unlike a fitted `LearnedDistribution`, it prints its parameters in **full**, because
they are what the user specified. Its docstring is in `dist.c`.

**Representation.** The head stays symbolic until a distribution builtin reads it.
`ml_read_dist` (`src/ml/dist.c`) recognises `NormalDistribution`: with no arguments
it is the standard normal; with two it reads `mu` and `sigma` as real scalars into a
small `MlDist` record (kind `ML_D_NORMAL`). A non-positive `sigma` makes the read
*decline* rather than produce NaNs. `PDF` then evaluates `exp(-z^2/2)/(sigma
sqrt(2 pi))` with `z = (x - mu)/sigma`, and `RandomVariate` draws from the same
random stream as `RandomReal`.

**Usage & limits.** A specification, not a value — it carries no density on its own;
pass it to `PDF[dist, x]` (which threads over a list of points) or
`RandomVariate[dist, n]`. Parameters and the evaluation point must be numeric for
`PDF`/`RandomVariate`; a symbolic point leaves `PDF` unevaluated.
