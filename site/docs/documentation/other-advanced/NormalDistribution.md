# NormalDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NormalDistribution[mu, sigma] represents a normal distribution; NormalDistribution[] is the standard normal. Unlike a fitted model it prints its parameters in full, because they are what the user specified rather than an implementation detail.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A specified distribution prints its parameters in full

```mathematica
In[1]:= NormalDistribution[2, 3]
Out[1]= NormalDistribution[2, 3]
```

The standard-normal density at 0 is 1/Sqrt[2 Pi]

```mathematica
In[2]:= PDF[NormalDistribution[0, 1], 0]
Out[2]= 0.398942
```

And away from the mean

```mathematica
In[3]:= PDF[NormalDistribution[0, 1], 1]
Out[3]= 0.241971
```

NormalDistribution[] is the standard normal

```mathematica
In[4]:= PDF[NormalDistribution[], 0]
Out[4]= 0.398942
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`NormalDistribution[mu, sigma]` is a specification, not a value: it stays symbolic
and prints its parameters in full (unlike a fitted `LearnedDistribution`, which
prints elided). Pass it to `PDF[dist, x]` for the density `exp(-z^2/2)/(sigma
Sqrt[2 Pi])` with `z = (x - mu)/sigma`, or to `RandomVariate[dist, n]` for samples
(from the same stream as `RandomReal`, so `SeedRandom` makes them reproducible).
`NormalDistribution[]` is the standard normal. A non-positive `sigma`, or a
non-numeric evaluation point, leaves the consuming call unevaluated rather than
producing NaNs.
