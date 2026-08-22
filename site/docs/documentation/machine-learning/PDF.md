# PDF

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PDF[dist, x] gives the probability density of dist at x, and threads over a list of x. Supports NormalDistribution and UniformDistribution.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]
Out[1]= {0.981398, -0.56572, 1.34033, 0.402313, -0.964221}

In[2]:= SeedRandom[1]; s = RandomVariate[NormalDistribution[5., 2.], 20000]; {Mean[s], StandardDeviation[s]}
Out[2]= {5.00479, 1.99627}

In[3]:= PDF[NormalDistribution[], {-1., 0., 1.}]
Out[3]= {0.241971, 0.398942, 0.241971}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
