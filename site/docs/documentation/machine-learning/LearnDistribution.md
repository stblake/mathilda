# LearnDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LearnDistribution[data] fits a distribution to data and returns a LearnedDistribution, usable with PDF. Method -> "Multinormal" is the default; Method -> "GaussianMixture" fits a mixture, choosing the component count by BIC. Multinormal fits a mean vector and a sample covariance (n-1 divisor, matching Variance). Rows are observations and columns are variables; a flat list is n observations of one variable. A singular covariance -- collinear columns, or fewer observations than dimensions -- returns unevaluated, because no density exists rather than because of an error. Method -> "ContingencyTable" is for NOMINAL data instead of numeric: it stores a probability per distinct outcome, which in one dimension is a categorical distribution. Outcomes may be any expressions -- strings, symbols, or equal-length lists of them -- compared structurally, and are kept in first-appearance order. Probabilities are empirical frequencies with no smoothing, so PDF of an outcome never observed is exactly 0; smoothing would require knowing how many outcomes were possible but unseen, which for arbitrary expressions is unknowable. Ragged outcomes decline.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {PDF[d, "r"], PDF[d, "b"], PDF[d, "g"]}
Out[1]= {0.75, 0.25, 0.0}
```

### Scope (3)

```mathematica
In[2]:= m = LearnDistribution[{1., 2., 3., 4., 5., 6.}]
Out[2]= LearnedDistribution["Multinormal", <>]

In[3]:= {PDF[m, {3.5}], PDF[NormalDistribution[Mean[{1.,2.,3.,4.,5.,6.}], StandardDeviation[{1.,2.,3.,4.,5.,6.}]], 3.5]}
Out[3]= {0.213244, 0.213244}

In[4]:= PDF[LearnDistribution[{{1.,2.},{2.,3.},{3.,5.},{4.,4.},{5.,7.},{6.,8.}}], {{3.5, 4.8}, {50., 50.}}]
Out[4]= {0.113186, 3.6193e-169}
```

### Options (1)

```mathematica
In[5]:= d = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]
Out[5]= LearnedDistribution["ContingencyTable", <>]
```

## Options & behaviour

**`"GaussianMixture"`** fits a mixture and chooses the component count by BIC — one
component for unimodal data, two for bimodal, with the fitted means landing on the modes.

**Its variance floor is the squared median nearest-neighbour distance, and it is
load-bearing.** A mixture's likelihood is *unbounded above*: a component collapsing onto
a single point drives its density, and hence the likelihood, to infinity. With a floor set
merely "small", the BIC search buys arbitrarily many near-singular spikes — measured in the
clustering path before its floor existed, six components for eight points. The median
nearest-neighbour distance says the honest thing instead: structure finer than the spacing
between samples is not resolvable. The *median* rather than the mean, so one tight pair
cannot drag the floor toward zero and reopen the same hole.

**A one-component mixture relates to the Multinormal fit exactly**, not approximately, and
the relationship is worth stating because it looks like a discrepancy:

EM maximises the likelihood, so its covariance uses the ML (`n`) divisor, while
`"Multinormal"` uses the unbiased (`n-1`) divisor to agree with `Variance`; the mixture
then adds the ridge. Both estimators are standard and both are correct for what they are.

**Verified against an independent implementation.** A one-dimensional fit reaches its
density through a Cholesky factor and a Mahalanobis distance, while
`PDF[NormalDdistribution[mu, sigma], x]` evaluates the scalar closed form; the two share
no code and agree exactly, including 3.5σ into the tail where a wrong log-determinant
would show as a small relative error rather than an obvious one. The density also
integrates to `1.0` over ±6σ, which pins the *normalisation* absolutely — an agreement
test alone would pass two densities that shared a normalisation error.

## Implementation notes

- Rows are observations, columns are variables; a flat list is `n` observations of one
  variable, which is a perfectly good univariate normal (unlike `PrincipalComponents`,
  which declines a single variable).
- The covariance uses the **`n - 1` divisor**, matching `Variance` and
  `StandardDeviation` — so a one-variable fit agrees with `StandardDeviation` squared.
- **A singular covariance returns unevaluated**: collinear columns, or fewer
  observations than dimensions, mean *no density exists*. A pseudo-inverse would invent
  one.
- **A fitted distribution prints elided**, the deliberate opposite of a *specified*
  distribution like `NormalDistribution[mu, sigma]`. A specified distribution's
  parameters are what the user wrote, so they are the information; a fitted one's are
  derived. `FullForm` reveals them either way.
- For a multinormal, `PDF[dist, {x1, …, xd}]` is **one** point — because the argument is
  itself a list — while a *matrix* threads to one density per row. That is the opposite
  reading from the scalar case, and it has to be.

**Attributes:** `Protected`.

## References

**See also:** [LearnedDistribution](../../other-advanced/LearnedDistribution/), [PDF](../../machine-learning/PDF/), [PrincipalComponents](../../machine-learning/PrincipalComponents/), [Variance](../../data-structures/Variance/), [StandardDeviation](../../data-structures/StandardDeviation/), [FullForm](../../expression-information/FullForm/)

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
