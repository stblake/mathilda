# LearnDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LearnDistribution[data] fits a distribution to data and returns a LearnedDistribution, usable with PDF. Method -> "Multinormal" is the default; Method -> "GaussianMixture" fits a mixture, choosing the component count by BIC. Multinormal fits a mean vector and a sample covariance (n-1 divisor, matching Variance). Rows are observations and columns are variables; a flat list is n observations of one variable. A singular covariance -- collinear columns, or fewer observations than dimensions -- returns unevaluated, because no density exists rather than because of an error. Method -> "ContingencyTable" is for NOMINAL data instead of numeric: it stores a probability per distinct outcome, which in one dimension is a categorical distribution. Outcomes may be any expressions -- strings, symbols, or equal-length lists of them -- compared structurally, and are kept in first-appearance order. Probabilities are empirical frequencies with no smoothing, so PDF of an outcome never observed is exactly 0; smoothing would require knowing how many outcomes were possible but unseen, which for arbitrary expressions is unknowable. Ragged outcomes decline.`**

## Examples (6)

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

### Options (2)

```mathematica
In[5]:= d = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]
Out[5]= LearnedDistribution["ContingencyTable", <>]

In[6]:= Last[LearnDistribution[{1.4, 1.4, 1.3, 1.5, 1.4, 1.7, 1.4, 1.5, 1.4, 1.5, 1.5, 1.6, 1.4, 1.1, 1.2, 4.7, 4.5, 4.9, 4., 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4., 4.7, 3.6, 6., 5.1, 5.9, 5.6, 5.8, 6.6, 4.5, 6.3, 5.8, 6.1, 5.1, 5.3, 5.5, 5., 5.1}, Method -> "GaussianMixture"]]
Out[6]= 2
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

**The spacing is measured between *distinct* points**, because a repeated value is not a
sample spacing of zero. Rounded data repeats most of its values — the first 15 iris petal
lengths per species, to 0.1 cm, have 45 values but only 28 distinct ones — and counting a
duplicate's zero distance to its twin put the median at 0, the floor at `1e-300`, and the
fit at nine components, two of them a single point with weight 1/45. Measured between
distinct values the floor is the data's real resolution (0.1² here) and BIC picks two
components, weights `{1/3, 2/3}` with means 1.42 and 4.91 — the same model Mathematica 15
learns. The floor is also held to at least `1e-4` of the average per-coordinate variance
(a component standard deviation of 1% of the data's), and the component count is capped
by the number of distinct points, as `FindClusters`' mixture path caps it.

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

**See also:** [LearnedDistribution](../../other-advanced/LearnedDistribution/), [PDF](../../machine-learning/PDF/), [PrincipalComponents](../../machine-learning/PrincipalComponents/), [Variance](../../data-structures/Variance/), [StandardDeviation](../../data-structures/StandardDeviation/), [FullForm](../../expression-information/FullForm/), [FindClusters](../../lists-and-iteration/FindClusters/)

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)
