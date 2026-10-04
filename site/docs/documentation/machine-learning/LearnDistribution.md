# LearnDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LearnDistribution[data] fits a distribution to data and returns a LearnedDistribution, usable with PDF. Method -> "Multinormal" is the default; Method -> "GaussianMixture" fits a mixture, choosing the component count by BIC. Multinormal fits a mean vector and a sample covariance (n-1 divisor, matching Variance). Rows are observations and columns are variables; a flat list is n observations of one variable. A singular covariance -- collinear columns, or fewer observations than dimensions -- returns unevaluated, because no density exists rather than because of an error. Method -> "ContingencyTable" is for NOMINAL data instead of numeric: it stores a probability per distinct outcome, which in one dimension is a categorical distribution. Outcomes may be any expressions -- strings, symbols, or equal-length lists of them -- compared structurally, and are kept in first-appearance order. Probabilities are empirical frequencies with no smoothing, so PDF of an outcome never observed is exactly 0; smoothing would require knowing how many outcomes were possible but unseen, which for arbitrary expressions is unknowable. Ragged outcomes decline.`**

## Examples (11)

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

### Applications (5)

A univariate normal from a flat list

```mathematica
In[7]:= m = LearnDistribution[{1., 2., 3., 4., 5., 6.}]
Out[7]= LearnedDistribution["Multinormal", <>]
```

Evaluate the fitted density at one point

```mathematica
In[8]:= PDF[m, {3.5}]
Out[8]= 0.213244
```

Nominal outcomes, not numbers

```mathematica
In[9]:= d = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]
Out[9]= LearnedDistribution["ContingencyTable", <>]
```

An outcome never observed has probability exactly 0

```mathematica
In[10]:= {PDF[d, "r"], PDF[d, "b"], PDF[d, "g"]}
Out[10]= {0.75, 0.25, 0.0}
```

BIC chooses two components for bimodal data

```mathematica
In[11]:= Last[LearnDistribution[{1.4, 1.4, 1.3, 1.5, 1.4, 1.7, 1.4, 1.5, 1.4, 1.5, 1.5, 1.6, 1.4, 1.1, 1.2, 4.7, 4.5, 4.9, 4., 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4., 4.7, 3.6}, Method -> "GaussianMixture"]]
Out[11]= 2
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

**Algorithm.** `builtin_learn_distribution` fits one of three families and returns a
`LearnedDistribution`, which prints elided (its parameters are derived) — the deliberate
opposite of a *specified* distribution such as `NormalDistribution[μ, σ]`.

- **`"Multinormal"`** (default). `ml_column_mean` and a sample covariance with the
  `n − 1` divisor (matching `Variance`, so a one-variable fit agrees with
  `StandardDeviation²`), factored by `ml_chol`. A singular covariance — collinear
  columns, or fewer observations than dimensions — declines, because no density exists.

- **`"GaussianMixture"`**. BIC model selection over `k = 1 … kmax`, each `k` fitted by
  `ml_gmm_fit` (EM in `src/ml/gmm.c`): deterministic farthest-first initialisation, a
  **log-space E-step** via log-sum-exp (a `dim`-factor density underflows for an
  outlying point, and a linear-space zero would hand it uniform responsibilities —
  silently refusing to distinguish exactly the informative points), and an M-step whose
  covariance floor is added as a **ridge** `floor · I` on the diagonal rather than a
  clamp. The variance floor itself (`ml_nn_floor`) is the **squared median
  nearest-neighbour distance between *distinct* points** — load-bearing, because a
  mixture likelihood is unbounded above and a merely-"small" floor lets BIC buy
  arbitrarily many near-singular spikes; counting duplicates' zero distances once drove
  the floor to `1e-300` and the fit to nine spurious components. A scale guard holds it
  to at least `1e-4` of the mean per-coordinate variance. `kmax` is bounded by
  `n/(dim+1)`, the distinct-point count, and 10; `BIC = paramCount · ln n − 2 · loglik`.

- **`"ContingencyTable"`** — nominal, not numeric. It runs *before* the numeric reader,
  builds a label vocabulary (`ml_labels_build`, `expr_eq`, first-appearance order) over
  outcomes that may be any expressions or equal-length lists of them, and stores
  frequency `count/n` per outcome. Probabilities are **empirical with no smoothing** —
  an unseen outcome is exactly `0`, since smoothing would require knowing the size of an
  unknowable outcome space. Ragged outcomes decline.

**Data structures.** An `MlGmm` (weights, means, covariances, Cholesky factors,
log-determinants) for the mixture; the stored `LearnedDistribution[method, payload, dim,
extra]` payload is a plain `List`.

**Complexity / limits.** Multinormal `O(n · dim² + dim³)`; the mixture search is
`O(kmax · iterations · n · k · dim²)` with `ml_nn_floor` an `O(n² · dim)` preprocess.
A flat list is accepted as `n` observations of one variable (a univariate normal),
unlike `PrincipalComponents`.

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

- A. P. Dempster, N. M. Laird and D. B. Rubin, *Maximum likelihood from incomplete data via the EM algorithm*, J. Roy. Statist. Soc. B **39** (1977) 1-38.
- G. Schwarz, *Estimating the dimension of a model*, Ann. Statist. **6** (1978) 461-464 (BIC).
- C. M. Bishop, *Pattern Recognition and Machine Learning* (Springer, 2006), §9.2 (Gaussian mixtures).
- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`LearnDistribution` fits a distribution and returns a `LearnedDistribution`, usable with
`PDF`. It prints elided — the deliberate opposite of a *specified* distribution such as
`NormalDistribution[μ, σ]`, whose parameters the user wrote and which therefore prints in
full; `FullForm` reveals the fitted parameters either way.

The default `"Multinormal"` fits a mean and a sample covariance with the `n − 1` divisor,
so a one-variable fit agrees with `StandardDeviation²`. A flat list is read as `n`
observations of one variable (unlike `PrincipalComponents`, which declines a single
variable). A singular covariance returns unevaluated, because no density exists.

`"GaussianMixture"` fits a mixture and chooses the component count by BIC — two
components, here, with the fitted means landing on the two modes. Its variance floor is
the squared median nearest-neighbour distance between *distinct* points, which stops BIC
buying arbitrarily many near-singular spikes. `"ContingencyTable"` is for nominal data:
probabilities are empirical frequencies with no smoothing, so an unseen outcome is
exactly `0`.
