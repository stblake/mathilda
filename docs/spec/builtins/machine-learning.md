# Machine Learning

Kernels for this category live in [`src/ml/`](../../../src/ml/) rather than inside the
builtin that first needed them, and that is deliberate. Each is written against a
row-major `n × dim` buffer of machine doubles rather than against an `Expr`, so a
second consumer can reach it: `src/ml/gmm.c`'s EM fit serves both
`FindClusters[…, Method -> "GaussianMixture"]` and (in prospect) `LearnDistribution`,
and `src/ml/pca.c`'s column statistics serve `Standardize`, `PrincipalComponents`,
and any later feature scaling. `src/list/find_clusters.c`, by contrast, holds 55
static functions and exports exactly one symbol, which is why `ClusteringComponents`,
`ClusteringTree`, `Dendrogram` and `NearestNeighborGraph` cannot yet be built on it.

**Results are plain `List`s.** The machine bridge's `na_build_matrix` returns a
*visible* `NDArray`, whose head is `NDArray`, so a result built that way compares
`False` against the literal list a user would write — while `Inverse`, `Dot` and
`LinearSolve` compare `True`. These builtins construct `List`s and let the evaluator's
own packing gate decide whether the result is held as a buffer, keeping the surface
consistent with the rest of the system.

## Standardize

Shifts each column to zero mean and rescales it to unit standard deviation.
Attributes: `Protected`.

- `Standardize[data]`

**Features**:
- **Columns are variables, rows are observations.** A flat list is treated as `n`
  observations of *one* variable, not one observation of `n`.
- The divisor is `n - 1` (the sample standard deviation), matching
  `StandardDeviation` — so `Standardize[x]` agrees with
  `(x - Mean[x])/StandardDeviation[x]` written out by hand. A mismatch here would be
  invisible on the mean but not on the scale.
- **A constant column becomes exactly `0`, not `Indeterminate`.** Zero variance
  carries no information, so "no deviation from the mean" is the honest value;
  dividing by the zero standard deviation would propagate `Indeterminate` through
  every reduction over the row.

```mathematica
In[1]:= Standardize[{1., 2., 3., 4.}]
Out[1]= {-1.1619, -0.387298, 0.387298, 1.1619}

In[2]:= Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]
Out[2]= {{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}

In[3]:= Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]
Out[3]= {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}
```

## PrincipalComponents

Gives the rows of a matrix in principal-component coordinates, components ordered by
decreasing variance. Attributes: `Protected`.

- `PrincipalComponents[matrix]`
- `PrincipalComponents[matrix, Method -> "Covariance" | "Correlation"]`

**Features**:
- Rows are observations, columns are variables. A flat list declines: one variable has
  no components to rotate.
- **The transform is an orthogonal rotation, so total variance is preserved** — it is
  redistributed into the leading components, not created or destroyed. Rank-deficient
  input therefore puts exactly zero variance in the trailing components: five points
  on a line give a second coordinate of `0`.
- `Method -> "Correlation"` standardises each variable to unit variance first, which
  is what you want when the columns have incommensurable units — otherwise the
  variable with the largest raw scale dominates for no statistical reason. The default
  `"Covariance"` does not. Note that "variance explained" means a different thing
  under each, so eigenvalues are not comparable across the two.
- An unrecognised `Method` returns unevaluated rather than silently choosing, so a
  typo is visible instead of quietly changing the statistics.
- **Eigenvector signs are canonical.** An eigenvector is defined only up to sign, and
  LAPACK's `dsyev` and the in-house Jacobi fallback do not agree on which they return;
  each component is flipped so its largest-magnitude loading is positive, so the
  output does not depend on whether the binary was linked against LAPACK.

```mathematica
In[1]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
Out[1]= {{-2.82843, 0.0}, {-1.41421, 0.0}, {0.0, 0.0}, {1.41421, 0.0}, {2.82843, 0.0}}

In[2]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}];
        {Variance[Map[First, p]], Variance[Map[Last, p]]}
Out[2]= {13.4817, 0.518295}
```

## DimensionReduce

Reduces each row of a matrix to a lower-dimensional representation.
Attributes: `Protected`.

- `DimensionReduce[data, k]`
- `DimensionReduce[data, k, Method -> m]`

**Features**:
- The three methods are **one algorithm with three ways of forming the symmetric
  matrix to decompose**, which is why they share the eigendecomposition rather than
  each carrying its own linear algebra:

| Method | What it decomposes |
|---|---|
| `"PrincipalComponentsAnalysis"` (default) | the covariance of the **centred** columns |
| `"LatentSemanticAnalysis"` | the Gram matrix `X'X`, **without** centring — a truncated SVD |
| `"MultidimensionalScaling"` | the double-centred squared-distance matrix (classical Torgerson scaling) |

- **Skipping the centring is the entire difference between PCA and LSA.** A
  term-document matrix is sparse and non-negative, and centring destroys both
  properties along with the meaning of a zero entry — which is why LSA does not.
- **Reducing to `k` gives exactly the first `k` principal components**, not a
  separately-fitted `k`-component model: `DimensionReduce[data, 2]` equals
  `Map[Take[#, 2] &, PrincipalComponents[data]]`.
- **Classical MDS on Euclidean distances is the same embedding as PCA**, reached by a
  different route (an `n × n` double-centred distance matrix rather than a
  `dim × dim` covariance). That agreement is used as a cross-check on both in the test
  suite; it also means MDS earns its keep only when the distances come from somewhere
  other than the coordinates.
- `"MultidimensionalScaling"` is **capped at 2000 rows**, its matrix being `n × n` —
  the same order of ceiling, for the same reason, as `FindClusters`' `"Spectral"`.
- **Asking for more dimensions than the data supports returns unevaluated** rather
  than padding with zeros, since padding would look like a successful reduction to a
  caller checking only the shape. An unknown `Method`, a non-positive `k`, an omitted
  `k`, and a flat list all decline too.

**Not implemented**: Wolfram's `DimensionReduce` can also choose `k` itself and can
return a `DimensionReducerFunction` applicable to *new* data. The second is the
substantive gap — a reusable reducer is a trained model, and that representation is
being designed with the `Predict` family rather than invented twice.

```mathematica
In[1]:= d = {{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}};
        DimensionReduce[d, 2]
Out[1]= {{-5.28583, 0.302136}, {-1.6734, -0.58098}, {0.141634, 3.22424},
         {1.80921, -4.66364}, {5.00839, 1.71825}}

In[2]:= DimensionReduce[d, 2, Method -> "LatentSemanticAnalysis"]
Out[2]= {{3.53625, 1.0399}, {7.03311, -0.290957}, {9.24375, 3.24698},
         {9.903, -4.78819}, {13.8764, 1.13662}}
```

## Trained models

`Predict` returns a **`PredictorFunction`** — a fitted object that can be stored in a
variable and applied to new inputs. The representation is deliberately plain: an
ordinary expression whose head is the symbol `PredictorFunction` and whose arguments
are ordinary lists and strings.

Two richer designs were considered and rejected, which is worth knowing because the
obvious one is the rejected one:

- **A new opaque node type**, following `CompiledFunction`. That is the codebase's
  precedent for a callable carrying binary state — but a compiled function is a VM
  program (bytecode, a register file, a reference count), whereas a fitted linear model
  is a short vector of numbers. A new node type would mean new cases in `expr_copy`,
  `expr_free`, `expr_eq`, `expr_hash`, `expr_compare` and the printer, to store what
  existing expression types already hold. Machine precision is not an argument for it
  either: a packed list already holds a dense buffer and is still a list.
- **An `Association` payload**, the most Wolfram-ish shape. It buys nothing positional
  arguments do not, while adding a dependency to every read. Named properties are
  reached through the application instead, which is where you would ask for them:
  `p["Coefficients"]`, `p["Method"]`, `p["FeatureCount"]`.

**Fitted models print abbreviated**, showing their method and eliding their parameters:

```mathematica
In[1]:= Predict[{1. -> 3., 2. -> 5., 3. -> 7.}]
Out[1]= PredictorFunction["LinearRegression", <>]

In[2]:= FullForm[%]
Out[2]= PredictorFunction["LinearRegression", List[1.0, 2.0], 1, 0]
```

This is not cosmetic. A `"NearestNeighbors"` predictor carries its *entire training set*
as its parameter block, so the unabridged form is unreadable at any real size and scrolls
the useful answer off the screen. `InterpolatingFunction` already elides its data the same
way, and for the same reason. Eliding is safe precisely because `FullForm` still reveals
everything — the information is one keystroke away rather than gone.

Applying a model uses the evaluator's existing composite-head dispatch — the same
mechanism behind `Function[…][args]` and `<|…|>[key]` — so a trained model introduces
no new evaluation concept.

## Predict

Fits a predictor to data and returns a `PredictorFunction`. Attributes: `Protected`.

- `Predict[data]`
- `Predict[data, Method -> "LinearRegression"]`
- `Predict[data, Method -> "NearestNeighbors"]`
- `Predict[data, Method -> {"NearestNeighbors", "NeighborsNumber" -> k}]`

**Features**:
- Data is either a list of rules `{features -> value, …}` or a **matrix whose last
  column is the response** — the shape data usually arrives in. Both give the same
  model; which was used is not remembered.
- A one-feature model accepts a bare scalar as well as a one-element list, so
  `p[3.]` works for a single-variable regression.
- **A singular system returns unevaluated.** Perfectly collinear features, or fewer
  observations than parameters, have no unique fit — a pseudo-inverse would return one
  of infinitely many answers and *look* like a successful fit.
- `"LinearRegression"` is the only method implemented and the default. Any other
  declines rather than silently linear-regressing.
- A wrongly-shaped input leaves the application unevaluated rather than guessing.
- **`"NearestNeighbors"` has no fitted parameters — the training set *is* the model.**
  That is why it costs nothing to fit and everything to apply. The prediction is the
  mean response of the `k` nearest training rows; `k` defaults to 3, enough to average
  away one noisy response and few enough to stay local, and is clamped to the number of
  training rows. It answers `"NeighborCount"` and `"TrainingData"` as well as the shared
  properties.
- `"NeighborsNumber"` on a `"LinearRegression"` is **refused, not ignored** — silently
  accepting a meaningless option would hide a real mistake. So is an unrecognised
  sub-option, and a non-positive `k`.

**On verification**: the 1-neighbour predictor is cross-checked against the existing
independent `Nearest` builtin, on both sides of a midpoint, so neighbour *selection* is
validated against separate code. `Nearest` supports only the one-neighbour scalar form
here (its `k` form and point form decline), so that check covers selection at `k = 1`
and says nothing about the averaging, which is asserted directly instead.

```mathematica
In[1]:= p = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]
Out[1]= PredictorFunction["LinearRegression", {1.0, 2.0}, 1]

In[2]:= {p[10.], p["Coefficients"], p["Method"]}
Out[2]= {21.0, {1.0, 2.0}, "LinearRegression"}

In[3]:= Predict[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.},
                 {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Out[3]= {1.0, 2.0, 3.0}
```

## LinearModelFit

- `LinearModelFit[data]`

Fits a linear model with an intercept and returns a `PredictorFunction` carrying its
coefficients — the same ones `Predict` finds. Attributes: `Protected`.

**Not implemented**: Wolfram's `LinearModelFit` returns a `FittedModel` whose
properties are regression diagnostics (`RSquared`, standard errors, ANOVA). Those are
a separate piece of work and are **not** approximated here.

## DimensionReduction

Returns a reusable **`DimensionReducerFunction`**. Attributes: `Protected`.

- `DimensionReduction[data, k]`

This is the counterpart to `DimensionReduce[data, k]`, and Wolfram's split of the two is
kept: `DimensionReduce` gives you the *reduced training data*, `DimensionReduction` gives
you a *reducer you can apply to data it was never trained on*.

**Features**:
- **New rows are centred on the training column means**, not on the incoming batch's.
  That is the entire point of a reusable reducer — it is what makes projections
  comparable across batches — and it is worth stating because the alternative bug is
  hard to see: centring a single new point against itself gives all zeros, which looks
  correct on data that happens to sit near the origin.
- Applies to a single feature vector or to a matrix of them, so a batch needs no `Map`.
- On its own training data it reproduces `DimensionReduce` exactly.
- Answers `"Method"`, `"FeatureCount"` and `"ReducedDimension"`.
- Only `"PrincipalComponentsAnalysis"` is available as a reducer so far; MDS has no
  out-of-sample extension without an explicit one (Nyström), and LSA's would need the
  term-document vocabulary carried along.

A reducer is the *same* representation as a `PredictorFunction` — positional and
method-tagged — which is why adding it needed no new evaluation machinery.

```mathematica
In[1]:= tr = {{100., 200.}, {102., 205.}, {104., 210.}, {106., 215.}, {108., 220.}};
        r = DimensionReduction[tr, 1]
Out[1]= DimensionReducerFunction["PrincipalComponentsAnalysis",
          {{104.0, 210.0}, {0.371391, 0.928477}}, 2, 1]

In[2]:= r[{110., 225.}]           (* a point NOT in the training set *)
Out[2]= {16.1555}

In[3]:= r[{104., 210.}]           (* exactly at the training mean *)
Out[3]= {0.0}
```

## Distributions

`NormalDistribution[mu, sigma]` and `UniformDistribution[{lo, hi}]` are distribution
objects — ordinary expressions whose head names the family and whose arguments are its
parameters. `NormalDistribution[]` and `UniformDistribution[]` give the standard cases.
Attributes: `Protected`.

**Distribution objects print in full, unlike fitted models**, and that contrast is
deliberate. A distribution is *specified* by its parameters — the user wrote them — so
they are the information. A fitted model's parameters are derived and are an
implementation detail, so it elides them. Same mechanism, opposite convention on
visibility.

## RandomVariate

- `RandomVariate[dist]` — one draw
- `RandomVariate[dist, n]` — a list of `n`

**Features**:
- **Draws come from the same stream as `RandomReal`**, so `SeedRandom` makes them
  reproducible. A sampler with its own generator would silently ignore `SeedRandom`
  while `RandomReal` honoured it — reproducibility half-working is worse than not
  working.
- Normal deviates use Box–Muller in its polar form, which needs no `sin`/`cos`. There
  was no Gaussian deviate anywhere in the tree before this.
- **A non-positive standard deviation, or an inverted range, returns unevaluated** —
  not `NaN`, which would propagate silently through a whole sample and surface much
  later as a strange plot. `RandomVariate[dist, 0]` is a valid empty request.

## PDF

- `PDF[dist, x]`, threading over a list of `x`.

Verified against the closed form computed symbolically: `PDF[NormalDistribution[], 0]`
equals `1/Sqrt[2 Pi]`, and the general case matches the longhand
`Exp[-((x-mu)/sigma)^2/2]/(sigma Sqrt[2 Pi])`. A uniform density is flat inside its
support and exactly zero outside.

```mathematica
In[1]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]
Out[1]= {0.981398, -0.56572, 1.34033, 0.402313, -0.964221}

In[2]:= SeedRandom[1]; s = RandomVariate[NormalDistribution[5., 2.], 20000];
        {Mean[s], StandardDeviation[s]}
Out[2]= {5.00479, 1.99627}

In[3]:= PDF[NormalDistribution[], {-1., 0., 1.}]
Out[3]= {0.241971, 0.398942, 0.241971}
```

## LearnDistribution

Fits a distribution to data and returns a `LearnedDistribution`, usable with `PDF`.
Attributes: `Protected`.

- `LearnDistribution[data]`
- `LearnDistribution[data, Method -> "Multinormal"]`
- `LearnDistribution[data, Method -> "GaussianMixture"]`

**Features**:
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

```
cov_mixture = ((n-1)/n) · cov_multinormal + floor
1.41        = (40/41)   · 1.435           + 0.01
```

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

```mathematica
In[1]:= m = LearnDistribution[{1., 2., 3., 4., 5., 6.}]
Out[1]= LearnedDistribution["Multinormal", <>]

In[2]:= {PDF[m, {3.5}], PDF[NormalDistribution[Mean[{1.,2.,3.,4.,5.,6.}],
                            StandardDeviation[{1.,2.,3.,4.,5.,6.}]], 3.5]}
Out[2]= {0.213244, 0.213244}

In[3]:= PDF[LearnDistribution[{{1.,2.},{2.,3.},{3.,5.},{4.,4.},{5.,7.},{6.,8.}}],
            {{3.5, 4.8}, {50., 50.}}]
Out[3]= {0.113186, 3.6193e-169}
```

## Clustering

`FindClusters` and its ten methods are documented under
[`lists-and-iteration.md`](lists-and-iteration.md), where the rest of the list
operations live. All ten cluster n-dimensional points; only string input is declined
above one dimension, having no coordinates.
