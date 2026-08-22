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
- `LearnDistribution[data, Method -> "ContingencyTable"]` — nominal, not numeric

### `"ContingencyTable"` — nominal outcomes

The one method here that does not take numbers. It stores a probability per distinct
outcome, which in one dimension is a categorical distribution:

```
In[1]:= d = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]
Out[1]= LearnedDistribution["ContingencyTable", <>]

In[2]:= {PDF[d, "r"], PDF[d, "b"], PDF[d, "g"]}
Out[2]= {0.75, 0.25, 0.}
```

Outcomes may be **any expressions** — strings, symbols, numbers, or equal-length lists of
them — compared structurally, so `"a"` and `a` are two outcomes rather than one. Order is
first appearance, the same contract every other vocabulary in `src/ml` uses. Several nominal
variables are a list per observation, which is the cross-tabulation case:
`LearnDistribution[{{"a","x"}, {"a","y"}, {"b","x"}, {"a","x"}}, …]` gives `{a,x}` a
probability of `0.5` and the other two `0.25` each.

This branch runs **before** the numeric reader, because that reader reads numbers and the
whole point of this method is outcomes that are not numbers — pushing them through it would
decline exactly the data the method exists for. The label vocabulary in `src/ml/encode.h`
does the work unchanged; this is its third consumer.

**Probabilities are empirical frequencies with no smoothing, and the consequence is stated
rather than hidden:** an outcome never observed has probability **exactly 0**. Add-one
smoothing would need a claim about the size of the outcome space, and for arbitrary
expressions there is no way to know how many nominal values were possible but unseen. A hard
zero for "never observed" is the honest answer.

**Ragged outcomes decline** — differing lengths, or a mix of list and non-list observations,
would be a different table per row, which is not a distribution.

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

## SmoothKernelDistribution

A kernel density estimate, returned as a `LearnedDistribution` usable with `PDF`.
Attributes: `Protected`.

- `SmoothKernelDistribution[data]`
- `SmoothKernelDistribution[data, h]` — bandwidth as one number, or one per dimension

**Features**:
- **The sample *is* the model.** Nothing is fitted except the bandwidth, which is why a
  KDE costs nothing to build and everything to evaluate — the same trade as a
  nearest-neighbour predictor.
- The kernel is a **product Gaussian** with a per-dimension bandwidth. A
  full-covariance kernel would need a bandwidth *matrix*, and estimating one from the
  same sample it smooths is a substantially harder problem than the default rule solves.
- Default bandwidth is the multivariate **normal-reference** rule,
  `h_a = sigma_a (4/((dim+2) n))^(1/(dim+4))`, which in one dimension is exactly
  Silverman's `1.06 sigma n^(-1/5)` — the constant being `(4/3)^(1/5) = 1.059224`.
- Being a *normal-reference* rule, the default **oversmooths strongly multimodal data**.
  That is a known property of the rule rather than a defect here, and it is why the
  explicit-bandwidth form exists.
- A **constant column** has no scale, so the rule gives a zero bandwidth and there is no
  density: it returns unevaluated rather than dividing by zero.

**An exact identity is used as the verification.** A KDE is the empirical distribution
*convolved* with the kernel, so its mean is exactly the sample mean and its variance is
exactly the ML sample variance plus `h²`. That one identity pins the bandwidth rule, the
kernel's own variance, and the normalisation simultaneously — it holds to `4e-15`. On 4000
draws from a standard normal the estimate also tracks the true density to within about
0.017, asserted as a *bound* rather than an equality because a KDE genuinely is an
approximation (error `O(n^(-2/5))` for this rule).

```mathematica
In[1]:= k = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}]]
Out[1]= LearnedDistribution["SmoothKernel", <>]

In[2]:= PDF[k, {2.0}]
Out[2]= 0.243738
```

## Classify

Trains a classifier and returns a `ClassifierFunction`. Attributes: `Protected`.

- `Classify[data]`
- `Classify[data, Method -> "NearestNeighbors"]`
- `Classify[data, Method -> "NaiveBayes"]`
- `Classify[data, Method -> "LogisticRegression"]`
- `Classify[data, Method -> "DecisionTree"]`
- `Classify[data, Method -> "RandomForest"]`
- `Classify[data, Method -> {"NearestNeighbors", "NeighborsNumber" -> k}]`

**A class may be any expression.** This is the substantive addition of this family:
every earlier one took numeric responses, and a class is a string, a symbol, a number, or
whatever the user names it. The distinct classes form a **label vocabulary** stored in the
model, and the numeric part of the model carries class *indices* — so the vocabulary is the
single place a class is named. Comparison is structural, so `"a"` and `a` are two classes
rather than one, the same distinction the pattern matcher makes everywhere else.

**Class order is first appearance in the training data**, not sorted. What matters is that
it is deterministic; first appearance is also stable under adding classes later, matches how
`FindClusters` numbers clusters, and reads naturally beside the data that produced it. The
*prediction* does not depend on it.

**Features**:
- Data must be a list of rules `{features -> class, …}`. A matrix with the class in its last
  column is *not* accepted, and that is not an omission — a class need not be a number, so
  the numeric matrix reader would refuse the whole thing. `Predict` accepts a matrix
  precisely because its response is numeric.
- `k` defaults to **1**, unlike the k-NN *predictor* where it defaults to 3. A regression
  averages, so a little smoothing helps; a classifier votes, and at `k = 1` it reproduces
  its training labels exactly.
- `classifier[x, "Probabilities"]` gives one rule per class with the vote shares, which sum
  to 1 by construction.
- Also answers `"Classes"`, `"Method"`, `"FeatureCount"` and `"NeighborCount"`.
- Ties in the vote go to the lowest class index, i.e. first appearance — deterministic.

**`"NaiveBayes"`** fits, per class, a mean and a per-feature variance plus the class
prior, and classifies by the largest log posterior. "Naive" is the independence
assumption — the joint density is the *product* of per-feature densities, i.e. a diagonal
covariance — which is why it needs no Cholesky and works with far fewer points per class
than a full-covariance Multinormal.

**Its variance floor is load-bearing.** A class whose feature takes one value everywhere
has zero variance there and therefore *infinite* density at that value, which would win
every comparison involving that feature. The floor is a fraction of the feature's
**overall** variance across all classes rather than a fixed epsilon, so it is
scale-invariant: the same feature measured in millimetres and in kilometres gets
proportionate floors, where a fixed epsilon would be enormous for one and negligible for
the other. Per-class variances use the ML (`n`) divisor, which the floor is what makes safe
for a single-member class.

`"NeighborsNumber"` is **refused** on a Bayes classifier rather than ignored, and
`"NeighborCount"` is not one of its properties.

**Verified against a closed form.** With one feature and equal priors the decision is just
"which prior-weighted Gaussian density is larger", which `PDF[NormalDistribution[…]]`
computes by a completely separate path. The two agree at eight points including 4.9 and 5.1
— either side of the boundary at 5.0 — so the test exercises the decision rather than two
obvious regions.

**`"LogisticRegression"`** fits a logistic model by iteratively reweighted least squares
(Newton's method on the log-likelihood).

**Two classes** give a single fit and a single coefficient vector. **More than two** are
fitted **one-vs-rest**: K binary models, class *k* against everything else, and the reported
class is the arg-max of the K fitted probabilities.

One-vs-rest was chosen over a softmax for two reasons, and the second is the deciding one.
It reuses the identical IRLS iteration K times rather than needing a different one. And a
softmax's parameters are identified only up to an additive constant per feature, so the
stored coefficients would not be unique — meaning no test could pin them, which is exactly
the kind of assertion that caught real bugs in the other four families.

What one-vs-rest does **not** give is calibrated probabilities. The K sigmoids come from K
separate fits with nothing tying them together, so they are normalised to sum to 1 on read.
That normalisation is a presentation convention, not a likelihood — but because it is
monotone it cannot move the arg-max, so **the class is the better-founded of the two
answers**, and it is the one the tests pin hardest. If every sigmoid underflows to zero the
normaliser would be zero, and that case declines rather than inventing a uniform answer.

The two shapes are told apart by the payload itself — a flat list of `dim + 1` reals is the
two-class fit, a list of K such lists is one-vs-rest — so no extra tag is stored. A **single**
class declines: it is not a classification problem.

`Classify[{1. -> "a", 5. -> "b", 9. -> "c"}, Method -> "LogisticRegression"]` used to
decline, and a test pinned that refusal so it would be noticed when the gap closed. It now
fits, and the test pins the answer instead.

**A small ridge on the non-intercept coefficients, and it is load-bearing.** On linearly
separable data the unpenalised likelihood is **unbounded**: driving the coefficients to
infinity drives every fitted probability to 0 or 1, so plain Newton diverges and never
converges. The ridge makes the penalised objective strictly concave, so the fit is finite and
unique even when the data *is* separable; the iteration count is capped as a backstop. The
intercept is left unpenalised, which is standard — shrinking it would bias the predicted base
rate. This is the same shape of problem as a mixture's unbounded likelihood, handled the same
honest way rather than left to hang.

**Verified by an exact identity rather than an accuracy figure.** The fitted boundary is where
`intercept + coef·x = 0`, and the probability there must be *exactly* 0.5 — because that is
what the logistic of zero is. That single assertion ties the **fit** and the **application**
together: any error in how coefficients are stored, read back, or recombined shows up in it.
On the test data the boundary lands at exactly 5.0, the midpoint of the gap, with coefficients
`{-33.27, 6.65}` — finite, which is the ridge working.

```mathematica
In[1]:= c = Classify[{{0.,0.} -> "red", {1.,0.} -> "red", {0.,1.} -> "red",
                      {10.,10.} -> "blue", {11.,10.} -> "blue", {10.,11.} -> "blue"}]
Out[1]= ClassifierFunction["NearestNeighbors", <>]

In[2]:= {c[{0.5, 0.5}], c[{10.5, 10.5}], c["Classes"]}
Out[2]= {"red", "blue", {"red", "blue"}}

In[3]:= c[{0.5, 0.5}, "Probabilities"]
Out[3]= {"red" -> 1.0, "blue" -> 0.0}
```

## Clustering

`FindClusters` and its ten methods are documented under
[`lists-and-iteration.md`](lists-and-iteration.md), where the rest of the list
operations live. All ten cluster n-dimensional points; only string input is declined
above one dimension, having no coordinates.

### `"DecisionTree"` — a CART classification tree

The first method in `src/ml` that reuses **nothing** from the other families: no distance, no
density, no linear algebra, no label-index arithmetic beyond counting. A tree needs an
impurity criterion, a recursive splitter, a stopping rule and a node representation, none of
which any earlier family had a use for — which is why it was recorded as deferred for as long
as it was rather than being wedged into an iteration. The kernel lives in `src/ml/tree.c`,
buffer-level, because `RandomForest` is this same fit run over bootstrap resamples.

**Gini rather than entropy**, and not as a coin toss. Both rank splits almost identically in
practice, but Gini is a sum of squares where entropy is a sum of `x log x` — no logarithm per
candidate threshold and no special case at `p = 0`. The splitter evaluates the criterion
O(`dim` × `n`) times per node, so the cheaper one with no domain edge is the better default.

**Thresholds are midpoints between consecutive distinct values**, never the values themselves.
A threshold sitting on a training value makes the `<=` boundary depend on floating-point
equality with it, so an unseen point equal to that value lands by luck; the midpoint puts the
boundary in the gap where nothing lies. Between `1.` and `9.` the split is at exactly `5.`,
and a test pins that.

**It grows until every leaf is pure or unsplittable** (depth 32, min-split 2). That is a
deliberate default: it makes "reproduces every training label" an *exact* property to assert
rather than an accuracy figure to hope for — the role `k = 1` plays for nearest neighbours. It
also **overfits**, which is the honest trade. Pruning needs a validation split or a complexity
parameter, and inventing either silently would be worse than growing the tree that was asked
for.

**Determinism is a requirement here, not a courtesy.** Ties on impurity decrease break by
lower feature index then lower threshold, and the sort comparator falls back to the point
index so `qsort`'s permitted instability cannot decide anything. Without all of that no test
could pin a tree at all — the same reasoning that ruled out a softmax for multi-class logistic
regression, whose parameters are likewise not unique. `Classify[d, …] === Classify[d, …]` is
asserted.

The payload is the vocabulary plus **two matrices with one row per node**: the split
(`feature`, `threshold`, `left`, `right`, with `feature = -1` marking a leaf) and the class
counts. Storing counts at *every* node rather than only at leaves is what lets the class and
`"Probabilities"` come out of the same array, with no separate leaf table to keep in step; an
internal node's histogram is its subtree's distribution, so nothing in the layout is dead
weight.

```
In[1]:= Part[Classify[{1. -> "lo", 9. -> "hi"}, Method -> "DecisionTree"], 2, 2]
Out[1]= {{0, 5.0, 1, 2}, {-1, 0.0, 0, 0}, {-1, 0.0, 0, 0}}
```

**Data the tree cannot separate is handled rather than fatal.** Identical feature rows with
different classes share a leaf holding the majority, and its `"Probabilities"` are that leaf's
real histogram — three points at one coordinate with classes `a`, `b`, `b` give `1/3` and
`2/3`. "No split improves impurity" is a stopping rule, and it is what would otherwise recurse
forever when every feature is constant. A single training point is a valid one-leaf tree.

### `"RandomForest"` — fifty trees, bagged and feature-sampled

Fifty CART trees, each grown on a bootstrap resample of the rows, each node choosing among
only `sqrt(dim)` features. The class is the arg-max of the averaged predicted distribution.

**Both sources of randomness are needed, and the second matters more.** Bootstrap resampling
alone leaves the trees highly correlated: if every tree may pick the single most informative
feature at its root, they all do, and averaging near-identical trees buys nothing. Restricting
each node's candidate set is what decorrelates them. `sqrt(dim)` is the standard choice for
classification.

**A forest is reproducible under `SeedRandom`.** Both the bootstrap draw and the per-node
feature sample come from `random_uniform_01` — the same stream `RandomReal` uses — so
`SeedRandom[42]` before two fits gives two identical forests, and a different seed gives a
different one. Both are asserted. This is load-bearing rather than a nicety: a fit that varied
run to run could not be pinned by any test, which is the same concern that made multi-class
logistic regression one-vs-rest instead of a softmax.

**Per-tree distributions are normalised before averaging, not pooled as raw counts.** Pooling
would let one tree with a large leaf outvote several trees with small ones, so each tree would
carry a weight set by its own depth rather than one vote. Averaging distributions is the
standard choice and the defensible one.

A pleasant consequence: the forest's probabilities are **genuinely mixed** near a boundary,
because the trees disagree there. A single tree on separable data always answers `1/0`, so its
sum-to-one assertion needed contrived contradictory data to mean anything; the forest's does
not.

```
In[1]:= SeedRandom[42]; f = Classify[data, Method -> "RandomForest"];
In[2]:= f[{2.5, 2.5}, "Probabilities"]
Out[2]= {"a" -> 0.641, "b" -> 0.174, "c" -> 0.185}
```

**Fifty trees, fixed.** Wolfram tunes the count against held-out data. Doing that here would
mean carving a validation split out of the user's data without being asked, and guessing
silently is worse than a documented constant.
