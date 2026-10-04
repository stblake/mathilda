# Classify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Classify[data] trains a classifier and returns a ClassifierFunction. Data is a list of rules {features -> class, ...}; a class may be any expression -- a string, a symbol, a number -- and the distinct classes are numbered by first appearance. Method -> "NearestNeighbors" is the only method implemented and is the default, with NeighborsNumber defaulting to 1: a classifier votes rather than averages, so at k = 1 it reproduces its training labels exactly. Apply the result to a feature vector for a class, or with "Probabilities" for the vote shares. It also answers "Classes", "Method", "FeatureCount" and "NeighborCount". Method -> "NaiveBayes" fits a Gaussian per class with a diagonal covariance; Method -> "LogisticRegression" fits a logistic model by iteratively reweighted least squares with a small ridge on the non-intercept coefficients -- the ridge is load-bearing, because on linearly separable data the unpenalised likelihood is unbounded and the coefficients would diverge. Two classes give a single fit; more than two are fitted one-vs-rest, one binary model per class, and the class is the arg-max of the fitted probabilities. Those probabilities are normalised to sum to 1, which is a convention rather than a likelihood -- being monotone it cannot change the arg-max, so the class is the better-founded of the two answers. A single class declines: it is not a classification problem.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= c = Classify[{{0.,0.} -> "red", {1.,0.} -> "red", {0.,1.} -> "red", {10.,10.} -> "blue", {11.,10.} -> "blue", {10.,11.} -> "blue"}]
Out[1]= ClassifierFunction["NearestNeighbors", <>]

In[2]:= {c[{0.5, 0.5}], c[{10.5, 10.5}], c["Classes"]}
Out[2]= {"red", "blue", {"red", "blue"}}

In[3]:= c[{0.5, 0.5}, "Probabilities"]
Out[3]= {"red" -> 1.0, "blue" -> 0.0}
```

### Applications (8)

Nearest-neighbour classifier, the default

```mathematica
In[4]:= c = Classify[{{0., 0.} -> "red", {1., 0.} -> "red", {0., 1.} -> "red", {10., 10.} -> "blue", {11., 10.} -> "blue", {10., 11.} -> "blue"}]
Out[4]= ClassifierFunction["NearestNeighbors", <>]
```

Classify two query points

```mathematica
In[5]:= {c[{0.5, 0.5}], c[{10.5, 10.5}]}
Out[5]= {"red", "blue"}
```

Vote shares, one rule per class

```mathematica
In[6]:= c[{0.5, 0.5}, "Probabilities"]
Out[6]= {"red" -> 1.0, "blue" -> 0.0}
```

Classes in first-appearance order

```mathematica
In[7]:= c["Classes"]
Out[7]= {"red", "blue"}
```

A Gaussian naive Bayes classifier

```mathematica
In[8]:= nb = Classify[{1. -> "lo", 2. -> "lo", 8. -> "hi", 9. -> "hi"}, Method -> "NaiveBayes"]
Out[8]= ClassifierFunction["NaiveBayes", <>]
```

Classify either side of the gap

```mathematica
In[9]:= {nb[1.5], nb[8.5]}
Out[9]= {"lo", "hi"}
```

Fifty bagged trees, reproducible under SeedRandom

```mathematica
In[10]:= SeedRandom[42]; f = Classify[{{0., 0.} -> "a", {1., 0.} -> "a", {0., 1.} -> "a", {5., 5.} -> "b", {6., 5.} -> "b", {5., 6.} -> "b", {10., 0.} -> "c", {11., 0.} -> "c", {10., 1.} -> "c"}, Method -> "RandomForest"]
Out[10]= ClassifierFunction["RandomForest", <>]
```

A three-class prediction

```mathematica
In[11]:= {f[{0.5, 0.5}], f[{5.5, 5.5}], f[{10.5, 0.5}]}
Out[11]= {"a", "b", "c"}
```

## Options & behaviour

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

## Algorithm

classify.c -- Classify and ClassifierFunction.

A ClassifierFunction is the FOURTH head on the model representation designed in src/ml/predict.h, after PredictorFunction, DimensionReducerFunction and LearnedDistribution. It needed no change to that design: the payload shape varies by method, which is exactly what a positional method-tagged representation is for.

What IS new is the label vocabulary (src/ml/encode.h). Every earlier family took numeric responses; a class is an arbitrary expression, so the vocabulary is the bridge between "the user's classes" and "indices an algorithm can count with". It lives in its own module because a ContingencyTable and a categorical FEATURE encoder will both need it.

## Implementation notes

**Algorithm.** `builtin_classify` reads `{features -> class, …}` only (`ml_read_labelled`
— a matrix is refused, since a class need not be numeric), builds a **label vocabulary**
(`src/ml/encode.c`: distinct classes by `expr_eq`, first-appearance order) and maps each
response to a class index. Five methods share the positional `ClassifierFunction[method,
payload, dim, k]` representation:

- **`"NearestNeighbors"`** (default, `k = 1` — a classifier votes, so `k = 1` reproduces
  the training labels exactly; clamped to `n`). The payload is the vocabulary plus one
  row per example, its features followed by its class *index*. Application takes the
  majority vote of the `k` nearest by squared Euclidean distance (partial insertion
  selection), ties to the lowest index.
- **`"NaiveBayes"`** — per class, a mean and a per-feature variance (ML `n` divisor) plus
  the prior; classifies by the largest log posterior, `log prior + Σ` Gaussian log
  densities (log space, then one softmax for `"Probabilities"`). The variance floor is a
  fraction (`1e-6`) of the feature's *overall* variance, so it is scale-invariant — a
  zero-variance feature in a class would otherwise give infinite density.
- **`"LogisticRegression"`** — `logit_irls` runs Newton's method on the log-likelihood
  (iteratively reweighted least squares, 100-iteration cap) with a small ridge `1e-6` on
  the non-intercept coefficients, left unpenalised on the intercept. The ridge is
  load-bearing: on linearly separable data the unpenalised likelihood is unbounded and
  plain Newton diverges. Two classes give a single fit; more than two are **one-vs-rest**
  (K binary fits, arg-max of the K probabilities), chosen over a softmax because a
  softmax's parameters are identified only up to a per-feature constant and so could not
  be pinned by a test. The K sigmoids are normalised to sum to 1 on read — a monotone
  presentation convention that cannot move the arg-max.
- **`"DecisionTree"`** — a CART tree (`src/ml/tree.c`): Gini impurity, thresholds at
  midpoints between consecutive distinct values, grown until every leaf is pure or
  unsplittable (depth 32, min-split 2), with a fully deterministic tie-break (lower
  feature, lower threshold, point index). The payload is the vocabulary plus two
  one-row-per-node matrices — the split `(feature, threshold, left, right)` and the
  class counts at *every* node, so class and `"Probabilities"` come from one array.
- **`"RandomForest"`** — 50 CART trees, each on a bootstrap resample, each node choosing
  among `√dim` features; per-tree distributions are normalised then averaged, arg-max
  wins. Both the bootstrap draw and the per-node feature sample come from
  `random_uniform_01`, so a forest is **reproducible under `SeedRandom`**.

**Data structures.** An `MlLabels` vocabulary; an `MlTree` of parallel arrays
(`feature`, `thresh`, `left`, `right`, `dist`) during fitting, serialised to `List`s for
the payload. Application (`ml_classifier_apply`) walks the payload directly, bounding its
step count by the node count because a hand-typed classifier's `left`/`right` are
untrusted.

**Complexity / limits.** k-NN apply `O(n · dim)`; NaiveBayes `O(n · dim)`; logistic
`O(iterations · K · (n · p² + p³))`; a tree node is `O(dim · n log n)` from the per-feature
sort, the forest ×50. Properties `"Classes"`, `"Method"`, `"FeatureCount"` and (k-NN
only) `"NeighborCount"`; `[x, "Probabilities"]` gives one rule per class. A single class
declines; `"NeighborsNumber"` on any non-k-NN method is refused, not ignored.

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

**Attributes:** `Protected`.

## References

**See also:** [ClassifierFunction](../../other-advanced/ClassifierFunction/), [FindClusters](../../lists-and-iteration/FindClusters/), [Predict](../../machine-learning/Predict/)

- L. Breiman, J. H. Friedman, R. A. Olshen and C. J. Stone, *Classification and Regression Trees* (Wadsworth, 1984).
- L. Breiman, *Random Forests*, Machine Learning **45** (2001) 5-32.
- P. McCullagh and J. A. Nelder, *Generalized Linear Models*, 2nd ed. (Chapman & Hall, 1989), §4.4 (iteratively reweighted least squares).
- T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009).
- Source: [`src/ml/classify.c`](https://github.com/stblake/mathilda/blob/main/src/ml/classify.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)

## Notes & additional examples

### Notes

Data is a list of rules `{features -> class, …}`; a matrix with the class in its last
column is *not* accepted, because a class need not be numeric. A class may be any
expression — string, symbol, number — compared structurally, so `"a"` and `a` are two
classes. The distinct classes form a vocabulary numbered by first appearance, which is
deterministic; the prediction does not depend on the order.

`k` defaults to **1**, unlike the `Predict` nearest-neighbour *regressor* where it
defaults to 3: a classifier votes rather than averages, so at `k = 1` it reproduces its
training labels exactly. `classifier[x, "Probabilities"]` gives one rule per class with
the vote shares, which sum to 1 by construction.

Six methods are available — `"NearestNeighbors"` (default), `"NaiveBayes"`,
`"LogisticRegression"`, `"DecisionTree"`, `"RandomForest"` — and the two with randomness,
the forest and (through feature sampling) its trees, draw from the same stream as
`RandomReal`, so `SeedRandom` before the fit makes a forest reproducible. The classifier
also answers `"Classes"`, `"Method"`, `"FeatureCount"` and `"NeighborCount"`.
