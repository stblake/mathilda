---
references:
  - "L. Breiman, J. H. Friedman, R. A. Olshen and C. J. Stone, *Classification and Regression Trees* (Wadsworth, 1984)."
  - "L. Breiman, *Random Forests*, Machine Learning **45** (2001) 5-32."
  - "P. McCullagh and J. A. Nelder, *Generalized Linear Models*, 2nd ed. (Chapman & Hall, 1989), §4.4 (iteratively reweighted least squares)."
  - "T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009)."
source: src/ml/classify.c
---
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
