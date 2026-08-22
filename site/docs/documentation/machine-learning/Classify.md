# Classify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Classify[data] trains a classifier and returns a ClassifierFunction. Data is a list of rules {features -> class, ...}; a class may be any expression -- a string, a symbol, a number -- and the distinct classes are numbered by first appearance. Method -> "NearestNeighbors" is the only method implemented and is the default, with NeighborsNumber defaulting to 1: a classifier votes rather than averages, so at k = 1 it reproduces its training labels exactly. Apply the result to a feature vector for a class, or with "Probabilities" for the vote shares. It also answers "Classes", "Method", "FeatureCount" and "NeighborCount". Method -> "NaiveBayes" fits a Gaussian per class with a diagonal covariance; Method -> "LogisticRegression" fits a logistic model by iteratively reweighted least squares with a small ridge on the non-intercept coefficients -- the ridge is load-bearing, because on linearly separable data the unpenalised likelihood is unbounded and the coefficients would diverge. Two classes give a single fit; more than two are fitted one-vs-rest, one binary model per class, and the class is the arg-max of the fitted probabilities. Those probabilities are normalised to sum to 1, which is a convention rather than a likelihood -- being monotone it cannot change the arg-max, so the class is the better-founded of the two answers. A single class declines: it is not a classification problem.`**

## Examples (3)

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

- Source: [`src/ml/classify.c`](https://github.com/stblake/mathilda/blob/main/src/ml/classify.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
