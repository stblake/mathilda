---
source: src/ml/classify.c
---
**Definition.** `ClassifierFunction[method, parameters, featureCount, k]` is the fitted
classifier object that `Classify` returns. It is an **inert head** (no builtin of its
own, just a docstring and the `Protected` attribute); it is produced by `Classify` and
consumed by the prediction machinery in `src/ml/predict.c` when you apply it.

**Representation.** A four-argument `EXPR_FUNCTION`: a `method` string (e.g.
`"NearestNeighbors"`, `"LogisticRegression"`, `"NaiveBayes"`), a `parameters` block
holding the fitted model (for nearest-neighbours this is the class-label list followed
by the labelled training rows), the integer `featureCount`, and `k`, the number of
classes. Because the object can also be typed by hand, `predict.c` treats its contents
as untrusted and validates them before use.

**Usage & limits.** `Protected`. Apply it to a feature vector to get the predicted
class; apply it to a feature vector and `"Probabilities"` to get one `class -> p` rule
per class. It also answers property queries such as `"Method"` and `"FeatureCount"`.
`Classify` declines a single-class training set (not a classification problem), so a
`ClassifierFunction` always has `k >= 2`. It is reusable on inputs it was not trained
on — that is the whole point of returning a function rather than labels.
