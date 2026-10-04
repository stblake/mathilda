# ClassifierFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClassifierFunction[method, parameters, featureCount, k] is the fitted classifier Classify returns. Apply it to a feature vector to get a class, or to a feature vector and "Probabilities" to get one rule per class.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Classify returns a ClassifierFunction

```mathematica
In[1]:= Head[Classify[{1 -> "low", 2 -> "low", 9 -> "high", 10 -> "high"}]]
Out[1]= ClassifierFunction
```

Apply it to a new point

```mathematica
In[2]:= c = Classify[{1 -> "low", 2 -> "low", 9 -> "high", 10 -> "high"}]; c[1.5]
Out[2]= "low"
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/ml/classify.c`](https://github.com/stblake/mathilda/blob/main/src/ml/classify.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)

## Notes & additional examples

### Notes

`ClassifierFunction[method, parameters, featureCount, k]` is the fitted classifier that
`Classify` returns — a reusable object, not a user-written head. Apply it to a feature
vector to get the predicted class, or to a vector and `"Probabilities"` to get one
`class -> p` rule per class; it also answers `"Method"` and `"FeatureCount"`.

The four stored parts are the method string (e.g. `"NearestNeighbors"`), the fitted
parameters (for nearest-neighbours, the class-label list followed by the labelled
training rows), the input feature count, and `k`, the number of classes. `Classify`
declines a single-class training set — that is not a classification problem — so a
`ClassifierFunction` always has `k >= 2`. It is `Protected`.
