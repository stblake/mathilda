### Worked examples

```mathematica
In[1]:= Head[Classify[{1 -> "low", 2 -> "low", 9 -> "high", 10 -> "high"}]]  (* Classify returns a ClassifierFunction *)
```

```mathematica
In[1]:= c = Classify[{1 -> "low", 2 -> "low", 9 -> "high", 10 -> "high"}]; c[1.5]  (* apply it to a new point *)
```

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
