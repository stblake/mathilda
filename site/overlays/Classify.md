### Worked examples

```mathematica
In[1]:= c = Classify[{{0., 0.} -> "red", {1., 0.} -> "red", {0., 1.} -> "red", {10., 10.} -> "blue", {11., 10.} -> "blue", {10., 11.} -> "blue"}]  (* nearest-neighbour classifier, the default *)
```

```mathematica
In[1]:= {c[{0.5, 0.5}], c[{10.5, 10.5}]}  (* classify two query points *)
```

```mathematica
In[1]:= c[{0.5, 0.5}, "Probabilities"]  (* vote shares, one rule per class *)
```

```mathematica
In[1]:= c["Classes"]  (* classes in first-appearance order *)
```

```mathematica
In[1]:= nb = Classify[{1. -> "lo", 2. -> "lo", 8. -> "hi", 9. -> "hi"}, Method -> "NaiveBayes"]  (* a Gaussian naive Bayes classifier *)
```

```mathematica
In[1]:= {nb[1.5], nb[8.5]}  (* classify either side of the gap *)
```

```mathematica
In[1]:= SeedRandom[42]; f = Classify[{{0., 0.} -> "a", {1., 0.} -> "a", {0., 1.} -> "a", {5., 5.} -> "b", {6., 5.} -> "b", {5., 6.} -> "b", {10., 0.} -> "c", {11., 0.} -> "c", {10., 1.} -> "c"}, Method -> "RandomForest"]  (* fifty bagged trees, reproducible under SeedRandom *)
```

```mathematica
In[1]:= {f[{0.5, 0.5}], f[{5.5, 5.5}], f[{10.5, 0.5}]}  (* a three-class prediction *)
```

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
