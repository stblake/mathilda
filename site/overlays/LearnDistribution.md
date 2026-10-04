### Worked examples

```mathematica
In[1]:= m = LearnDistribution[{1., 2., 3., 4., 5., 6.}]  (* a univariate normal from a flat list *)
```

```mathematica
In[1]:= PDF[m, {3.5}]  (* evaluate the fitted density at one point *)
```

```mathematica
In[1]:= d = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]  (* nominal outcomes, not numbers *)
```

```mathematica
In[1]:= {PDF[d, "r"], PDF[d, "b"], PDF[d, "g"]}  (* an outcome never observed has probability exactly 0 *)
```

```mathematica
In[1]:= Last[LearnDistribution[{1.4, 1.4, 1.3, 1.5, 1.4, 1.7, 1.4, 1.5, 1.4, 1.5, 1.5, 1.6, 1.4, 1.1, 1.2, 4.7, 4.5, 4.9, 4., 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4., 4.7, 3.6}, Method -> "GaussianMixture"]]  (* BIC chooses two components for bimodal data *)
```

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
