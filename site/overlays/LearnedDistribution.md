### Worked examples

```mathematica
In[1]:= d = LearnDistribution[{1.0, 2.0, 3.0, 4.0, 5.0}]  (* fit a distribution: a Multinormal here *)
```

```mathematica
In[1]:= Head[d]  (* the fitted object is a LearnedDistribution *)
```

```mathematica
In[1]:= PDF[d, {3.0}]  (* usable with PDF; the point is a length-1 vector for a 1-D fit *)
```

### Notes

`LearnedDistribution[method, parameters, dimension, extra]` is what
`LearnDistribution` (and `SmoothKernelDistribution`) return. It prints elided —
`LearnedDistribution["Multinormal", <>]` — because its parameters are derived, not
user-supplied (contrast `NormalDistribution`, which prints in full); `FullForm`
reveals them. Pass it to `PDF[dist, x]`, where `x` is a length-`dimension` vector —
so a 1-D fit is queried at `{v}`, not a bare `v`. Densities are computed in log
space so tail points do not underflow.
