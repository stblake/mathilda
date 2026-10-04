# LearnedDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LearnedDistribution[method, parameters, dimension, extra] is the fitted distribution LearnDistribution returns. Use it with PDF. Unlike a SPECIFIED distribution such as NormalDistribution[mu, sigma], it prints elided, because its parameters are derived rather than user-supplied.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Fit a distribution: a Multinormal here

```mathematica
In[1]:= d = LearnDistribution[{1.0, 2.0, 3.0, 4.0, 5.0}]
Out[1]= LearnedDistribution["Multinormal", <>]
```

The fitted object is a LearnedDistribution

```mathematica
In[2]:= Head[d]
Out[2]= LearnedDistribution
```

Usable with PDF; the point is a length-1 vector for a 1-D fit

```mathematica
In[3]:= PDF[d, {3.0}]
Out[3]= 0.252313
```

## Implementation notes

**Definition.** `LearnedDistribution[method, parameters, dimension, extra]` is the
fitted distribution object that `LearnDistribution` and `SmoothKernelDistribution`
return. Use it with `PDF`. Unlike a *specified* distribution such as
`NormalDistribution[mu, sigma]`, it prints **elided** —
`LearnedDistribution["method", <>]` — because its parameters are derived rather than
user-supplied. It is `Protected`, has no builtin of its own, and its docstring is in
`dist.c`.

**Representation.** A four-slot inert head: a method-tag string (`"Multinormal"`,
`"SmoothKernel"`, `"GaussianMixture"`), the fitted `parameters` (mean and covariance
rows; or kernel bandwidths and sample rows; or mixture weights and components — a
positional layout that differs per method), the data `dimension`, and an `extra`
count (number of kernel samples or mixture components). `ml_read_learned`
(`src/ml/dist.c`) reads whichever shape the method tag names into an `MlDist` record;
`PDF` then evaluates the density in **log space** — log-sum-exp across mixture
components, a Cholesky Mahalanobis form for the multinormal — so a point several
standard deviations out does not underflow to a flat zero.

**Usage & limits.** Produced by fitting, not written by hand; pass it to `PDF[dist,
x]`, where `x` is a length-`dimension` vector — so for a 1-D fit the point is `{v}`,
not a bare scalar. The elided print hides the derived parameters by design;
`FullForm` still shows them.

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`LearnedDistribution[method, parameters, dimension, extra]` is what
`LearnDistribution` (and `SmoothKernelDistribution`) return. It prints elided —
`LearnedDistribution["Multinormal", <>]` — because its parameters are derived, not
user-supplied (contrast `NormalDistribution`, which prints in full); `FullForm`
reveals them. Pass it to `PDF[dist, x]`, where `x` is a length-`dimension` vector —
so a 1-D fit is queried at `{v}`, not a bare `v`. Densities are computed in log
space so tail points do not underflow.
