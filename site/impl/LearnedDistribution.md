---
source: src/ml/dist.c
---
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
