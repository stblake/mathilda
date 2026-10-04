# SmoothKernelDistribution

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SmoothKernelDistribution[data] gives a kernel density estimate as a LearnedDistribution, usable with PDF. The kernel is a product Gaussian with a per-dimension bandwidth from the multivariate normal-reference rule, which in one dimension is Silverman's 1.06 sigma n^(-1/5). SmoothKernelDistribution[data, h] sets the bandwidth explicitly, as one number or one per dimension. Being a normal-reference rule the default oversmooths strongly multimodal data -- a known property of the rule, and the reason the explicit form exists. A constant column has no scale and returns unevaluated.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= k = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}]]
Out[1]= LearnedDistribution["SmoothKernel", <>]

In[2]:= PDF[k, {2.0}]
Out[2]= 0.243738
```

### Applications (4)

A KDE over 41 points on 0..4

```mathematica
In[3]:= k = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}]]
Out[3]= LearnedDistribution["SmoothKernel", <>]
```

Evaluate the estimated density

```mathematica
In[4]:= PDF[k, {2.0}]
Out[4]= 0.243738
```

Set the bandwidth explicitly

```mathematica
In[5]:= k2 = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}], 0.5]
Out[5]= LearnedDistribution["SmoothKernel", <>]
```

A different bandwidth gives a different estimate

```mathematica
In[6]:= PDF[k2, {2.0}]
Out[6]= 0.243893
```

## Options & behaviour

**An exact identity is used as the verification.** A KDE is the empirical distribution
*convolved* with the kernel, so its mean is exactly the sample mean and its variance is
exactly the ML sample variance plus `h²`. That one identity pins the bandwidth rule, the
kernel's own variance, and the normalisation simultaneously — it holds to `4e-15`. On 4000
draws from a standard normal the estimate also tracks the true density to within about
0.017, asserted as a *bound* rather than an equality because a KDE genuinely is an
approximation (error `O(n^(-2/5))` for this rule).

## Implementation notes

**Algorithm.** `builtin_smooth_kernel` builds a kernel density estimate in which the
**sample is the model** — nothing is fitted except the bandwidth, so it costs nothing to
build and everything to evaluate, the same trade a nearest-neighbour predictor makes.
The default per-dimension bandwidth is the multivariate **normal-reference** rule
`h_a = σ_a · (4 / ((dim + 2) n))^(1/(dim + 4))`, where `σ_a` comes from `ml_column_sd`
(the `n − 1` divisor, matching `Variance`); in one dimension this is exactly Silverman's
`1.06 σ n^(−1/5)`, the constant being `(4/3)^(1/5) = 1.0592`. A second argument sets the
bandwidth explicitly, as one number for every dimension or one per dimension.

The estimate is stored as `LearnedDistribution["SmoothKernel", {bandwidths, sample
rows…}, dim, n]` and evaluated by `ml_kde_pdf` (in `src/ml/dist.c`): the mean of
**product-Gaussian** kernels — a diagonal kernel with a per-dimension bandwidth, since a
full-covariance kernel would need a bandwidth *matrix* that is far harder to estimate
from the same sample it smooths — summed by log-sum-exp so a point several bandwidths
from every sample does not read as a flat zero.

**Data structures.** The payload `List` carries the bandwidth vector as its first row
and the raw sample as the rest; a transient `MlDist` borrows into it on evaluation.

**Complexity / limits.** `O(n · dim)` to build, `O(n · dim)` per density query. Being a
normal-reference rule the default **oversmooths strongly multimodal data** — a known
property of the rule, and the reason the explicit-bandwidth form exists. A constant
column has no scale, so the rule gives a zero bandwidth and the call declines rather
than dividing by zero.

- **The sample *is* the model.** Nothing is fitted except the bandwidth, which is why a
  KDE costs nothing to build and everything to evaluate — the same trade as a
  nearest-neighbour predictor.
- The kernel is a **product Gaussian** with a per-dimension bandwidth. A
  full-covariance kernel would need a bandwidth *matrix*, and estimating one from the
  same sample it smooths is a substantially harder problem than the default rule solves.
- Default bandwidth is the multivariate **normal-reference** rule,
  `h_a = sigma_a (4/((dim+2) n))^(1/(dim+4))`, which in one dimension is exactly
  Silverman's `1.06 sigma n^(-1/5)` — the constant being `(4/3)^(1/5) = 1.059224`.
- Being a *normal-reference* rule, the default **oversmooths strongly multimodal data**.
  That is a known property of the rule rather than a defect here, and it is why the
  explicit-bandwidth form exists.
- A **constant column** has no scale, so the rule gives a zero bandwidth and there is no
  density: it returns unevaluated rather than dividing by zero.

**Attributes:** `Protected`.

## References

**See also:** [LearnedDistribution](../../other-advanced/LearnedDistribution/), [PDF](../../machine-learning/PDF/)

- B. W. Silverman, *Density Estimation for Statistics and Data Analysis* (Chapman & Hall, 1986), §3.4 (the normal-reference rule).
- D. W. Scott, *Multivariate Density Estimation* (Wiley, 1992).
- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

The sample *is* the model: nothing is fitted except the bandwidth, so a KDE costs
nothing to build and everything to evaluate — the same trade a nearest-neighbour
predictor makes. The kernel is a product Gaussian with a per-dimension bandwidth.

The default bandwidth is the multivariate normal-reference rule, which in one dimension
is exactly Silverman's `1.06 σ n^(−1/5)`. Being a normal-reference rule it oversmooths
strongly multimodal data — a known property of the rule rather than a defect — which is
why the explicit-bandwidth form shown above exists. A constant column has no scale, so
the rule gives a zero bandwidth and the call returns unevaluated rather than dividing by
zero.
