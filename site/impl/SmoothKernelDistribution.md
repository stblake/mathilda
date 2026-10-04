---
references:
  - "B. W. Silverman, *Density Estimation for Statistics and Data Analysis* (Chapman & Hall, 1986), §3.4 (the normal-reference rule)."
  - "D. W. Scott, *Multivariate Density Estimation* (Wiley, 1992)."
source: src/ml/dist.c
---
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
