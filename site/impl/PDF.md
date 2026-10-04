---
source: src/ml/dist.c
---
**Algorithm.** `builtin_pdf` tries the nominal path first: `ct_pdf` recognises a
`LearnedDistribution["ContingencyTable", …]` and looks the outcome up structurally
(`expr_eq`), returning an exact `0` for an outcome never observed — this runs before the
numeric reader because the argument is an outcome, not a number. Otherwise `ml_read_dist`
classifies the distribution and `ml_pdf_at` (or a point kernel) evaluates it:

- **Normal** — the closed form `exp(−½ z²) / (σ √(2π))` with `z = (x − μ)/σ`.
- **Uniform** — `1/(b − a)` on the closed interval `[a, b]`, zero strictly outside.
- **Multinormal / GaussianMixture / SmoothKernel** — evaluated in **log space** and
  exponentiated once. A `dim`-factor Gaussian density underflows to zero for a point a
  few standard deviations out, so `ml_multinormal_pdf` works through a Mahalanobis
  distance against the stored Cholesky factor plus its log-determinant, while
  `ml_mixture_pdf` and `ml_kde_pdf` combine their terms by log-sum-exp.

The argument reading differs by kind, and it must: for a **scalar** distribution a
`List` of `x` threads to a list of densities (what a caller plotting a density wants),
whereas for a **point** distribution the argument is itself a vector, so a list is *one*
observation and a *matrix* threads to one density per row.

**Data structures.** A transient `MlDist` struct holding borrowed pointers into one
owned buffer decoded from the distribution object; Cholesky factors and log-determinants
are recomputed on read rather than stored twice.

**Complexity / limits.** Scalar `O(1)` per point; Multinormal `O(dim²)`; mixture
`O(k · dim²)`; KDE `O(n · dim)` per query. Verified against symbolically-computed closed
forms — `PDF[NormalDistribution[], 0]` equals `1/Sqrt[2 Pi]` — and, for the point
kernels, against an independent Cholesky-based path that also integrates to `1`.
