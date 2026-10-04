---
references:
  - "A. P. Dempster, N. M. Laird and D. B. Rubin, *Maximum likelihood from incomplete data via the EM algorithm*, J. Roy. Statist. Soc. B **39** (1977) 1-38."
  - "G. Schwarz, *Estimating the dimension of a model*, Ann. Statist. **6** (1978) 461-464 (BIC)."
  - "C. M. Bishop, *Pattern Recognition and Machine Learning* (Springer, 2006), §9.2 (Gaussian mixtures)."
source: src/ml/dist.c
---
**Algorithm.** `builtin_learn_distribution` fits one of three families and returns a
`LearnedDistribution`, which prints elided (its parameters are derived) — the deliberate
opposite of a *specified* distribution such as `NormalDistribution[μ, σ]`.

- **`"Multinormal"`** (default). `ml_column_mean` and a sample covariance with the
  `n − 1` divisor (matching `Variance`, so a one-variable fit agrees with
  `StandardDeviation²`), factored by `ml_chol`. A singular covariance — collinear
  columns, or fewer observations than dimensions — declines, because no density exists.

- **`"GaussianMixture"`**. BIC model selection over `k = 1 … kmax`, each `k` fitted by
  `ml_gmm_fit` (EM in `src/ml/gmm.c`): deterministic farthest-first initialisation, a
  **log-space E-step** via log-sum-exp (a `dim`-factor density underflows for an
  outlying point, and a linear-space zero would hand it uniform responsibilities —
  silently refusing to distinguish exactly the informative points), and an M-step whose
  covariance floor is added as a **ridge** `floor · I` on the diagonal rather than a
  clamp. The variance floor itself (`ml_nn_floor`) is the **squared median
  nearest-neighbour distance between *distinct* points** — load-bearing, because a
  mixture likelihood is unbounded above and a merely-"small" floor lets BIC buy
  arbitrarily many near-singular spikes; counting duplicates' zero distances once drove
  the floor to `1e-300` and the fit to nine spurious components. A scale guard holds it
  to at least `1e-4` of the mean per-coordinate variance. `kmax` is bounded by
  `n/(dim+1)`, the distinct-point count, and 10; `BIC = paramCount · ln n − 2 · loglik`.

- **`"ContingencyTable"`** — nominal, not numeric. It runs *before* the numeric reader,
  builds a label vocabulary (`ml_labels_build`, `expr_eq`, first-appearance order) over
  outcomes that may be any expressions or equal-length lists of them, and stores
  frequency `count/n` per outcome. Probabilities are **empirical with no smoothing** —
  an unseen outcome is exactly `0`, since smoothing would require knowing the size of an
  unknowable outcome space. Ragged outcomes decline.

**Data structures.** An `MlGmm` (weights, means, covariances, Cholesky factors,
log-determinants) for the mixture; the stored `LearnedDistribution[method, payload, dim,
extra]` payload is a plain `List`.

**Complexity / limits.** Multinormal `O(n · dim² + dim³)`; the mixture search is
`O(kmax · iterations · n · k · dim²)` with `ml_nn_floor` an `O(n² · dim)` preprocess.
A flat list is accepted as `n` observations of one variable (a univariate normal),
unlike `PrincipalComponents`.
