---
references:
  - "G. Marsaglia and T. A. Bray, *A convenient method for generating normal variables*, SIAM Review **6** (1964) 260-264 (the polar method)."
  - "D. E. Knuth, *The Art of Computer Programming*, Vol. 2, 3rd ed. (Addison-Wesley, 1997), §3.4.1."
source: src/ml/dist.c
---
**Algorithm.** `builtin_randomvariate` reads the distribution with `ml_read_dist` (which
declines a non-positive standard deviation or an inverted range up front, rather than
producing `NaN`s that would surface much later) and draws with `ml_draw`. A uniform
deviate is `a + (b − a) · random_uniform_01()`; a normal deviate is `μ + σ ·
ml_normal_deviate()`. `ml_normal_deviate` is Box–Muller in its **polar (Marsaglia)
form** — draw a point in the unit square, reject it unless `0 < s = u² + v² < 1`, and
return `u · √(−2 ln s / s)`; this needs no `sin`/`cos`, and the second deviate `v · √(…)`
is cached in a one-slot spare for the next call.

The decisive design point is that both draws go through `random_uniform_01`, the **same
stream `RandomReal` uses**, so `SeedRandom` makes them reproducible. A sampler with its
own generator would silently ignore `SeedRandom` while `RandomReal` honoured it, and
reproducibility that half-works is worse than none.

`RandomVariate[dist]` returns one scalar; `RandomVariate[dist, n]` returns a `List` of
`n` (`n` must be a non-negative integer, and `n = 0` is a valid empty request).

**Data structures.** A transient `MlDist`, a `double` output buffer, and the two static
`bm_have` / `bm_spare` Box–Muller cache slots (reset by `ml_dist_reset_cache`).

**Complexity / limits.** `O(n)` draws; the polar rejection accepts with probability
`π/4 ≈ 0.785` per pair. Only `NormalDistribution` and `UniformDistribution` are sampled;
each is also usable with no arguments for the standard case.
