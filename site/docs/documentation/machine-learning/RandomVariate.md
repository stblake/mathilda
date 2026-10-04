# RandomVariate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomVariate[dist] draws one value from dist; RandomVariate[dist, n] draws a list of n. Supports NormalDistribution[mu, sigma] and UniformDistribution[{lo, hi}], each also usable with no arguments for the standard case. Draws come from the same stream as RandomReal, so SeedRandom makes them reproducible. A non-positive standard deviation, or an inverted range, returns unevaluated rather than producing NaNs.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]
Out[1]= {0.981398, -0.56572, 1.34033, 0.402313, -0.964221}

In[2]:= SeedRandom[7]; RandomVariate[UniformDistribution[{0., 10.}], 4]
Out[2]= {0.553604, 1.72116, 7.17576, 4.2721}

In[3]:= RandomVariate[NormalDistribution[], 0]
Out[3]= {}
```

### Applications (3)

Reproducible once the stream is seeded

```mathematica
In[4]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]
Out[4]= {0.981398, -0.56572, 1.34033, 0.402313, -0.964221}
```

Four uniform draws on 0..10

```mathematica
In[5]:= SeedRandom[7]; RandomVariate[UniformDistribution[{0., 10.}], 4]
Out[5]= {0.553604, 1.72116, 7.17576, 4.2721}
```

An empty request is valid

```mathematica
In[6]:= RandomVariate[NormalDistribution[], 0]
Out[6]= {}
```

## Implementation notes

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

- **Draws come from the same stream as `RandomReal`**, so `SeedRandom` makes them
  reproducible. A sampler with its own generator would silently ignore `SeedRandom`
  while `RandomReal` honoured it — reproducibility half-working is worse than not
  working.
- Normal deviates use Box–Muller in its polar form, which needs no `sin`/`cos`. There
  was no Gaussian deviate anywhere in the tree before this.
- **A non-positive standard deviation, or an inverted range, returns unevaluated** —
  not `NaN`, which would propagate silently through a whole sample and surface much
  later as a strange plot. `RandomVariate[dist, 0]` is a valid empty request.

**Attributes:** `Protected`.

## References

**See also:** [RandomReal](../../random-number-generation/RandomReal/), [SeedRandom](../../random-number-generation/SeedRandom/)

- G. Marsaglia and T. A. Bray, *A convenient method for generating normal variables*, SIAM Review **6** (1964) 260-264 (the polar method).
- D. E. Knuth, *The Art of Computer Programming*, Vol. 2, 3rd ed. (Addison-Wesley, 1997), §3.4.1.
- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`RandomVariate[dist]` draws one value and `RandomVariate[dist, n]` a list of `n`.
Normal, uniform, and their argument-free standard cases are supported.

Draws come from the **same stream as `RandomReal`**, so `SeedRandom` makes them
reproducible — a sampler with its own generator would silently ignore `SeedRandom` while
`RandomReal` honoured it, and reproducibility that half-works is worse than none. The
seeded lines above therefore print the same values on every run. Normal deviates use
Box–Muller in its polar form, which needs no `sin`/`cos`.

A non-positive standard deviation, or an inverted range, returns unevaluated rather than
producing `NaN`s that would propagate silently through a whole sample; `n = 0` is a valid
empty request.
