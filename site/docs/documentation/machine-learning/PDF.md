# PDF

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PDF[dist, x] gives the probability density of dist at x, and threads over a list of x. Supports NormalDistribution and UniformDistribution.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]
Out[1]= {0.981398, -0.56572, 1.34033, 0.402313, -0.964221}

In[2]:= SeedRandom[1]; s = RandomVariate[NormalDistribution[5., 2.], 20000]; {Mean[s], StandardDeviation[s]}
Out[2]= {5.00479, 1.99627}

In[3]:= PDF[NormalDistribution[], {-1., 0., 1.}]
Out[3]= {0.241971, 0.398942, 0.241971}
```

### Applications (5)

The standard normal at its peak

```mathematica
In[4]:= PDF[NormalDistribution[], 0]
Out[4]= 0.398942
```

A scalar distribution threads over a list of points

```mathematica
In[5]:= PDF[NormalDistribution[], {-1., 0., 1.}]
Out[5]= {0.241971, 0.398942, 0.241971}
```

Zero outside the support, flat inside

```mathematica
In[6]:= PDF[UniformDistribution[{0., 2.}], {-1., 1., 3.}]
Out[6]= {0.0, 0.5, 0.0}
```

A fitted multinormal

```mathematica
In[7]:= m = LearnDistribution[{{1., 2.}, {2., 3.}, {3., 5.}, {4., 4.}, {5., 7.}, {6., 8.}}]
Out[7]= LearnedDistribution["Multinormal", <>]
```

A matrix threads to one density per row

```mathematica
In[8]:= PDF[m, {{3.5, 4.8}, {50., 50.}}]
Out[8]= {0.113186, 3.6193e-169}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/ml/dist.c`](https://github.com/stblake/mathilda/blob/main/src/ml/dist.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
- Tests: [`tests/test_ml_dist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_dist.c)

## Notes & additional examples

### Notes

`PDF[NormalDistribution[], 0]` is `1/Sqrt[2 Pi]` numerically, which is how the
closed-form density is verified. A uniform density is `1/(b − a)` on its closed support
and exactly zero outside.

The reading of a list argument depends on the distribution, and it has to. For a
**scalar** distribution a list of `x` threads to a list of densities — the shape a caller
plotting a density wants. For a **point** distribution (a fitted multinormal, mixture or
kernel estimate) the argument is itself a coordinate vector, so a list is *one*
observation, and only a *matrix* threads to one density per row. The second row above,
far out in the tail, comes back as a vanishingly small but non-zero density because the
point kernels work in log space.
