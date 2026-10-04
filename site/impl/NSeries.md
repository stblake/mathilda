---
source: src/numerical_calculus/nseries.c
references:
  - "J. N. Lyness and G. Sande, *Algorithm 413: ENTCAF and ENTCRE: evaluation of normalized Taylor coefficients of an analytic function*, Comm. ACM **14** (1971) 669–675."
  - "F. Bornemann, *Accuracy and stability of computing high-order derivatives of analytic functions by Cauchy integrals*, Found. Comput. Math. **11** (2011) 1–63."
---
**Algorithm.** `builtin_nseries` evaluates `NSeries[f, {x, x0, n}]`, the
numerical Taylor/Laurent expansion with terms `(x−x0)^{−n}` through `(x−x0)^n`,
returned as a `SeriesData`. `f` is sampled at `N` equispaced points on a circle
of radius `r` (`Radius`, default 1) centred at `x0`, `z_j = x0 + r
e^{2πi j/N}`, and a **discrete Fourier transform of the samples recovers the
Laurent coefficients by Cauchy's integral formula**, `c_k = (1/N) Σ_j f(z_j)
e^{−2πi j k/N}`, `a_e = c_{e mod N} · r^{−e}`. The upper DFT bins (`k = N−m`)
supply the negative-power coefficients, so one transform yields both the
principal and the analytic part — which is why it expands poles and **essential
singularities** (e.g. `Sin[x + 1/x]`) that symbolic `Series` cannot. It is exact
when `f` is analytic on an annulus containing the circle, and fails only if the
disk contains a branch cut. `N = 2^{⌈log2 n⌉ + 2}` (4× oversampling so the
leading aliased term sits below the round-off floor), clamped to `N > 2n`. Each
sample binds `x` Block-style and bumps the evaluation clock.

It is a **direct `O(N²)` DFT, deliberately not an FFT**: `N` is at most a few
hundred and each sample requires a full symbolic evaluation of `f` that
dominates the runtime, and one code path must serve both the machine and the
MPFR computations. The retained `2n+1` coefficients are refined by **adaptive
doubling** of `N` (up to 4 times) until they stop moving within the combined
tolerance; the change between successive `N` is the aliasing/round-off error
estimate, and a shortfall warns `NSeries::accgl`.

**Data structures.** The machine path keeps the samples in a `double _Complex*`
buffer and accumulates the DFT directly; the MPFR path keeps paired `mpfr_t*`
real/imag sample buffers, with angles from `mpfr_const_pi`/`mpfr_sin_cos` and the
scale `r^{−e}` via `mpfr_pow_si`. `WorkingPrecision` (`MachinePrecision` or a
digit count) selects the path. The result is a `SeriesData[x, x0, {a_{-n} … a_n},
−n, n+1, 1]` (denominator 1 — no fractional powers).

**Complexity / limits.** `O(N²)` arithmetic per pass but `N` full symbolic
samples dominate, with up to 4 doublings. Round-off breaks the conjugate
symmetry of a real-coefficient function, leaving tiny spurious residuals that are
**not** recognised as zero — `Chop` when needed; for a Laurent series the
`SeriesData` neglects higher-order poles, and no effort is made to justify the
coefficients' precision. Options: `Radius` (default 1, must be positive),
`WorkingPrecision` (default `MachinePrecision`), `AccuracyGoal` (default
`MachinePrecision`), `PrecisionGoal` (default Automatic). The spec requires `x` a
symbol and `n` a non-negative integer; `NSeries::badopt`/`ivar`/`nnum`/`accgl`
route through `mth_message`.
