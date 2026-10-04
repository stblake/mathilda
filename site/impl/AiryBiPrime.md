---
references:
  - "DLMF §9.7.7-9.7.8 (dominant asymptotic series for Bi, Bi') and §9.2.10 (the Bi-from-Ai connection formula)."
source: src/special_functions/airybi.c
---
**Algorithm.** `builtin_airybiprime` evaluates `Bi'(z)` through exactly the same
unified core as `AiryBi`: `airy_bi_core` returns `Bi(z)` and `Bi'(z)` together, and
`AiryBiPrime` simply selects the derivative component. The core first takes the exact
special values — `AiryBiPrime[0] = 3^(1/6)/Gamma[1/3]`, `AiryBiPrime[+Infinity] =
Infinity` (Bi' is the dominant, growing solution), with `Indeterminate` passed through
and `-Infinity` deliberately left unevaluated (Bi' oscillates there with growing
`~|z|^{1/4}` amplitude, so it has no limit). For a machine real, an arbitrary-precision
`EXPR_MPFR`, or a complex argument with an inexact part, it routes by `r = |z|`,
`theta = arg z` and the requested precision `P` between three algorithms:

1. **Maclaurin series** (small/moderate `|z|`) — from `Bi'' = z Bi` the coefficients
   obey `b_n = b_{n-3}/(n(n-1))`; the core adds `(2/3) r^{3/2} / ln2` guard bits to
   absorb the cancellation against the `~exp((2/3) r^{3/2})` partial sums.
2. **Dominant asymptotic series** (large `|z|`, central sector, DLMF 9.7.8) —
   `Bi'(z) ~ z^{1/4} exp(zeta)/sqrt(pi) Sum v_k / zeta^k` with `zeta = (2/3) z^{3/2}`,
   summed to the smallest-term truncation; valid out to the anti-Stokes line
   `|arg z| = pi/3`, covering the whole positive half-plane.
3. **Connection to Ai** (large `|z|` elsewhere, DLMF 9.2.10) — two rotated Ai'
   evaluations recover the negative-axis oscillation.

The derivative rule `D[AiryBiPrime[z], z] = z AiryBi[z]` lives in `calculus/deriv.c`,
and `AiryBiPrime[0]` seeds the generic Taylor-via-`D` series of `AiryBi` at 0.

**Data structures.** A file-local complex-MPFR toolkit (`acx`, pairs of `mpfr_t` — no
MPC library is assumed), mirroring the `acx`/`ecx`/`gcx` toolkits in `airyai.c` /
`erf.c` / `gamma.c`. A machine real uses the 53-bit path. `AiryBiPrime` is registered
as a unary (`REG_U`) real ND kernel (`ndk_AiryBiPrime_r` → `sf_machine_airy_bi_prime`
in `src/ndkernels.c`), so it is packed-aware: a packed or visible `NDArray` of reals
runs element-wise through the machine kernel.

**Complexity / limits.** `O(1)` per element at machine precision; arbitrary precision
costs grow with the series length and the guard bits. Attributes: `Listable`,
`NumericFunction`, `Protected`. `CompileDiagnostics` reports `Compiled -> True` at both
scalar (`{{x, _Real}}`) and rank-1 (`{{v, _Real, 1}}`) shapes. Non-numeric arguments
stay symbolic.
