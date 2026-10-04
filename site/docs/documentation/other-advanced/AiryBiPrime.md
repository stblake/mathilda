# AiryBiPrime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AiryBiPrime[z]`**

gives the derivative Bi'(z) of the Airy function AiryBi.

**`AiryBiPrime[0] = 3^(1/6)/Gamma[1/3], AiryBiPrime[+Infinity] = Infinity. Real`**

<details>
<summary>Notes</summary>

and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[AiryBiPrime\[z\], z\] = z AiryBi\[z\]. Listable.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= AiryBiPrime[0]
Out[1]= 3^(1/6)/Gamma[1/3]

In[2]:= N[AiryBiPrime[0], 40]
Out[2]= 0.44828835735382635791482371039882839086621

In[3]:= D[AiryBiPrime[z], z]
Out[3]= z AiryBi[z]

In[4]:= AiryBiPrime[1.0 + 1.0 I]
Out[4]= 0.0756628 + 0.783701*I
```

## Algorithm

Mathilda -- the Airy function Bi.

```text
  AiryBi[z]   Airy function Bi(z), the solution of  y'' = z y  that grows
              exponentially as z -> +Infinity along the real axis. Bi is an
              *entire* function of z (no branch cuts), the companion of Ai.
```

Evaluation is layered so each kind of argument takes the most accurate and cheapest route:

```text
  exact special values   ->  AiryBi[0] = 1/(3^(1/6) Gamma[2/3]),
                             AiryBi[+Infinity] = Infinity, AiryBi[-Infinity] = 0
  machine real           ->  unified complex-MPFR core at 53 bits, real part
  arbitrary real         ->  unified complex-MPFR core at mpfr_get_prec bits
  complex (any precision) ->  unified complex-MPFR core, Complex[..] result
  everything else        ->  stays symbolic (return NULL)
```

The unified core `airy_bi_core` evaluates Bi(z) and Bi'(z) together in a file-local complex-MPFR toolkit (`acx`, pairs of mpfr_t -- no MPC library is available; this mirrors the `acx`/`ecx`/`pcx`/`gcx` toolkits in airyai.c/ erf.c/polylog.c/gamma.c). It routes between three algorithms on r = |z|, theta = arg z, and the requested output precision P:

```text
  - Maclaurin series (small/moderate |z|, accurate everywhere). From
    Bi'' = z Bi the Taylor coefficients satisfy b_0 = Bi(0), b_1 = Bi'(0),
    b_2 = 0 and b_n = b_{n-3} / (n (n-1)) for n >= 3 -- identical recurrence
    to Ai, different seed constants. The partial sums reach magnitude
    ~exp((2/3) r^{3/2}) before cancelling for complex / negative arguments,
    so the core adds  (2/3) r^{3/2} / ln2  guard bits to absorb that exactly.

  - Dominant asymptotic series (large |z|, central sector), DLMF 9.7.7/9.7.8.
    With zeta = (2/3) z^{3/2}
        Bi(z)  ~ exp(zeta)/(sqrt(pi) z^{1/4}) Sum u_k / zeta^k,
        Bi'(z) ~ z^{1/4} exp(zeta)/sqrt(pi) Sum v_k / zeta^k,
    summed to the optimal (smallest-term) truncation. The u_k, v_k are the
    SAME coefficients as Ai's asymptotic series, but with no (-1)^k sign and
    prefactor 1/sqrt(pi) (not 1/(2 sqrt(pi))). The single dominant series is
    accurate only where the neglected recessive companion ~exp(-2 Re zeta)
    is below 2^-P, i.e. Re zeta = (2/3) r^{3/2} cos(3 theta/2) > (P ln2)/2.
    Bi's anti-Stokes line is |arg z| = pi/3, so this keeps the whole positive
    half-plane (including the exponentially large positive axis) at full
    precision.

  - Connection to Ai (large |z|, otherwise -- near and left of |arg z| = pi/3,
    covering the oscillatory negative real axis). DLMF 9.2.10:
        Bi(z)  = e^{ i pi/6} Ai(z e^{ 2 pi i/3}) + e^{-i pi/6} Ai(z e^{-2 pi i/3}),
        Bi'(z) = e^{i5pi/6} Ai'(z e^{ 2 pi i/3}) + e^{-i5pi/6} Ai'(z e^{-2 pi i/3}).
    The two rotated points have |arg| <= pi and |w| = |z| (large), so they are
    evaluated by a file-local Ai asymptotic kernel (direct series + Ai's own
    2 pi/3 connection wrapper, DLMF 9.2.12). The Bi oscillation on the
    negative axis emerges naturally from the two rotated Ai evaluations.
```

D[AiryBi[z], z] = AiryBiPrime[z] (see calculus/deriv.c); the Maclaurin series at 0 is produced by the generic Taylor-via-D path once AiryBi[0] / AiryBiPrime[0] have closed-form values.

AiryBiPrime[z] = Bi'(z) is a full numeric evaluator in its own right: because

```text
`airy_bi_core` returns Bi(z) and Bi'(z) together, AiryBiPrime reuses the very
```

same Maclaurin / asymptotic / connection machinery and simply selects the derivative component. Its exact values are AiryBiPrime[0] = 3^(1/6)/Gamma[1/3] and AiryBiPrime[+Infinity] = Infinity (Bi' is the dominant, growing solution); at -Infinity Bi' has no limit (oscillation with growing ~|z|^(1/4) amplitude) and is left unevaluated.

Attributes (both heads): Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- DLMF §9.7.7-9.7.8 (dominant asymptotic series for Bi, Bi') and §9.2.10 (the Bi-from-Ai connection formula).
- Source: [`src/special_functions/airybi.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/airybi.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`AiryBiPrime[z]` is the derivative `Bi'(z)`. Its exact origin value is
`3^(1/6)/Gamma[1/3]`, and a further derivative satisfies the Airy equation in
the form `D[AiryBiPrime[z], z] == z AiryBi[z]`. Complex arguments evaluate to
machine precision and, under `N[..., n]`, to arbitrary MPFR precision;
`AiryBiPrime` is Listable.
