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
```

```mathematica
In[1]:= N[AiryBiPrime[0], 40]
Out[1]= 0.44828835735382635791482371039882839086616
```

```mathematica
In[1]:= D[AiryBiPrime[z], z]
Out[1]= z AiryBi[z]
```

```mathematica
In[1]:= AiryBiPrime[1.0 + 1.0 I]
Out[1]= 0.0756628 + 0.783701*I
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

Attributes (both heads): Listable, NumericFunction, Protected, ReadProtected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
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
