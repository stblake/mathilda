# AiryAiPrime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AiryAiPrime[z]`**

gives the derivative Ai'(z) of the Airy function AiryAi.

**`AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]), AiryAiPrime[+Infinity] = 0. Real`**

<details>
<summary>Notes</summary>

and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[AiryAiPrime\[z\], z\] = z AiryAi\[z\]. Listable.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= AiryAiPrime[0]
Out[1]= -1/(3^(1/3) Gamma[1/3])

In[2]:= N[AiryAiPrime[0], 40]
Out[2]= -0.25881940379280679840518356018920396347907

In[3]:= D[AiryAiPrime[z], z]
Out[3]= z AiryAi[z]

In[4]:= AiryAiPrime[1.0 + 1.0 I]
Out[4]= -0.130628 + 0.163068*I
```

## Algorithm

Mathilda -- the Airy function Ai.

```text
  AiryAi[z]   Airy function Ai(z), the solution of  y'' = z y  that tends to
              zero as z -> +Infinity along the real axis. Ai is an *entire*
              function of z (no branch cuts).
```

Evaluation is layered so each kind of argument takes the most accurate and cheapest route:

```text
  exact special values   ->  AiryAi[0] = 1/(3^(2/3) Gamma[2/3]),
                             AiryAi[+-Infinity] = 0
  machine real           ->  unified complex-MPFR core at 53 bits, real part
  arbitrary real         ->  unified complex-MPFR core at mpfr_get_prec bits
  complex (any precision) ->  unified complex-MPFR core, Complex[..] result
  everything else        ->  stays symbolic (return NULL)
```

The unified core `airy_ai_core` evaluates Ai(z) and Ai'(z) together in a file-local complex-MPFR toolkit (`acx`, pairs of mpfr_t -- no MPC library is available; this mirrors the `ecx`/`pcx`/`gcx` toolkits in erf.c/polylog.c/ gamma.c). It routes between two algorithms on r = |z| and the requested output precision P:

```text
  - Maclaurin series (small/moderate |z|). From Ai'' = z Ai the Taylor
    coefficients satisfy a_0 = Ai(0), a_1 = Ai'(0), a_2 = 0 and
    a_n = a_{n-3} / (n (n-1)) for n >= 3. The partial sums reach magnitude
    ~exp((2/3) r^{3/2}) before cancelling for complex / negative arguments,
    so the core adds  (2/3) r^{3/2} / ln2  guard bits to absorb that exactly.

  - Asymptotic series (large |z|), DLMF 9.7.5/9.7.6. With zeta = (2/3) z^{3/2}
        Ai(z)  ~ exp(-zeta)/(2 sqrt(pi) z^{1/4}) Sum (-1)^k u_k / zeta^k,
        Ai'(z) ~ -z^{1/4} exp(-zeta)/(2 sqrt(pi)) Sum (-1)^k v_k / zeta^k,
    summed to the optimal (smallest-term) truncation. The single series is
    accurate for |arg z| <= 2 pi / 3; closer to the negative real axis
    (2 pi / 3 < |arg z| <= pi) the core uses the connection relation
    (DLMF 9.2.12)
        Ai(z) = -[ w Ai(w z) + conj(w) Ai(conj(w) z)],  w = e^{2 pi i / 3},
    which maps the argument into two points with |arg| <= 2 pi / 3 where the
    direct series is accurate; the oscillation on the negative real axis then
    emerges naturally from the sum of the two rotated evaluations.
```

D[AiryAi[z], z] = AiryAiPrime[z] (see calculus/deriv.c); the Maclaurin series at 0 is produced by the generic Taylor-via-D path once AiryAi[0] / AiryAiPrime[0] have closed-form values.

AiryAiPrime[z] = Ai'(z) is a full numeric evaluator in its own right: because

```text
`airy_ai_core` returns Ai(z) and Ai'(z) together, AiryAiPrime reuses the very
```

same Maclaurin / asymptotic / connection machinery and simply selects the derivative component. Its exact values are AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]) and AiryAiPrime[+Infinity] = 0; at -Infinity Ai' has no limit (oscillation with growing ~|z|^(1/4) amplitude) and is left unevaluated.

Attributes (both heads): Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`AiryAiPrime[z]` is the derivative `Ai'(z)`. Its exact origin value is
`-1/(3^(1/3) Gamma[1/3])`, and differentiating once more recovers the Airy
equation in the form `D[AiryAiPrime[z], z] == z AiryAi[z]`. Complex arguments
evaluate to machine precision (and to arbitrary precision under `N[..., n]`);
`AiryAiPrime` is Listable.
