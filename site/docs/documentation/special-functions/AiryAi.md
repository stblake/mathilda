# AiryAi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AiryAi[z]`**

gives the Airy function Ai(z), the solution of y'' = z y that decays as z -\> +Infinity.

**`AiryAi[0] = 1/(3^(2/3) Gamma[2/3]), AiryAi[+-Infinity] = 0. An entire`**

<details>
<summary>Notes</summary>

function of z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[AiryAi\[z\], z\] = AiryAiPrime\[z\]. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362
```

### Applications (5)

```mathematica
In[3]:= AiryAi[0]
Out[3]= 1/(3^(2/3) Gamma[2/3])

In[4]:= N[AiryAi[1], 40]
Out[4]= 0.13529241631288141552414742351546630617494

In[5]:= D[AiryAi[z], {z, 2}]
Out[5]= z AiryAi[z]

In[6]:= Series[AiryAi[z], {z, 0, 4}]
Out[6]= 1/(3^(2/3) Gamma[2/3]) + -1/(3^(1/3) Gamma[1/3]) z + 1/6/(3^(2/3) Gamma[2/3]) z^3 + -1/12/(3^(1/3) Gamma[1/3]) z^4 + O[z]^5

In[7]:= N[AiryAi[2] AiryBiPrime[2] - AiryAiPrime[2] AiryBi[2], 30]
Out[7]= 0.3183098861837906715377675267449
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

**Algorithm.** `builtin_airyai` takes the exact special values first — `AiryAi[0] = 1/(3^(2/3) Gamma[2/3])` (built as a symbolic expression and evaluated), `AiryAi[±Infinity] = 0`, `AiryAi[Indeterminate] = Indeterminate` — and otherwise routes numeric arguments (machine real, arbitrary-precision real, or complex with an inexact part) to the unified core `airy_ai_core`, which computes `Ai(z)` and `Ai'(z)` together. The core selects between two algorithms on `r = |z|` and the requested precision `P`: a **Maclaurin series** for small/moderate `|z|` (from `Ai'' = z Ai`, coefficients `a_0 = Ai(0)`, `a_1 = Ai'(0)`, `a_2 = 0`, `a_n = a_{n-3}/(n(n-1))`), with `(2/3) r^{3/2}/ln2` guard bits to absorb the partial-sum cancellation; and the **asymptotic series** DLMF 9.7.5/9.7.6 (`zeta = (2/3) z^{3/2}`) summed to optimal truncation for large `|z|` with `|arg z| ≤ 2π/3`. Closer to the negative real axis the core uses the DLMF 9.2.12 connection relation `Ai(z) = -[w Ai(wz) + conj(w) Ai(conj(w) z)]`, `w = e^{2πi/3}`, so the oscillation on the negative axis emerges from two rotated evaluations. `AiryAiPrime` reuses the same core and selects the derivative component; `D[AiryAi[z], z] = AiryAiPrime[z]` lives in `calculus/deriv.c`. A genuinely symbolic or exact non-zero argument (e.g. `AiryAi[2]`) stays unevaluated.

**Data structures.** `Expr`; the numeric core works in a file-local complex-MPFR toolkit `acx` (pairs of `mpfr_t`, no MPC library), each op at an explicit working precision and alias-safe. Machine results overflowing a `double` are promoted to MPFR. The ND kernel is a real `REG_U` registration (`NDKU_AiryAi` → `sf_machine_airy_ai` in `src/special_functions/sf_machine.c`), so a packed or visible real `NDArray` runs element-wise; the kernel is real-only and declines complex buffers. `Compile[]` lowers `AiryAi` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; the MPFR core scales its term count and guard bits with `P` and `|z|`. Entire function (no branch cuts). The exact value at `0` is closed form; symbolic/exact non-zero arguments stay symbolic. `AiryAiPrime` carries its own exact values (`-1/(3^(1/3) Gamma[1/3])` at 0, `0` at `+Infinity`); `-Infinity` is left unevaluated for the derivative since `Ai'` oscillates with growing amplitude there.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [N](../../arithmetic/N/), [AiryAiPrime](../../other-advanced/AiryAiPrime/)

- DLMF §9.2, §9.7 — Airy functions: power series, asymptotic expansions, and the connection formula Ai(z) = -[w Ai(wz) + conj(w) Ai(conj(w) z)].
- Source: [`src/special_functions/airyai.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/airyai.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

`AiryAi[z]` is the recessive solution of the Airy equation `y'' == z y` decaying
as `z -> +Infinity`; the second-derivative identity `D[AiryAi[z], {z, 2}] == z AiryAi[z]`
makes this explicit. The exact origin value is `1/(3^(2/3) Gamma[2/3])`, and the
Maclaurin series shows the characteristic missing `z^2` term (every third
coefficient vanishes). The last example is the Airy Wronskian
`Ai(z) Bi'(z) - Ai'(z) Bi(z) == 1/Pi`, recovered numerically as
`0.318309886...`. Real and complex arguments evaluate at machine or MPFR
precision; `AiryAi` is Listable.
