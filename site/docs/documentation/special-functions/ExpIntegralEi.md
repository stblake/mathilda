# ExpIntegralEi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpIntegralEi[z]`**

gives the exponential integral Ei(z), the principal value of -Integral\_{-z}^Infinity e^-t/t dt, with a branch cut on (-Infinity, 0).

**`ExpIntegralEi[0] = -Infinity, ExpIntegralEi[Infinity] = Infinity,`**

**`ExpIntegralEi[-Infinity] = 0, ExpIntegralEi[+-I Infinity] = +-I Pi. Real and`**

**`D[ExpIntegralEi[z], z] = E^z/z. Listable.`**

<details>
<summary>Notes</summary>

complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= ExpIntegralEi[0]
Out[1]= -Infinity

In[2]:= D[ExpIntegralEi[z], z]
Out[2]= E^z/z

In[3]:= N[ExpIntegralEi[1], 40]
Out[3]= 1.8951178163559367554665209343316342690171

In[4]:= N[ExpIntegralEi[I], 30]
Out[4]= 0.3374039229009681346626462038893 + 2.516879397162079634172675005462*I
```

## Algorithm

Mathilda -- the exponential integral Ei.

```text
  ExpIntegralEi[z]   Ei(z) = -PV Int_{-z}^Inf e^-t/t dt
```

Ei has a branch cut along the negative real axis (-Inf, 0); the principal value is taken on the cut. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values    ->  -Infinity, Infinity, 0, +-I Pi, Indeterminate
  machine real x > 0       ->  MPFR mpfr_eint (correctly rounded)
  machine real x < 0       ->  real convergent series with ln|x|
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the convergent series in MPFR with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent series (DLMF 6.6.2) is

```text
  Ei(z) = gamma + Log(z) + Sum_{k>=1} z^k / (k k!),
```

valid for all z != 0. On the real negative axis the real part ln|x| is used (principal value on the cut); for genuinely complex z the principal Log(z) supplies the +-i Pi jump across the cut. The partial sums of the series can reach magnitude ~e^|z| before the answer (~e^Re z) emerges, so the MPFR path adds (|z| + |Re z|)/ln2 guard bits to absorb that cancellation exactly. For x > 0 MPFR's native mpfr_eint is used directly (fast and correctly rounded, which also covers very-high-precision arguments).

Attributes: Listable, NumericFunction, Protected.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |
| integrate 1/(1+Exp[x]) | 2.43 s | 0.109 s | 3.31 s |
| integrate Log[x]^3 | 1.59 s | 1.85 s | 10.6 s |
| integrate x Exp[x^2] | 0.685 s | 0.087 s | 4.42 s |

## Implementation notes

- Exact special values: `ExpIntegralEi[0] = -Infinity`,
  `ExpIntegralEi[Infinity] = Infinity`, `ExpIntegralEi[-Infinity] = 0`,
  `ExpIntegralEi[I Infinity] = I Pi`, `ExpIntegralEi[-I Infinity] = -I Pi`;
  `ComplexInfinity` and `Indeterminate` map to `Indeterminate`.
- Exact non-zero arguments stay symbolic (`ExpIntegralEi[2]`, `ExpIntegralEi[1/2]`);
  numeric values follow from a `Real`/MPFR argument or from `N`.
- Numeric evaluation (machine *and* arbitrary precision):
  - Real x > 0 → MPFR `mpfr_eint` (correctly rounded, fast even at very high
    precision: `ExpIntegralEi[1.8] = 4.24987`, `ExpIntegralEi[2.] = 4.95423`,
    `N[ExpIntegralEi[2], 50] = 4.9542343560018901633795051302270352755180535624200`).
  - Real x < 0 → the on-cut convergent series Ei(x) = γ + ln|x| + Σ xᵏ/(k·k!),
    in MPFR with `|x|/ln2` guard bits to absorb the partial-sum cancellation, and
    returns a **real** principal value (`ExpIntegralEi[-1.] = -0.219384`,
    `ExpIntegralEi[-5.] = -0.00114830`).
  - **Complex** → the same series with the principal `Log(z)`, evaluated in MPFR
    with `(|z| + |Re z|)/ln2` guard bits, so machine-precision complex results are
    fully accurate, e.g. `ExpIntegralEi[2. + I] = 4.06998 + 3.40094 I`,
    `N[ExpIntegralEi[2 + I], 30] = 4.06998094789392774228769025521 + 3.40094396980012162163040462603 I`.
    Approaching the cut from above gives +I Pi, from below −I Pi
    (`ExpIntegralEi[-1. + 10^-10 I] ≈ -0.219384 + 3.14159 I`). A `double complex`
    series is the `USE_MPFR=0` fallback.
  - **Large |z|** (real or complex) → the convergent series is infeasible
    (it would need ~`|z|/ln2` guard bits and ~`2|z|` terms), so once `|z|`
    exceeds roughly `prec·ln2` the **asymptotic expansion** takes over:
    `Ei(z) ~ (e^z/z) Σ_{k≥0} k!/z^k + i π sign(Im z)`, summed to its smallest
    term (DLMF 6.12.2; the `i π sign(Im z)` constant is the branch jump, so
    `ExpIntegralEi[±I Infinity] = ±I Pi` is recovered in the limit). The two
    regimes overlap and agree (`ExpIntegralEi[-50 I]` from either route gives
    `-0.00562839 − 3.12241 I`). This also fixes an `mpfr_init2` abort: the old
    guard-bit count `(long)(|z|/ln2)` overflowed for astronomically large `|z|`
    (e.g. `N[ExpIntegralEi[-10^60 I] + I Pi, 20]`).
- Derivative: `D[ExpIntegralEi[z], z] = E^z/z` (chain rule applies, e.g.
  `D[ExpIntegralEi[x^2], x] = (2 E^x^2)/x`); the origin Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `ExpIntegralEi::argx` and stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_expintegralei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expintegralei.c)

## Notes & additional examples

### Notes

`ExpIntegralEi[z]` is the exponential integral `Ei(z)`, with a branch cut on
`(-Infinity, 0)` and derivative `E^z/z`. On the imaginary axis it ties to the
cosine/sine integrals via `Ei(I) = Ci(1) + I (Pi/2 + Si(1))`. Real and complex
arguments evaluate at machine or arbitrary (MPFR) precision. Listable.
