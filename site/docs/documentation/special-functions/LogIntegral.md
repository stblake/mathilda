# LogIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LogIntegral[z]`**

gives the logarithmic integral li(z), the principal value of Integral\_0^z dt/ln t, equal to ExpIntegralEi\[Log\[z\]\], with a branch cut on (-Infinity, 1). LogIntegral\[0\] = 0, LogIntegral\[1\] = -Infinity,

**`LogIntegral[Infinity] = Infinity. Real and complex inputs evaluate`**

<details>
<summary>Notes</summary>

numerically at machine or arbitrary (MPFR) precision; D\[LogIntegral\[z\], z\] = 1/Log\[z\]. Listable.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[LogIntegral[2], 40]
Out[1]= 1.0451637801174927848445888891946131365227
```

```mathematica
In[1]:= D[LogIntegral[z], z]
Out[1]= 1/Log[z]
```

```mathematica
In[1]:= N[LogIntegral[10^6], 30]
Out[1]= 78627.54915946218191986291074769
```

```mathematica
In[1]:= N[LogIntegral[1000], 20]
Out[1]= 177.609657990152226688
```

## Algorithm

Mathilda -- the logarithmic integral li.

```text
  LogIntegral[z]   li(z) = PV Int_0^z dt / ln t
```

li has a branch cut along (-Infinity, +1); the principal value is taken on the cut. The implementation rests on the identity

```text
  li(z) = Ei(Log z),
```

where Ei is ExpIntegralEi and Log is the principal logarithm. This lets us reuse ExpIntegralEi's fully-tested numeric stack (mpfr_eint / the real and complex convergent series with cancellation guard bits) without duplicating any of it, and the principal Log automatically supplies the +-i Pi jump that places the branch cut on (-Infinity, +1).

Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, -Infinity, Infinity, Indeterminate
  numeric (inexact) z       ->  evaluate ExpIntegralEi[Log[z]]
  everything else           ->  stays symbolic (return NULL)
```

Exact non-special numbers (e.g. LogIntegral[2], LogIntegral[1/2]) stay symbolic, matching the Wolfram Language; only inexact input or an explicit N[...] (which rewrites the argument to an MPFR number) evaluates numerically.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

- Exact special values: `LogIntegral[0] = 0`, `LogIntegral[1] = -Infinity`,
  `LogIntegral[Infinity] = Infinity`; `ComplexInfinity` and `Indeterminate` map
  to `Indeterminate`.
- Exact non-special arguments stay symbolic (`LogIntegral[2]`,
  `LogIntegral[1/2]`); numeric values follow from a `Real`/MPFR argument or from
  `N`.
- Numeric evaluation (machine *and* arbitrary precision) routes through
  `ExpIntegralEi[Log[z]]`:
  - Real z > 1 (`Log z > 0`) → MPFR `mpfr_eint`, correctly rounded and fast even
    at very high precision: `LogIntegral[20.] = 9.9053`, `LogIntegral[2.] = 1.04516`,
    `N[LogIntegral[2], 50] = 1.0451637801174927848445888891946131365226155781512`.
  - Real 0 < z < 1 (`Log z < 0`) → the on-cut convergent series, returning a
    **real** principal value: `LogIntegral[0.5] = -0.378671`,
    `LogIntegral[1.2] = -0.933787`.
  - **Complex** (and real z < 0, whose principal `Log` is complex) → the complex
    series with guard bits, so machine-precision complex results are fully
    accurate, e.g. `LogIntegral[2. + I] = 1.41126 + 1.22471 I`,
    `N[Re[LogIntegral[2 + I]], 30] = 1.41125904201780100568439320706`.
- Derivative: `D[LogIntegral[z], z] = 1/Log[z]` (chain rule applies, e.g.
  `D[LogIntegral[x^2], x] = (2 x)/Log[x^2]`); the Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `LogIntegral::argx` and stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[ExpIntegralEi](../../special-functions/ExpIntegralEi/), [Log](../../elementary-functions/Log/), [N](../../arithmetic/N/), [D](../../calculus/D/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_cherry_sigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_sigma.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_integrate_risch_transcendental.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_risch_transcendental.c)

## Notes & additional examples

### Notes

`LogIntegral[z]` is the logarithmic integral `li(z)`, the principal value of `Integral_0^z dt/Log[t]`, equal to `ExpIntegralEi[Log[z]]`, with a branch cut on `(-Infinity, 1)`. Its derivative is `1/Log[z]`. `li(x)` is the leading term of the prime-counting approximation `PrimePi[x] ~ li(x)`; for example `li(10^6)` is about `78627.5`, close to `PrimePi[10^6] = 78498`. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision. `LogIntegral[1] = -Infinity`. Listable.
