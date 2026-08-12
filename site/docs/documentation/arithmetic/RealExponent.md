# RealExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RealExponent[x] gives Log[10, |x|] -- the base-10 real exponent of x.`**

**`RealExponent[x, b] gives Log[b, |x|] in the specified base b.`**

<details>
<summary>Notes</summary>

Accepts Integer, BigInt, Rational, Real, and (with USE\_MPFR) MPFR inputs, plus symbolic numeric values such as Pi, E, or Pi^Pi.  Result is a machine Real unless an MPFR input lifts it to MPFR at that precision.  Exact zero gives -Infinity; machine 0. gives Log\[b, $MinMachineNumber\] (~ -307.65 in base 10); MPFR 0 with precision p digits gives -p / Log10\[b\].  Threads over lists.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= RealExponent[123.456]
Out[1]= 2.09151

In[2]:= RealExponent[123.456, 2]
Out[2]= 6.94785

In[3]:= RealExponent[N[Pi, 32]]
Out[3]= 0.497149872694133854351268288290899

In[4]:= RealExponent[Pi, E]
Out[4]= 1.14473

In[5]:= RealExponent[987654321/123456789]
Out[5]= 0.90309

In[6]:= RealExponent[{1, 2, 3, 4, 5}]
Out[6]= {0.0, 0.30103, 0.477121, 0.60206, 0.69897}

In[7]:= Table[RealExponent[Pi, b], {b, {2, 3, 5, 7, 10}}]
Out[7]= {1.6515, 1.04198, 0.711261, 0.588275, 0.49715}

In[8]:= RealExponent[0]
Out[8]= -Infinity
```

### Applications (6)

```mathematica
In[1]:= RealExponent[1234.5]
Out[1]= 3.09149
```

In an explicit base the result is the exact logarithm — the base-2 exponent of a
pure power of two is an integer, which makes `RealExponent` a quick way to count
the decimal digits of a huge integer (`Floor[log10] + 1`):

```mathematica
In[1]:= RealExponent[2^100, 2]
Out[1]= 100.0

In[2]:= Floor[RealExponent[2^1000]] + 1
Out[2]= 302
```

It accepts symbolic numeric values and lifts to MPFR precision when given an MPFR
argument — for example `Log10[E]` to 40 digits:

```mathematica
In[1]:= RealExponent[N[Pi^Pi]]
Out[1]= 1.56184

In[2]:= RealExponent[N[E, 40]]
Out[2]= 0.4342944819032518276511289189166050822944
```

`RealExponent` threads over lists:

```mathematica
In[1]:= RealExponent[{10, 100, 1000}]
Out[1]= {1.0, 2.0, 3.0}
```

## Algorithm

real.c

RealDigits builtin -- positional-notation digit expansion.

```text
  RealDigits[x]               default base 10, length set by Precision[x].
  RealDigits[x, b]            base b, length set by Precision[x] / Log10[b].
  RealDigits[x, b, len]       exactly `len` digits, MSD-first.
  RealDigits[x, b, len, n]    `len` digits, first one = coefficient of b^n.

Result form is `{ digits-list, exp }`.  The first element of digits-list is
the coefficient of b^(exp-1).  Sign of x is discarded.  Exact rationals
```

with non-terminating base-b expansions return a list ending in a nested

```text
list of the recurring block.  Inexact (machine or MPFR) reals get
```

Indeterminate for any requested digit beyond the available precision.

```text
  x can be: Integer, BigInt, Rational[n,d], Real (machine), or
            EXPR_MPFR (arbitrary precision; USE_MPFR builds only).
```

The general algorithm scales |x| by base^(-low) where `low` is the lowest digit position we need, floors to an integer N, and reads off the base-b

```text
digits of N (padding with leading zeros as needed).  This single GMP /
MPFR shift handles every numeric type uniformly.  For the special case of
```

an exact rational with no explicit `len`, a remainder-tracked long division detects terminating vs recurring expansions and emits the nested-list form.

```text
Implementation only supports integer bases b >= 2.  Non-integer bases
```

(e.g. GoldenRatio) emit a `::ibase` diagnostic and leave the call unevaluated -- adding them requires MPFR floor-iteration and has been deferred.

## Implementation notes

`builtin_real_exponent` returns `RealExponent[x]` / `RealExponent[x, b]` — essentially `⌊Log_b|x|⌋`, the exponent of the leading digit. It rejects true (non-zero-imaginary) `Complex` inputs (`RealExponent::realx`/`::ibase`) and bad arg counts (`RealExponent::argt`). Symbolic constants (Pi, E, …) and either argument are numericalised to a recognised numeric kind at a working precision lifted to cover any MPFR input (`+32` guard bits, so the downstream `Log` keeps precision), then the floor of the base-b logarithm of `|x|` is taken.

- `Protected`, `Listable`. Threads over lists in any argument position.
- Accepts `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` inputs.  Symbolic numeric
  arguments (`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`,
  or any numeric-valued composite such as `Pi^Pi` or `1/Pi`) are
  numericalized at the combined working precision before computation.
  Plain symbols with no numeric value are left unevaluated.
- Output is a machine `Real` unless one of the inputs already carries
  MPFR precision, in which case the result is `MPFR` at the higher of
  the input precisions.  This matches Mathematica's contagion: an
  explicit `N[..., p]` lifts the exponent to the same precision.
- The base must be a real number `> 1`; non-positive, `<= 1`, or complex
  bases emit `RealExponent::ibase` and leave the call unevaluated.
- Complex arguments with non-zero imaginary part emit
  `RealExponent::realx` and leave the call unevaluated.
- Sign of `x` is discarded.
- Zero handling (Mathematica-compatible):
  - Exact zero (Integer 0, BigInt 0, Rational 0/n) → `-Infinity`.
  - Machine `0.` → `Log[b, $MinMachineNumber]` (`≈ -307.65` for base 10).
  - MPFR `0``p` of `p` digits → `-p / Log10[b]` (`-p` for base 10).

**Attributes:** `Listable`, `Protected`.

## See also

[Rational](../../arithmetic/Rational/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [EulerGamma](../../mathematical-constants/EulerGamma/), [Catalan](../../mathematical-constants/Catalan/), [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Degree](../../mathematical-constants/Degree/)

## References

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_real_exponent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_real_exponent.c)

## Notes & additional examples

### Notes

`RealExponent[x]` gives `Log[10, |x|]`; `RealExponent[x, b]` gives `Log[b, |x|]`.
It accepts `Integer`, `BigInt`, `Rational`, `Real`, MPFR, and symbolic numeric
values such as `Pi`, `E`, or `Pi^Pi`. The result is a machine `Real` unless an
MPFR input lifts it to that precision. Exact zero gives `-Infinity`; the sign of
`x` is ignored.
