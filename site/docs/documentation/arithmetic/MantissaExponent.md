# MantissaExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MantissaExponent[x] gives a list {m, e} containing the mantissa and exponent of the real number x, such that x = m * 10^e and 1/10 <= |m| < 1 (or m = 0 when x = 0).`**

**`MantissaExponent[x, b] gives the base-b mantissa and exponent; the mantissa then lies in 1/b <= |m| < 1.`**

<details>
<summary>Notes</summary>

Works for exact (Integer, Rational) and approximate (Real, MPFR) numeric inputs.  For exact inputs the mantissa is an exact Rational; for inexact inputs the mantissa carries the same precision as x.  Currently only integer bases \>= 2 are supported.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= MantissaExponent[3.4 10^30]
Out[1]= {0.34, 31}

In[2]:= MantissaExponent[456.1414]
Out[2]= {0.456141, 3}

In[3]:= MantissaExponent[123451]
Out[3]= {123451/1000000, 6}

In[4]:= MantissaExponent[1027, 2]
Out[4]= {1027/2048, 11}

In[5]:= MantissaExponent[2^100, 2]
Out[5]= {1/2, 101}

In[6]:= MantissaExponent[N[Pi, 30]]
Out[6]= {0.3141592653589793238462643383278, 1}

In[7]:= MantissaExponent[-3/2]
Out[7]= {-3/20, 1}
```

### Applications (5)

```mathematica
In[8]:= MantissaExponent[123.45]
Out[8]= {0.12345, 3}

In[9]:= MantissaExponent[7/3]
Out[9]= {7/30, 1}

In[10]:= MantissaExponent[1024, 2]
Out[10]= {1/2, 11}

In[11]:= MantissaExponent[N[Pi, 30]]
Out[11]= {0.3141592653589793238462643383278, 1}

In[12]:= MantissaExponent[255, 16]
Out[12]= {255/256, 2}
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

`builtin_mantissa_exponent` returns `{m, e}` with `x = m * b^e` and `1/b <= |m| < 1` (default base 10). It classifies the input via `rd_classify`. For **exact** inputs (Integer/BigInt/Rational) it works in a signed `mpq_t`: it finds the natural exponent `e` with `rd_rational_natural_exp`, scales numerator or denominator by `b^|e|`, canonicalises, and emits an exact `Rational` mantissa. For **machine reals** it computes `e = floor(log|x|/log b) + 1` then `m = x / b^e` with off-by-one corrections for log double-rounding; the **MPFR** path mirrors this at the input precision. `MantissaExponent[0]` is `{0, 0}`. Complex inputs emit `MantissaExponent::realx`; base `< 2` emits `::ibase`; non-integer bases are left unevaluated (only integer bases supported).

- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` numbers. The sign of `x` carries
  through to the mantissa.
- For exact inputs the mantissa is an exact `Rational` of the form
  `x / b^e`; for inexact inputs the mantissa is returned with the same
  precision as `x`.
- Bases must be integers `>= 2`. Non-integer bases leave the call
  unevaluated (general `Real` / symbolic bases are not yet supported).
- Complex arguments emit `MantissaExponent::realx` and leave the call
  unevaluated. Symbolic (non-numeric) arguments are left unevaluated
  silently.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/)

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_mantissa_exponent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mantissa_exponent.c)

## Notes & additional examples

### Notes

`MantissaExponent[x]` returns `{m, e}` with `x = m * 10^e` and `1/10 <= |m| < 1` (or `{0, 0}` when `x` is 0). `MantissaExponent[x, b]` uses base `b`, so `1/b <= |m| < 1`. Exact inputs keep an exact `Rational` mantissa (`7/3 -> {7/30, 1}`); inexact inputs keep their full working precision (the `N[Pi, 30]` mantissa carries 30 digits). Only integer bases `>= 2` are currently supported.
