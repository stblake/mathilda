# MantissaExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MantissaExponent[x] gives a list {m, e} containing the mantissa and exponent of the real number x, such that x = m * 10^e and 1/10 <= |m| < 1 (or m = 0 when x = 0).
MantissaExponent[x, b] gives the base-b mantissa and exponent; the mantissa then lies in 1/b <= |m| < 1.
Works for exact (Integer, Rational) and approximate (Real, MPFR) numeric inputs.  For exact inputs the mantissa is an exact Rational; for inexact inputs the mantissa carries the same precision as x.  Currently only integer bases >= 2 are supported.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

`builtin_mantissa_exponent` returns `{m, e}` with `x = m * b^e` and `1/b <= |m| < 1` (default base 10). It classifies the input via `rd_classify`. For **exact** inputs (Integer/BigInt/Rational) it works in a signed `mpq_t`: it finds the natural exponent `e` with `rd_rational_natural_exp`, scales numerator or denominator by `b^|e|`, canonicalises, and emits an exact `Rational` mantissa. For **machine reals** it computes `e = floor(log|x|/log b) + 1` then `m = x / b^e` with off-by-one corrections for log double-rounding; the **MPFR** path mirrors this at the input precision. `MantissaExponent[0]` is `{0, 0}`. Complex inputs emit `MantissaExponent::realx`; base `< 2` emits `::ibase`; non-integer bases are left unevaluated (only integer bases supported).

- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MantissaExponent[123.45]
Out[1]= {0.12345, 3}
```

```mathematica
In[1]:= MantissaExponent[7/3]
Out[1]= {7/30, 1}
```

```mathematica
In[1]:= MantissaExponent[1024, 2]
Out[1]= {1/2, 11}
```

```mathematica
In[1]:= MantissaExponent[N[Pi, 30]]
Out[1]= {0.3141592653589793238462643383278, 1}
```

```mathematica
In[1]:= MantissaExponent[255, 16]
Out[1]= {255/256, 2}
```

### Notes

`MantissaExponent[x]` returns `{m, e}` with `x = m * 10^e` and `1/10 <= |m| < 1` (or `{0, 0}` when `x` is 0). `MantissaExponent[x, b]` uses base `b`, so `1/b <= |m| < 1`. Exact inputs keep an exact `Rational` mantissa (`7/3 -> {7/30, 1}`); inexact inputs keep their full working precision (the `N[Pi, 30]` mantissa carries 30 digits). Only integer bases `>= 2` are currently supported.
