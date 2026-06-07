---
source: src/real.c
---
`builtin_mantissa_exponent` returns `{m, e}` with `x = m * b^e` and `1/b <= |m| < 1` (default base 10). It classifies the input via `rd_classify`. For **exact** inputs (Integer/BigInt/Rational) it works in a signed `mpq_t`: it finds the natural exponent `e` with `rd_rational_natural_exp`, scales numerator or denominator by `b^|e|`, canonicalises, and emits an exact `Rational` mantissa. For **machine reals** it computes `e = floor(log|x|/log b) + 1` then `m = x / b^e` with off-by-one corrections for log double-rounding; the **MPFR** path mirrors this at the input precision. `MantissaExponent[0]` is `{0, 0}`. Complex inputs emit `MantissaExponent::realx`; base `< 2` emits `::ibase`; non-integer bases are left unevaluated (only integer bases supported).
