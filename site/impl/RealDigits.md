---
source: src/real.c
---
**Algorithm.** `builtin_realdigits` returns the digit list of a real number, in the Mathematica `{digits, exponent}` form. It accepts `RealDigits[x]`, `RealDigits[x, b]`, `RealDigits[x, b, len]`, `RealDigits[x, b, len, p]` (1–4 args; wrong count emits `RealDigits::argb`). `x` is classified by `rd_classify`; concrete non-real `Complex` input emits `RealDigits::realx`, and symbolic constants (Pi, E, …) are numericalised only once enough precision context (base and length) is known. The base defaults to 10, must be an integer ≥ 2 (`RealDigits::ibase` otherwise) and fit in `unsigned long`. Digits are extracted by repeated scaled-floor / MPFR shifting in the requested base, honouring the optional length and starting-position arguments.

**Data structures.** GMP `mpz_t` for the base and integer parts; MPFR for the fractional digit extraction when built. Output is a `List` of digits paired with an integer exponent.
