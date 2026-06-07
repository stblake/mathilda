---
source: src/real.c
---
`builtin_real_exponent` returns `RealExponent[x]` / `RealExponent[x, b]` — essentially `⌊Log_b|x|⌋`, the exponent of the leading digit. It rejects true (non-zero-imaginary) `Complex` inputs (`RealExponent::realx`/`::ibase`) and bad arg counts (`RealExponent::argt`). Symbolic constants (Pi, E, …) and either argument are numericalised to a recognised numeric kind at a working precision lifted to cover any MPFR input (`+32` guard bits, so the downstream `Log` keeps precision), then the floor of the base-b logarithm of `|x|` is taken.
