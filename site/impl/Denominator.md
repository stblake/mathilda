---
source: src/rat.c
---
`builtin_denominator` calls `extract_num_den` and returns the denominator part (freeing the numerator). `extract_num_den` recognises `Rational[n, d]` (returns `d`); `Complex` (clears to a common integer denominator); `Power[b, e]`/`Exp[e]` with a superficially-negative exponent or a `Plus` exponent split into positive/negative pieces (the negative-exponent base becomes the denominator); and `Times`, which recurses into each factor and multiplies the collected denominators. A factor with no denominator contributes `1`. `Numerator` in the same file is the symmetric accessor.
