---
source: src/int.c
---
`builtin_digitcount` tallies base-`b` digit occurrences of `|n|` (default base 10). The 3-argument form `DigitCount[n, b, d]` returns the scalar count of digit `d` (`dc_count_one_digit`); the 1/2-argument form returns the histogram `{c(1), c(2), ..., c(b-1), c(0)}` — digit 0 last — built by `dc_build_histogram` into a `calloc`'d `int64` array. Validates arity (`DigitCount::argb`), numeric non-integer `n` (`::int`), base `>= 2` (`::base`, also for fractional bases), digit in `[0, base)` (`::digit`), and caps the list-form base to avoid OOM (`::ovfl`); symbolic `n` returns `NULL`.
