---
source: src/int.c
---
`builtin_integerdigits` returns the base-`b` digit list of `|n|` (default base 10, optional fixed length). It coerces `n` and `base` to GMP, then peels least-significant digits with a `mpz_tdiv_qr` divmod loop into a geometrically-grown `mpz_t` buffer (`O(log_b |n|)`), and emits them most-significant-first into a `List`. A length argument left-pads with zeros (or keeps only the low-order digits if shorter). Validates arity (`IntegerDigits::argb`), non-integer numeric `n` (`::int`), base `>= 2` (`::ibase`), and a non-negative machine-sized length (`::intnn`); symbolic `n` flows through as `NULL`.
