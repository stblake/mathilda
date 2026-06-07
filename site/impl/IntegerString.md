---
source: src/int.c
---
`builtin_integerstring` renders `|n|` as a base-`b` digit string (default 10, optional left-padded length). It calls GMP's `mpz_get_str`, which emits `'0'..'9','a'..'z'` for bases up to 62; the base is capped at 36 to match Mathematica's surface convention. A length argument left-pads with `'0'` or keeps only the low-order digits. Sign is discarded; `IntegerString[0]` is `"0"`. Validates arity (`IntegerString::argb`), numeric non-integer `n` (`::int`), base in `[2, 36]` (`::basf`), and non-negative machine length (`::intnn`); registered `Listable` so list inputs thread.
