---
source: src/int.c
---
`builtin_fromdigits` is the inverse of `IntegerDigits`/`IntegerString` (default base 10). For a `List` of digits it Horner-folds `value = value*base + digit` in GMP; symbolic digits or a symbolic base produce a polynomial in the base via the evaluator instead. For a `String` argument each character is mapped to a digit value in `[0, 36)` and folded with an integer base (`FromDigits["abc", b]`); a symbolic/non-integer base over a string is left unevaluated. Validates arity (`FromDigits::argb`) and integer base `>= 2` (`::ibase`).
