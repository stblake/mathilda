---
source: src/arithmetic.c
---
`builtin_divide` (`src/arithmetic.c`) evaluates `Divide[num, den]`. It special-cases numerics: if either operand is a `Real`, it computes a `double` quotient (coercing integer/bigint operands via `mpz_get_d`), emitting `Power::infy` and returning `ComplexInfinity` on a zero denominator. If both are exact integers/rationals (`is_rational`), it forms the reduced rational `n1*d2 / d1*n2` via `make_rational`, handling `x/0 -> ComplexInfinity` and `0/0 -> Indeterminate` with the matching diagnostics. Otherwise it falls back to the symbolic canonical form `Times[num, Power[den, -1]]`, which the evaluator then simplifies. (Divide carries no Hold attributes, so arguments arrive pre-evaluated.)
