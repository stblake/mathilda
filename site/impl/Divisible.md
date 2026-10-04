---
source: src/numbertheory/divisible.c
---
**Algorithm.** `builtin_divisible` requires exactly two arguments (otherwise it emits
`Divisible::argm`/`argt` and leaves the call unevaluated). When both `n` and `m` are
integer-like it answers directly with GMP's `mpz_divisible_p` — exact at any precision,
and honouring the convention that divisibility by `0` holds iff `n == 0`. Otherwise it
assembles the quotient `Times[n, Power[m, -1]]`, evaluates it with `eval_and_free`, and
returns `True` iff the result is an integer or a Gaussian integer (so `3 + I` is divisible
by `1 - I`, `3/2` by `1/2`, and `2 Pi` by `Pi/2`). If the quotient is non-integral and
both arguments are concrete numeric quantities it returns `False`; if either argument is
symbolic/non-numeric it returns `NULL`, leaving the call unevaluated so user rules and
pattern matching can apply.

**Data structures.** The integer path uses two `mpz_t`. The general path builds the
quotient as `Power`/`Times` `Expr` trees and runs them back through the evaluator;
`divisible_is_numeric_quantity` walks the result to decide numeric-vs-symbolic, recognising
exact numbers, the named constants (`Pi`, `E`, `EulerGamma`, …), `Complex`/`Rational`, and
any `NumericFunction` applied to numeric quantities. There is no ND/packed/`Compile` path —
`Divisible` returns a Boolean and is `Listable`, threaded by the evaluator before the
builtin runs.

**Complexity / limits.** The integer test is a single quasi-linear `mpz_divisible_p`, so
`Divisible[10^3000 + 1, 16001]` is exact and cheap; the general path costs one full
evaluation of the quotient. Sign is ignored through the multiple test, `Divisible[0, 0]`
is `True`, and a numeric but non-divisible pair is `False`, while symbolic arguments stay
unevaluated.
