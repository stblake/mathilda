---
source: src/power.c
---
**Algorithm.** `builtin_power` evaluates `Power[base, exp]`. `Power[x]` is x; `Power[b, e1, e2, ...]` is right-associated into `Power[b, Power[e1, e2, ...]]` (right-associative grouping). The two-argument core handles, in order: infinity/`Indeterminate` algebra (`0^Infinity -> 0`, `1^Infinity -> Indeterminate` with message, `Infinity^n` by sign of n, etc.); numeric exact folding (integer/rational/bigint powers via GMP, e.g. exact `2^10`, `(1/2)^3`); inexact Real/MPFR exponentiation; partial radical simplification of `integer^(p/q)` (pulling out perfect-power factors so `Sqrt[8] -> 2 Sqrt[2]`); `(b^m)^n -> b^(m·n)` and product/zero/one identities; and `Sqrt`-style rational-exponent canonicalisation. `Sqrt[x]` is a thin wrapper (`builtin_sqrt`) that rewrites to `Power[x, 1/2]`. Symbolic cases that cannot be reduced return `NULL`, leaving the call unevaluated. `Power` is `ONEIDENTITY | LISTABLE | NUMERICFUNCTION | PROTECTED` (note: not Flat/Orderless — exponentiation is neither associative nor commutative).

**Data structures.** `Expr*` trees; exact integer/bigint exponentiation uses GMP `mpz`, rationals via `make_rational`, and MPFR for high-precision reals. Radical factor extraction works on integer factorisation of the base.

**Complexity / limits.** Integer powers are `O(log exp)` GMP multiplies; radical canonicalisation costs a factorisation of the integer base.
