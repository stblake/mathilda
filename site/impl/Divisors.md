---
source: src/numbertheory/divisors.c
---
**Algorithm.** `builtin_divisors` separates the single positional argument from an optional
`GaussianIntegers` rule, then enumerates divisors from the prime factorisation. The ordinary
path (`divisors_ordinary` in `nt_gaussian.c`) factors `|n|`, walks the divisor lattice with
a mixed-radix exponent odometer (digit `i` ranges `0..e_i`, forming each divisor as a
product of prime powers), and `qsort`s the results into ascending order; `Divisors[1]` is
`{1}`. The Gaussian path (`divisors_gaussian`), used under `GaussianIntegers -> True` or for
a non-real input, returns one first-quadrant representative per associate class, sorted by
`(Re, Im)`.

**Data structures.** Divisors are built in a `mpz_t` array via `mpz_pow_ui`/`mpz_mul` over
the prime powers, then emitted as a `List` of `Integer`/`BigInt` (ordinary) or `Complex`
(Gaussian) `Expr`. There is no ND/packed/`Compile` kernel — the result is a
variable-length list, not a machine buffer; `Listable` threading over a list of arguments is
done by the evaluator.

**Complexity / limits.** Cost is factoring plus `O(d log d)` to sort the `d = prod_i (e_i +
1)` divisors. The divisor count is computed first, and the call is left unevaluated if it
overflows `size_t` (e.g. `Divisors[100!]` has ~10^28 divisors, intractable to materialise).
`Divisors[0]`, a non-integer `n`, and a factorisation whose bases cannot be confirmed prime
are all left unevaluated; `Divisors[]` emits `Divisors::argx`. The sign of `n` is ignored.
