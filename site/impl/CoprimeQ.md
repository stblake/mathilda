---
source: src/numbertheory/coprimeq.c
---
**Algorithm.** `builtin_coprimeq` first scans its arguments, lifting an optional
`GaussianIntegers -> True|False` rule out of the numeric operands — because `CoprimeQ`
is `Orderless` the option may sit at any position — and treating any other rule, or a
non-Boolean option value, as malformed (result `False`). With zero operands it returns
`False`, with one operand `True` (there are no pairs). Over the ordinary integers it
extracts every operand with `expr_to_mpz` and tests each unordered pair with `mpz_gcd`,
declaring the set coprime iff every pairwise GCD is `1`. When `GaussianIntegers -> True`
is set, or any operand is an exact non-real Gaussian integer, it switches to `Z[i]`:
`gaussian_pair_coprime` runs the Gaussian Euclidean algorithm with round-to-nearest
division (`coprimeq_round_div` computes `q = floor((2 num + den)/(2 den))`), the remainder
norm strictly decreasing each step, and reports coprime iff the final GCD has norm `1`
(a unit). As a `*Q` predicate it always returns a Boolean, so any operand that is not a
manifest integer or Gaussian integer forces `False`.

**Data structures.** Operands stay as borrowed `Expr*`; each is coerced into a GMP `mpz_t`
(integer path) or a pair of `mpz_t` real/imaginary parts (Gaussian path, via
`coprimeq_to_gaussian`, which accepts a rational integer or a `Complex[a, b]` with
integer-like parts). The Gaussian Euclidean step runs entirely in `mpz_t` scratch
registers. There is no ND/packed/`Compile` kernel: `CoprimeQ` is a Boolean predicate, and
`Listable` threading over list arguments is performed by the evaluator before this builtin
runs.

**Complexity / limits.** The integer path is `O(p^2)` GCDs for `p` operands, each `mpz_gcd`
quasi-linear in the operand size, so machine integers and arbitrary-precision bigints are
handled uniformly. The Gaussian path's Euclidean algorithm always terminates because the
remainder norm decreases. Sign is ignored (GCD uses magnitudes); `CoprimeQ[]` is `False`,
`CoprimeQ[n]` is `True`, and rationals, reals, symbols and malformed options all yield
`False` rather than leaving the call unevaluated.
