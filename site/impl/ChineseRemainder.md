---
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", §5.4."
source: src/numbertheory/chineseremainder.c
---
**Algorithm.** `builtin_chineseremainder` reconstructs the integer `x` from a list of
residues `{r1, ..., rk}` and moduli `{m1, ..., mk}` by a **streaming
extended-Euclidean CRT fold**. It carries a running pair `(x, M)` -- the current
solution and the LCM of the moduli seen so far -- initialised from the first
congruence, and folds in each subsequent `(ri, mi)` in turn.

*The fold step.* To merge `x ≡ a (mod M)` with `x ≡ ri (mod mi)`, it runs the extended
Euclidean algorithm on `M` and `mi` to obtain `g = GCD[M, mi]` and Bezout cofactors.
The merged congruence is solvable iff `g` divides `ri - a`; if it does not, the system
is inconsistent and the builtin returns `NULL`, leaving the `ChineseRemainder[...]`
surface unevaluated. When it is solvable, the update sets the new modulus to
`LCM[M, mi] = M*mi/g` and lifts `x` to the unique residue modulo that LCM. Handling the
general `g > 1` case here is exactly what lets the moduli be non-coprime.

*Range normalisation.* The final `x` is reduced into `0 <= x < M` (or, given the
optional offset `d`, into `d <= x < d + M` by adding the appropriate multiple of `M`).

**Data structures.** All arithmetic is on `mpz_t` GMP integers, so residues and moduli
of arbitrary size are exact; results are normalised back to `EXPR_INTEGER` /
`EXPR_BIGINT` via `expr_bigint_normalize`. Inputs are validated to be two equal-length
lists of integer-like values (plus an optional integer offset) before the fold begins.

**Complexity.** Linear in the number of congruences, each fold step costing one extended
GCD -- `O(k * M(n) log n)` for `k` congruences of `n`-bit moduli, dominated by GMP's
sub-quadratic GCD. No factorisation of the moduli is required.
