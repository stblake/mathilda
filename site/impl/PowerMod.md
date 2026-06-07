---
references:
  - "A. Tonelli, \"Bemerkung über die Auflösung quadratischer Congruenzen\", Göttinger Nachrichten, 1891; D. Shanks, \"Five number-theoretic algorithms\", 1973."
source: src/numbertheory.c
---
**Algorithm.** `builtin_powermod` computes `PowerMod[a, b, m]` on integer-like `a`, `m` (any sign on `m`; reduced to `|m|`). Two cases.

*Integer/BigInt exponent (GMP fast path).* For `b ≥ 0`, `mpz_powm(r, a, b, m)`. For `b < 0`, it first inverts `a` modulo `m` with `mpz_invert` (returning NULL/unevaluated if no inverse exists, i.e. `gcd(a,m) ≠ 1`), then raises the inverse to `-b`.

*Rational exponent `p/q` (modular root).* This asks for `x` with `x^q ≡ a^p (mod m)`. It first forms `c = a^p mod m` (inverting `a` when `p < 0`), then calls `modular_root(root, c, q, m)`, since GMP has no primitive modular r-th root. `modular_root` (1) brute-forces for tiny moduli (`m ≤ 1000000`, `modroot_brute`); (2) otherwise factors `m` via `internal_factorinteger`, (3) for each prime power `p^e` solves `x_0^r ≡ c (mod p)` — Tonelli-Shanks when `r = 2` (`tonelli_shanks`, with the `p ≡ 3 mod 4` shortcut), the closed form `x = c^(r^{-1} mod p-1)` when `gcd(r, p-1) = 1`, or brute force for small primes — then Hensel-lifts `x_0` to `mod p^e`, and (4) combines the per-prime-power roots by CRT. Any unsupported/no-solution case returns 0, leaving the surface `PowerMod[...]` echoed back unevaluated.

**Data structures.** All-`mpz_t` GMP integers throughout; results normalised back to `EXPR_INTEGER`/`EXPR_BIGINT` via `expr_bigint_normalize`.

**Complexity / limits.** Integer case is GMP's modular exponentiation, `O(log b)` modular multiplies. The modular-root case is bounded by `FactorInteger` on `m` plus Tonelli-Shanks (`O(log^2 p)` per prime) and Hensel lifting; the brute-force fast path is capped at `m ≤ 10^6`.
