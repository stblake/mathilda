# PowerMod

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PowerMod[a, b, m] gives a^b mod m.
PowerMod[a, -1, m] finds the modular inverse of a modulo m.
PowerMod[a, 1/r, m] finds a modular r-th root of a.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}

In[5]:= PowerMod[100, 1/2, 17 * 19 * 23]
Out[5]= 10

In[6]:= PowerMod[2, 1/2, 10^18 + 9]
Out[6]= 742174169206529574
```

## Implementation notes

**Algorithm.** `builtin_powermod` computes `PowerMod[a, b, m]` on integer-like `a`, `m` (any sign on `m`; reduced to `|m|`). Two cases.

*Integer/BigInt exponent (GMP fast path).* For `b ≥ 0`, `mpz_powm(r, a, b, m)`. For `b < 0`, it first inverts `a` modulo `m` with `mpz_invert` (returning NULL/unevaluated if no inverse exists, i.e. `gcd(a,m) ≠ 1`), then raises the inverse to `-b`.

*Rational exponent `p/q` (modular root).* This asks for `x` with `x^q ≡ a^p (mod m)`. It first forms `c = a^p mod m` (inverting `a` when `p < 0`), then calls `modular_root(root, c, q, m)`, since GMP has no primitive modular r-th root. `modular_root` (1) brute-forces for tiny moduli (`m ≤ 1000000`, `modroot_brute`); (2) otherwise factors `m` via `internal_factorinteger`, (3) for each prime power `p^e` solves `x_0^r ≡ c (mod p)` — Tonelli-Shanks when `r = 2` (`tonelli_shanks`, with the `p ≡ 3 mod 4` shortcut), the closed form `x = c^(r^{-1} mod p-1)` when `gcd(r, p-1) = 1`, or brute force for small primes — then Hensel-lifts `x_0` to `mod p^e`, and (4) combines the per-prime-power roots by CRT. Any unsupported/no-solution case returns 0, leaving the surface `PowerMod[...]` echoed back unevaluated.

**Data structures.** All-`mpz_t` GMP integers throughout; results normalised back to `EXPR_INTEGER`/`EXPR_BIGINT` via `expr_bigint_normalize`.

**Complexity / limits.** Integer case is GMP's modular exponentiation, `O(log b)` modular multiplies. The modular-root case is bounded by `FactorInteger` on `m` plus Tonelli-Shanks (`O(log^2 p)` per prime) and Hensel lifting; the brute-force fast path is capped at `m ≤ 10^6`.

- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Integer-exponent path uses GMP `mpz_powm` / `mpz_invert`; `a`, `b`, `m`

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on modular exponentiation and modular inverses.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on the extended Euclidean algorithm.
- A. Tonelli, "Bemerkung über die Auflösung quadratischer Congruenzen", Göttinger Nachrichten, 1891; D. Shanks, "Five number-theoretic algorithms", 1973.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PowerMod[3, 1000000, 7]
Out[1]= 4
```

```mathematica
In[1]:= PowerMod[2, 100, 101]
Out[1]= 1
```

```mathematica
In[1]:= PowerMod[7, -1, 11]
Out[1]= 8
```

```mathematica
In[1]:= PowerMod[2, -1, 4]
Out[1]= PowerMod[2, -1, 4]
```

### Notes

`PowerMod[a, b, m]` computes `a^b mod m` using binary modular exponentiation, so a
huge exponent like `PowerMod[3, 1000000, 7] = 4` is evaluated without forming the
full power. `PowerMod[2, 100, 101] = 1` illustrates Fermat's little theorem, since
`101` is prime and `2^100 ≡ 1 (mod 101)`. A negative exponent computes the modular
inverse via the extended Euclidean algorithm: `PowerMod[7, -1, 11] = 8` because
`7*8 = 56 ≡ 1 (mod 11)`. When the inverse does not exist (the base shares a factor
with the modulus, as with `2` and `4`), the call returns unevaluated rather than
producing a spurious result.
