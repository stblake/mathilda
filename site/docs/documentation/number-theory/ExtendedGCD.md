# ExtendedGCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExtendedGCD[n1, n2, ...] gives the extended GCD {g, {r1, r2, ...}}, where g == GCD[n1, ...] and g == r1 n1 + r2 n2 + ... -- Bezout's identity, the certificate that also yields modular inverses.`**

<details>
<summary>Notes</summary>

Computed by folding GMP's mpz\_gcdext pairwise; accepts machine and BigInt integers and threads over lists. Non-integer or inexact arguments leave ExtendedGCD unevaluated.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (1)

```mathematica
In[1]:= ExtendedGCD[3, {5, 15}]
Out[1]= {{1, {2, -1}}, {3, {1, 0}}}
```

### Applications (5)

```mathematica
In[2]:= ExtendedGCD[12, 18]
Out[2]= {6, {-1, 1}}

In[3]:= ExtendedGCD[15, 25, 35]
Out[3]= {5, {2, -1, 0}}

In[4]:= ExtendedGCD[17, 100]
Out[4]= {1, {-47, 8}}

In[5]:= PowerMod[17, -1, 100]
Out[5]= 53

In[6]:= ExtendedGCD[2^64, 3^40]
Out[6]= {1, {3997565229372176830, -6065478849745282079}}
```

## Implementation notes

**Algorithm.** `builtin_extendedgcd` computes a multi-argument Bézout identity, returning `{g, {r1, ..., rN}}` with `g = GCD[n1, ...]` and `g = sum r_i n_i`. It folds GMP's `mpz_gcdext` pairwise: at each step `gcd(running_g, a_i) = s*running_g + t*a_i`, every previously accumulated cofactor is scaled by `s` and the new cofactor `t` is appended. `running_g` starts at `0`, so the first step yields `|a0|` with cofactor `sign(a0)`, and GMP normalises `g` non-negative throughout — matching the conventional sign. `ExtendedGCD[]` is `{0, {}}`.

**Data structures.** A heap array of `count` `mpz_t` cofactor accumulators plus scalar `mpz_t` registers for `gcdext`; results pass through `expr_bigint_normalize`. Integer-only — machine ints auto-promote via `expr_to_mpz`. Inexact (Real/MPFR) arguments emit `ExtendedGCD::exact`; exact-but-non-integer rationals emit `ExtendedGCD::egcd`; both return `NULL`.

**Complexity / limits.** Dominated by the `gcdext` chain, `O(N · M(d) log d)` for `d`-digit inputs (`M` = GMP's multiplication cost).

- Integer-only; both machine integers and GMP bigints are accepted, and cofactors demote back to machine integers when they fit.
- Computed by folding GMP's `mpz_gcdext` pairwise: `gcd(running_g, n_i) = s·running_g + t·n_i`, scaling the accumulated cofactors by `s` and appending `t` each step. The running gcd stays non-negative, so `g` is exactly `GCD` and the sign convention is the conventional one.
- Threads element-wise over lists, e.g. `ExtendedGCD[3, {5, 15}]` → `{{1, {2, -1}}, {3, {1, 0}}}`.
- Edge cases: `ExtendedGCD[]` → `{0, {}}`; `ExtendedGCD[n]` → `{|n|, {±1}}`; zeros and negatives handled (`g` non-negative).
- Diagnostics: an inexact (Real) argument emits `ExtendedGCD::exact`; an exact non-integer (Rational) argument emits `ExtendedGCD::egcd`; symbolic arguments leave the call unevaluated silently.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [OneIdentity](../../expression-information/OneIdentity/), [GCD](../../number-theory/GCD/)

- D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, 1997 — the extended Euclidean algorithm and Bézout's identity (§4.5.2).
- H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer, 1993 — extended GCD and modular inversion.
- D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed. (Addison-Wesley, 1997), §4.5.2.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_extended_gcd.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extended_gcd.c)

## Notes & additional examples

### Bézout's identity

`ExtendedGCD[a, b]` returns `{g, {s, t}}` with `s a + t b = g = GCD[a, b]` — a *Bézout
certificate*. The coefficients `s, t` are exactly what invert one integer modulo another
(`s` is `a`'s inverse mod `b` when `g = 1`), which is why the extended Euclidean algorithm
sits underneath `PowerMod[a, -1, m]` and the Chinese Remainder Theorem.

### Notes

`ExtendedGCD[n1, ...]` returns `{g, {r1, r2, ...}}` where `g == GCD[n1, ...]` and `g == r1 n1 + r2 n2 + ...`. For Out[1], `6 == (-1)(12) + (1)(18)`; the multi-argument form folds `mpz_gcdext` pairwise. The Bezout coefficient gives the modular inverse: since `17 (-47) == 1 (mod 100)` and `-47 == 53 (mod 100)`, it matches `PowerMod[17, -1, 100]`. BigInt inputs such as `2^64` and `3^40` are handled exactly.
