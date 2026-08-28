---
references:
  - "A. J. Menezes, P. C. van Oorschot and S. A. Vanstone, *Handbook of Applied Cryptography*, CRC Press, 1996 — §3.6 covers the discrete logarithm problem and its algorithms (baby-step giant-step, Pohlig–Hellman, index calculus)."
  - "R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, 2nd ed., Springer, 2005 — §5.2, discrete logarithms."
---

### Worked examples

```mathematica
In[1]:= MultiplicativeOrder[2, 7]
Out[1]= 3
```

```mathematica
In[1]:= MultiplicativeOrder[10, 7]
Out[1]= 6
```

```mathematica
In[1]:= MultiplicativeOrder[7, 1000000007]
Out[1]= 500000003
```

```mathematica
In[1]:= MultiplicativeOrder[3, 998244353]
Out[1]= 998244352
```

```mathematica
In[1]:= MultiplicativeOrder[2, 11, {1, 10}]
Out[1]= 5
```

### The discrete logarithm problem

The multiplicative order of `k` modulo `n` is the size of the cyclic subgroup that `k`
generates in the units mod `n`. It is intimately tied to the **discrete logarithm problem**
(DLP): given a base `k` and a target `b`, find the exponent `x` with `k^x ≡ b (mod n)`. The
three-argument form `MultiplicativeOrder[k, n, {r1, ...}]` returns the least exponent whose
power lands in the residue set, so it *is* a (multi-target) discrete logarithm:

```mathematica
In[1]:= MultiplicativeOrder[3, 7, {5}]
Out[1]= 5
```

Here `3` is a primitive root of `7`, and the least `x` with `3^x ≡ 5 (mod 7)` is `5` — the
discrete logarithm of `5` to base `3`. Computing the *order* is easy: it divides `EulerPhi[n]`
and is recovered by stripping prime factors from it. Inverting it — the DLP itself — is
believed **hard**: no polynomial-time algorithm is known over a general prime field, and that
presumed hardness is the foundation of Diffie–Hellman key exchange and the ElGamal and DSA
signature schemes. The standard algorithms (baby-step giant-step, Pollard's rho for
logarithms, Pohlig–Hellman, and index calculus) are surveyed in the references below. Mathilda's
three-argument form simply walks the orbit `k^1, k^2, …`, so it is a teaching tool for small
moduli, not a practical logarithm at cryptographic sizes.

### Notes

`MultiplicativeOrder[k, n]` is the smallest `m > 0` with `k^m ≡ 1 (mod n)`. The
order `6` for `10` modulo `7` reflects that `1/7 = 0.142857...` has a repeating
block of length `6`. The two large-modulus cases use prime moduli: `3` is a
primitive root of the NTT prime `998244353`, so its order equals `n - 1`. The
three-argument form `MultiplicativeOrder[k, n, {r1, ...}]` finds the least `m`
with `k^m` congruent to one of the listed residues. All arithmetic is exact via
GMP. The result is unevaluated when `gcd(k, n) ≠ 1`.
