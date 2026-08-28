---
references:
  - "D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, 1997 — the extended Euclidean algorithm and Bézout's identity (§4.5.2)."
  - "H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer, 1993 — extended GCD and modular inversion."
---

### Bézout's identity

`ExtendedGCD[a, b]` returns `{g, {s, t}}` with `s a + t b = g = GCD[a, b]` — a *Bézout
certificate*. The coefficients `s, t` are exactly what invert one integer modulo another
(`s` is `a`'s inverse mod `b` when `g = 1`), which is why the extended Euclidean algorithm
sits underneath `PowerMod[a, -1, m]` and the Chinese Remainder Theorem.

### Worked examples

```mathematica
In[1]:= ExtendedGCD[12, 18]
Out[1]= {6, {-1, 1}}

In[2]:= ExtendedGCD[15, 25, 35]
Out[2]= {5, {2, -1, 0}}
```

```mathematica
In[1]:= ExtendedGCD[17, 100]
Out[1]= {1, {-47, 8}}

In[2]:= PowerMod[17, -1, 100]
Out[2]= 53
```

```mathematica
In[1]:= ExtendedGCD[2^64, 3^40]
Out[1]= {1, {3997565229372176830, -6065478849745282079}}
```

### Notes

`ExtendedGCD[n1, ...]` returns `{g, {r1, r2, ...}}` where `g == GCD[n1, ...]` and `g == r1 n1 + r2 n2 + ...`. For Out[1], `6 == (-1)(12) + (1)(18)`; the multi-argument form folds `mpz_gcdext` pairwise. The Bezout coefficient gives the modular inverse: since `17 (-47) == 1 (mod 100)` and `-47 == 53 (mod 100)`, it matches `PowerMod[17, -1, 100]`. BigInt inputs such as `2^64` and `3^40` are handled exactly.
