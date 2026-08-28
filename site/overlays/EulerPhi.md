---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the totient function and Euler's theorem."
---

### The totient and Euler's theorem

`EulerPhi[n]` counts the integers in `1, …, n` coprime to `n`, which is the *order* of the
group of units `(Z/nZ)*`. It is multiplicative, with `EulerPhi[n] = n ∏_{p|n} (1 - 1/p)`.
Its central role is *Euler's theorem*: `a^EulerPhi[n] ≡ 1 (mod n)` whenever `gcd(a, n) = 1` —
the identity that makes `PowerMod[a, -1, m]` and RSA decryption work.

### Worked examples

```mathematica
In[1]:= EulerPhi[36]
Out[1]= 12
```

```mathematica
In[1]:= Table[EulerPhi[n], {n, 1, 12}]
Out[1]= {1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4}
```

```mathematica
In[1]:= EulerPhi[2^61 - 1]
Out[1]= 2305843009213693950
```

```mathematica
In[1]:= Total[Map[EulerPhi, {1, 2, 3, 5, 6, 10, 15, 30}]]
Out[1]= 30
```

### Notes

`EulerPhi[n]` counts the integers in `1..n` coprime to `n`. The Mersenne
prime `2^61 - 1` is prime, so `phi = p - 1`. The last example is Gauss's
identity `Sum phi(d) = n` over the divisors `d` of `30`, recovering `30`
exactly.
