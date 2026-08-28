---
references:
  - "R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, 2nd ed., Springer, 2005 — Miller–Rabin (§3.5) and the Baillie–PSW test."
  - "R. Baillie and S. S. Wagstaff Jr., \"Lucas pseudoprimes\", *Math. Comp.* 35 (1980), 1391–1417."
---

### Primality testing

`PrimeQ` does not trial-divide. It runs a strong probabilistic test — a Baillie–PSW
pre-screen (a strong Fermat test composed with a Lucas test) backed by Miller–Rabin rounds
— for which no composite counterexample is known, and which is deterministic below `2^64`.
Deciding that a large number is *composite* is thus fast; *factoring* it (see
[`FactorInteger`](FactorInteger.md)) is the hard problem behind RSA.

### Worked examples

```mathematica
In[1]:= PrimeQ[97]
Out[1]= True
```

Mersenne numbers `2^p - 1` are handled instantly; `2^31 - 1` is the prime 2147483647, while `2^67 - 1` is composite (a famous factorisation by Frank Nelson Cole):

```mathematica
In[1]:= PrimeQ[2^31 - 1]
Out[1]= True

In[2]:= PrimeQ[2^67 - 1]
Out[2]= False
```

Carmichael numbers fool the naive Fermat test but not `PrimeQ`; 561 is correctly reported composite:

```mathematica
In[1]:= PrimeQ[561]
Out[1]= False
```

`GaussianIntegers -> True` tests primality in `Z[i]`. A rational prime `p ≡ 1 (mod 4)` splits and is *not* a Gaussian prime, whereas `p ≡ 3 (mod 4)` remains prime:

```mathematica
In[1]:= PrimeQ[5, GaussianIntegers -> True]
Out[1]= False

In[2]:= PrimeQ[3, GaussianIntegers -> True]
Out[2]= True
```

A Gaussian-integer argument is tested directly; `2 + 3 I` has norm 13 and is a Gaussian prime:

```mathematica
In[1]:= PrimeQ[2 + 3 I]
Out[1]= True
```

### Notes

`PrimeQ[n]` tests primality with GMP's `mpz_probab_prime_p` (25 Miller-Rabin
rounds atop a Baillie-PSW pre-screen), so it is definitive for `n < 2^64` and
has false-positive probability below `4^-25` otherwise. It is not deceived by
Carmichael numbers such as 561. With `GaussianIntegers -> True`, or when given a
Gaussian integer `a + b I`, primality is decided in the ring `Z[i]`: rational
primes `≡ 1 (mod 4)` factor and are reported composite, while those `≡ 3 (mod 4)`
stay prime.
