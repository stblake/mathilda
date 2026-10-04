---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the prime-factor counting functions nu and Omega (Chapter 22)."
---

### Worked examples

```mathematica
In[1]:= PrimeOmega[360]  (* Counted with multiplicity: 360 = 2^3 * 3^2 * 5 gives 3 + 2 + 1 *)
```

```mathematica
In[1]:= PrimeOmega[2^10]  (* A prime power contributes its full exponent *)
```

```mathematica
In[1]:= PrimeOmega[1]  (* The empty factorisation gives 0 *)
```

```mathematica
In[1]:= PrimeOmega[Range[10]]  (* Threaded element-wise over a range by Listable *)
```

```mathematica
In[1]:= PrimeOmega[72]  (* Compare with PrimeNu[72] = 2; 72 = 2^3 * 3^2 is not square-free, so Omega is larger *)
```

### Notes

`PrimeOmega[n]` is the total number of prime factors of `n` counted with multiplicity — the
sum of the exponents in the prime factorisation. Its additive companion
[`PrimeNu`](PrimeNu.md) counts *distinct* primes, and the two agree exactly when `n` is
square-free; in general `PrimeOmega[n] >= PrimeNu[n]`. `PrimeOmega` is the exponent
[`LiouvilleLambda`](LiouvilleLambda.md) raises `-1` to, so `LiouvilleLambda[n] ==
(-1)^PrimeOmega[n]`.

The sign of `n` is ignored and `PrimeOmega[1] = 0`. A non-real Gaussian-integer argument, or
`GaussianIntegers -> True`, counts Gaussian prime factors over `Z[i]`. `PrimeOmega[0]` is
left unevaluated.
