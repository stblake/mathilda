---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the prime-factor counting functions nu and Omega (Chapter 22)."
---

### Worked examples

```mathematica
In[1]:= PrimeNu[360]  (* Distinct prime factors {2, 3, 5}, regardless of exponents *)
```

```mathematica
In[1]:= PrimeNu[2^10]  (* A prime power has a single distinct factor *)
```

```mathematica
In[1]:= PrimeNu[1]  (* No prime factors *)
```

```mathematica
In[1]:= PrimeNu[Range[10]]  (* Threaded element-wise over a range by Listable *)
```

```mathematica
In[1]:= PrimeNu[72]  (* Equal to PrimeOmega only when n is square-free; 72 = 2^3 * 3^2 is not *)
```

### Notes

`PrimeNu[n]` is the number of *distinct* prime factors of `n`, independent of their
multiplicities. Its companion [`PrimeOmega`](PrimeOmega.md) counts prime factors *with*
multiplicity, so `PrimeNu[n] <= PrimeOmega[n]`, with equality exactly when `n` is
square-free (the condition [`MoebiusMu`](MoebiusMu.md) tests). The sign of `n` is ignored
and `PrimeNu[1] = 0`.

A non-real Gaussian-integer argument, or `GaussianIntegers -> True`, counts distinct Gaussian
prime factors over `Z[i]`. `PrimeNu[0]` is left unevaluated.
