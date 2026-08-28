---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — coprimality and the density `6/π²` of coprime pairs."
---

### Relatively prime integers

Two integers are *coprime* (relatively prime) when their greatest common divisor is `1` —
they share no prime factor. `CoprimeQ` tests this **pairwise**, so `CoprimeQ[n1, n2, ...]` is
`True` only when every pair is coprime. A classical density result: two integers drawn at
random are coprime with probability `6/π² = 1/ζ(2) ≈ 0.6079`.

### Worked examples

```mathematica
In[1]:= CoprimeQ[14, 15]
Out[1]= True
```

```mathematica
In[1]:= CoprimeQ[14, 21]
Out[1]= False
```

```mathematica
In[1]:= CoprimeQ[6, 35, 143]
Out[1]= True
```

`14` and `15` share no factor, so they are coprime; `14` and `21` share a `7`. The three
numbers `6, 35, 143 = 2·3, 5·7, 11·13` are built from disjoint primes, so they are pairwise
coprime. With `GaussianIntegers -> True` the test is carried out in `Z[i]`.
