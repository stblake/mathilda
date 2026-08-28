---
references:
  - "K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — primitive roots and the structure of `(Z/nZ)*` (Chapter 4)."
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008."
---

### Generators of the units

A primitive root of `n` is a generator of the multiplicative group `(Z/nZ)*`: its powers run
through every residue coprime to `n`, so its `MultiplicativeOrder` equals `EulerPhi[n]`. Such
a generator exists exactly when the group is cyclic — for `n = 2, 4, p^k, 2p^k` — and how
often a *small* number is a primitive root is the subject of Artin's still-open conjecture.

### Worked examples

```mathematica
In[1]:= PrimitiveRoot[7]
Out[1]= 3
```

It finds generators even for large primes; 5 generates the multiplicative group modulo the prime `10^9 + 7`:

```mathematica
In[1]:= PrimitiveRoot[10^9 + 7]
Out[1]= 5
```

Prime powers are supported directly:

```mathematica
In[1]:= PrimitiveRoot[3^5]
Out[1]= 2
```

The two-argument form returns the smallest primitive root not below `k`:

```mathematica
In[1]:= PrimitiveRoot[7, 5]
Out[1]= 5
```

When the multiplicative group is non-cyclic (e.g. `n = 8`), no primitive root exists and the call stays unevaluated:

```mathematica
In[1]:= PrimitiveRoot[8]
Out[1]= PrimitiveRoot[8]
```

### Notes

`PrimitiveRoot[n]` returns a generator of the multiplicative group of integers
coprime to `n`. Such a generator exists only when `n` is `2`, `4`, an odd prime
power `p^k`, or twice one (`2 p^k`); for all other `n` the group is non-cyclic
and the call is left unevaluated. `PrimitiveRoot[n, k]` returns the smallest
primitive root that is `>= k`.
