---
references:
  - "K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — the Legendre and Jacobi symbols and quadratic reciprocity (Chapter 5)."
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008."
---

### Quadratic residues and reciprocity

For an odd prime `m`, the Jacobi symbol `(n/m)` is the *Legendre symbol*: `+1` if `n` is a
non-zero quadratic residue modulo `m`, `-1` if it is a non-residue, and `0` if `m ∣ n`. By
*Euler's criterion*, `(n/m) ≡ n^((m-1)/2) (mod m)`. For composite (odd) `m` the Jacobi
symbol is the product of the Legendre symbols over the prime factors, and it obeys the *law
of quadratic reciprocity*, which is what makes it computable in `O(log² n)` steps without
factoring `m` — the same recursion the [`PowerMod`](PowerMod.md) modular square root relies
on. Mathilda returns the full Kronecker generalisation, so `m` may be even or non-positive.

### Worked examples

```mathematica
In[1]:= JacobiSymbol[2, 7]
Out[1]= 1
```

```mathematica
In[1]:= JacobiSymbol[3, 7]
Out[1]= -1
```

`(2/7) = +1` because `2 ≡ 3² (mod 7)` is a square; `(3/7) = -1` because `3` is a non-residue
modulo `7`.
