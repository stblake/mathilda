---
references:
  - "A. Ya. Khinchin, *Continued Fractions*, Dover, 1997 — the classical theory of convergents and best approximations."
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — Chapter X, continued fractions."
---

### Convergents and periodicity

The partial quotients of `ContinuedFraction[p/q]` are exactly the quotients of the Euclidean
algorithm on `p` and `q`, so a rational has a *finite* expansion. Truncating the expansion
and folding it back (`FromContinuedFraction`) gives the *convergents* — the best rational
approximations to `x` for their size, which is why `{3, 7, 15, 1}` recovers `355/113 ≈ π`. A
quadratic irrational such as `Sqrt[7]` has an *eventually periodic* expansion (Lagrange's
theorem), returned in the no-count form as a bracketed repeating block.

### Worked examples

```mathematica
In[1]:= ContinuedFraction[123/47]
Out[1]= {2, 1, 1, 1, 1, 1, 1, 3}
```

```mathematica
In[1]:= ContinuedFraction[Sqrt[2]]
Out[1]= {1, {2}}

In[2]:= ContinuedFraction[Sqrt[7]]
Out[2]= {2, {1, 1, 1, 4}}
```

```mathematica
In[1]:= ContinuedFraction[N[Pi, 40], 12]
Out[1]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1}
```

```mathematica
In[1]:= FromContinuedFraction[{1, 2, 2, 2, 2, 2}]
Out[1]= 99/70
```

### Notes

For an exact rational the expansion is finite and canonical (last term `>= 2`). For `Sqrt[d]` with `d` a non-square integer the no-count form returns the periodic block in braces, `{a1, {b1, ...}}`. Inexact inputs yield only as many terms as the precision determines; the `Pi` example exposes the famous large term `292`. `ContinuedFraction` is Listable.
