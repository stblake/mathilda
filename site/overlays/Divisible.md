---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — divisibility (Chapter I)."
---

### The divisibility relation

`Divisible[n, m]` tests the relation `m ∣ n` — whether `n` is an integer multiple of `m`,
equivalently whether `Mod[n, m] == 0`. It extends beyond the ordinary integers: over the
Gaussian integers `Z[i]`, `m ∣ n` when the quotient `n/m` is itself a Gaussian integer.

### Worked examples

```mathematica
In[1]:= Divisible[100, 4]
Out[1]= True
```

```mathematica
In[1]:= Divisible[100, 7]
Out[1]= False
```

```mathematica
In[1]:= Divisible[10 + 5 I, 1 + 2 I]
Out[1]= True
```

`4` divides `100` but `7` does not. The last example is a Gaussian divisibility:
`(10 + 5i)/(1 + 2i) = 4 - 3i` is a Gaussian integer, so `1 + 2i` divides `10 + 5i` in `Z[i]`.
