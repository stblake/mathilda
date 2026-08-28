---
references:
  - "G. E. Andrews, *The Theory of Partitions*, Cambridge University Press, 1998."
  - "H. Rademacher, \"On the partition function p(n)\", *Proc. London Math. Soc.* 43 (1937), 241–254 — the exact convergent series."
---

### The partition function

`PartitionsP[n]` counts the partitions of `n`. It grows super-polynomially — Hardy and
Ramanujan showed `p(n) ~ exp(π√(2n/3)) / (4n√3)` — so it cannot be found by enumeration at
large `n`. Mathilda evaluates the **Hardy–Ramanujan–Rademacher** convergent series, an exact
sum of analytic terms rounded to the nearest integer, giving `p(1000)` in an instant.
Ramanujan's congruences hold: `p(5k+4) ≡ 0 (mod 5)`, `p(7k+5) ≡ 0 (mod 7)`, `p(11k+6) ≡ 0
(mod 11)`.

### Worked examples

```mathematica
In[1]:= PartitionsP[10]
Out[1]= 42
```

```mathematica
In[1]:= PartitionsP[100]
Out[1]= 190569292
```

```mathematica
In[1]:= Table[Mod[PartitionsP[5 k + 4], 5], {k, 0, 5}]
Out[1]= {0, 0, 0, 0, 0, 0}
```

The last example exhibits Ramanujan's congruence `p(5k+4) ≡ 0 (mod 5)`: every value is
divisible by `5`.
