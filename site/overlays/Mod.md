---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on the division algorithm and modular reduction."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on residues and modular arithmetic."
---
### Worked examples

```mathematica
In[1]:= Mod[-17, 5]
Out[1]= 3
```

```mathematica
In[1]:= Mod[17, -5]
Out[1]= -3
```

```mathematica
In[1]:= Mod[2^100, 7]
Out[1]= 2
```

```mathematica
In[1]:= Mod[10, 3, 1]
Out[1]= 1
```

### Notes

`Mod[m, n]` returns a residue that takes the sign of the divisor `n`, so
`Mod[-17, 5] = 3` (non-negative) while `Mod[17, -5] = -3` (non-positive). This is
the mathematician's floored modulus, not C's truncated `%`. Reduction is exact on
bigints, so `Mod[2^100, 7] = 2` without overflow. The three-argument form
`Mod[m, n, d]` returns the representative in the offset range `[d, d+n)`, e.g.
`Mod[10, 3, 1] = 1`.
