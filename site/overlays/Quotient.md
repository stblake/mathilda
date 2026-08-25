---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on the division algorithm."
---
### Worked examples

```mathematica
In[1]:= Quotient[17, 5]
Out[1]= 3
```

```mathematica
In[1]:= Quotient[-17, 5]
Out[1]= -4
```

```mathematica
In[1]:= 5 Quotient[-17, 5] + Mod[-17, 5]
Out[1]= -17
```

```mathematica
In[1]:= Quotient[5 + 3 I, 2]
Out[1]= 2 + 2 I
```

### Notes

`Quotient[m, n]` floors the ratio toward `-Infinity`, so `Quotient[-17, 5] = -4`
(not `-3`) and the division identity `n Quotient[m, n] + Mod[m, n] == m` holds
exactly. For complex arguments it is Gaussian-integer division, rounding each part
of the ratio to the nearest integer, so `Quotient[5 + 3 I, 2] = 2 + 2 I` — the
Gaussian integer nearest the ratio `2.5 + 1.5 I`. The three-argument
`Quotient[m, n, d]` uses the same offset convention as the three-argument `Mod`.
