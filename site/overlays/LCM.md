---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on the Euclidean algorithm and least common multiples."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on GCD/LCM relations."
---
### Worked examples

```mathematica
In[1]:= LCM[3, 4, 5]
Out[1]= 60
```

```mathematica
In[1]:= LCM[12, 18, 30]
Out[1]= 180
```

```mathematica
In[1]:= LCM[4, 6]
Out[1]= 12
```

```mathematica
In[1]:= LCM[0, 5]
Out[1]= 0
```

### Notes

LCM uses the identity `lcm(a, b) = a*b / gcd(a, b)` and folds across all
arguments, so `LCM[3, 4, 5]` gives `60` and `LCM[12, 18, 30]` gives `180`. The
absorbing convention `LCM[0, n] = 0` holds, matching the fact that zero is the
only common multiple involving zero. Pairwise reduction keeps intermediate values
small, and results promote to GMP bigints when they exceed machine-word range.
