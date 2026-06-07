---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on arbitrary-precision integer addition."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on normal forms for sums."
---
### Worked examples

```mathematica
In[1]:= 2^100 + 3^50
Out[1]= 1267651318126217093349291975625
```

```mathematica
In[1]:= 1/2 + 1/3 + 1/6
Out[1]= 1
```

```mathematica
In[1]:= a + a + b + 2 a
Out[1]= 4 a + b
```

### Notes

Plus is `Flat`, `Orderless`, and `OneIdentity`, so nested sums are flattened and
terms are sorted into canonical order before combination. Integer addition
auto-promotes to GMP bigints on overflow, as in `2^100 + 3^50`. Exact rationals
are added with a common denominator and the result is reduced to lowest terms,
collapsing to an integer when the denominator divides out. Like terms are
collected by their symbolic factor, so `a + a + 2 a` becomes `4 a`.
