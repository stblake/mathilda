---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on arbitrary-precision integer multiplication."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on canonical forms for products."
---
### Worked examples

```mathematica
In[1]:= 2^64 * 3
Out[1]= 55340232221128654848
```

```mathematica
In[1]:= (1/3)*(3/7)*7
Out[1]= 1
```

```mathematica
In[1]:= x y x z
Out[1]= x^2 y z
```

### Notes

Times shares the `Flat`, `Orderless`, and `OneIdentity` attributes with Plus, so
factors are flattened and canonically ordered. Repeated symbolic factors are
gathered into powers, turning `x y x z` into `x^2 y z`. Rational factors are
multiplied and reduced to lowest terms, and a product of reciprocals can cancel
exactly to `1`. Integer products promote to GMP bigints automatically when they
exceed machine-word range.
