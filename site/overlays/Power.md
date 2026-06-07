---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on binary exponentiation."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on simplification of radical powers."
---
### Worked examples

```mathematica
In[1]:= 2^200
Out[1]= 1606938044258990275541962092341162602522202993782792835301376
```

```mathematica
In[1]:= (1/2)^-5
Out[1]= 32
```

```mathematica
In[1]:= 27^(2/3)
Out[1]= 9
```

```mathematica
In[1]:= 0^0
Out[1]= Indeterminate
```

### Notes

Integer powers use binary exponentiation and promote to GMP bigints, so `2^200`
is exact. A rational base with a negative integer exponent inverts and raises,
giving `(1/2)^-5 = 32`. Rational exponents trigger perfect-power extraction:
`27^(2/3)` reduces to `9`, while non-extractable cases such as `8^(1/3)` of a
non-cube stay symbolic. The indeterminate form `0^0` evaluates to
`Indeterminate` rather than `1`.
