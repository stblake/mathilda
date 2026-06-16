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

A Gaussian-integer base is raised exactly, keeping the real and imaginary parts as integers:

```mathematica
In[1]:= (3 + 4 I)^10
Out[1]= -9653287 + 1476984*I
```

A negative radicand is split into its real radical and the principal complex unit:

```mathematica
In[1]:= Sqrt[-12]
Out[1]= (2*I) Sqrt[3]
```

Irrational powers numericalise to arbitrary precision on request:

```mathematica
In[1]:= N[2^(1/2), 40]
Out[1]= 1.4142135623730950488016887242096980785697
```

### Notes

Integer powers use binary exponentiation and promote to GMP bigints, so `2^200`
is exact. A rational base with a negative integer exponent inverts and raises,
giving `(1/2)^-5 = 32`. Rational exponents trigger perfect-power extraction:
`27^(2/3)` reduces to `9`, while non-extractable cases such as `8^(1/3)` of a
non-cube stay symbolic. The indeterminate form `0^0` evaluates to
`Indeterminate` rather than `1`. Complex bases (Gaussian integers, negative
radicands) are handled in closed form, and irrational powers of numeric bases
evaluate to the requested precision under `N[...]`.
