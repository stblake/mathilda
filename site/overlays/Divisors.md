---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — divisors and perfect numbers (Chapter 16)."
---

### Worked examples

```mathematica
In[1]:= Divisors[28]  (* A perfect number equals the sum of its proper divisors *)
```

```mathematica
In[1]:= Divisors[60]  (* The twelve divisors, ascending *)
```

```mathematica
In[1]:= Total[Divisors[28]]  (* Summing all divisors of a perfect number gives twice the number *)
```

```mathematica
In[1]:= Length[Divisors[720]] == DivisorSigma[0, 720]  (* The number of divisors is sigma_0 *)
```

```mathematica
In[1]:= Divisors[{12, 15}]  (* Threaded over a list by Listable *)
```

```mathematica
In[1]:= Divisors[5, GaussianIntegers -> True]  (* Over Z[i]: one first-quadrant representative per associate class *)
```

### Notes

`Divisors[n]` returns the ascending list of positive integers dividing `n`, built from the
prime factorisation (the divisor lattice), so the sign of `n` is ignored and `Divisors[1]`
is `{1}`. The count of divisors is [`DivisorSigma`](DivisorSigma.md)`[0, n]` and their sum
is `DivisorSigma[1, n]`; a number is perfect when that sum is `2 n`, as with `28`.

`Divisors[n, GaussianIntegers -> True]`, or a non-real Gaussian-integer `n`, returns the
divisors in `Z[i]`, one representative per unit-associate class sorted by `(Re, Im)`.
`Divisors[0]`, a non-integer argument, and any call whose divisor count would overflow (such
as `Divisors[100!]`) are left unevaluated.
