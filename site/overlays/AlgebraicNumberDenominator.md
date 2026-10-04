### Worked examples

```mathematica
In[1]:= AlgebraicNumberDenominator[Sqrt[2]]  (* an algebraic integer: denominator 1 *)
```

```mathematica
In[1]:= AlgebraicNumberDenominator[1/5 + Sqrt[2]]  (* the smallest n making n a an algebraic integer *)
```

```mathematica
In[1]:= AlgebraicNumberDenominator[1/2]  (* a plain rational denominator *)
```

```mathematica
In[1]:= AlgebraicNumberDenominator[{1/5 + Sqrt[2], Sqrt[2], 1/2}]  (* threads over a list *)
```

### Notes

`AlgebraicNumberDenominator[a]` is the smallest positive integer `n` such that
`n a` is an algebraic integer. It is computed exactly, and is **not** simply the
leading coefficient of the minimal polynomial, which only bounds it from above:
`1/5 + Sqrt[2]` has minimal polynomial `25 x^2 - 10 x - 49` (leading coefficient
`25`), yet its denominator is `5`. The engine factors the leading coefficient
once and takes a per-prime valuation over the lower coefficients.

`Listable` and `Protected`; requires FLINT. For any algebraic integer the answer
is `1`.
