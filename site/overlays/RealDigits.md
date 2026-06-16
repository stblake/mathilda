### Worked examples

```mathematica
In[1]:= RealDigits[123.456]
Out[1]= {{1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 3}
```

For an exact rational with a non-terminating decimal expansion, the digit list
ends in a *nested* list giving the recurring block — here the period-6 cycle of
`1/7`, and the mixed `22/7 = 3.142857142857…`:

```mathematica
In[1]:= RealDigits[1/7]
Out[1]= {{{1, 4, 2, 8, 5, 7}}, 0}

In[2]:= RealDigits[22/7]
Out[2]= {{3, {1, 4, 2, 8, 5, 7}}, 1}
```

High-precision constants expose their digits directly. Thirty digits of `π`
from a 30-digit MPFR value, and the first ten significant base-10 digits of a
40-digit `π`:

```mathematica
In[1]:= RealDigits[N[Pi, 30]]
Out[1]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3, 8, 3, 2, 8}, 1}

In[2]:= RealDigits[N[Pi, 40], 10, 10]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3}, 1}
```

`RealDigits` also works in other bases — the binary expansion of `255` is eight
ones:

```mathematica
In[1]:= RealDigits[255, 2]
Out[1]= {{1, 1, 1, 1, 1, 1, 1, 1}, 8}
```

### Notes

`RealDigits[x]` gives `{digits, exp}` where the first digit is the coefficient of
`10^(exp - 1)`. `RealDigits[x, b]` uses base `b`; `RealDigits[x, b, len]` returns
`len` digits; `RealDigits[x, b, len, n]` starts from the coefficient of `b^n`.
For rationals with non-terminating expansions the digit list ends in a nested
list of the recurring block. For inexact reals, digits beyond the available
precision are returned as `Indeterminate`. The sign of `x` is discarded.
