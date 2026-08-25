---
status: Stable
---
### Worked examples

```mathematica
In[1]:= QuotientRemainder[17, 5]
Out[1]= {3, 2}
```

```mathematica
In[1]:= QuotientRemainder[-17, 5]
Out[1]= {-4, 3}
```

```mathematica
In[1]:= QuotientRemainder[17, -5]
Out[1]= {-4, -3}
```

### Notes

`QuotientRemainder[m, n]` returns the quotient and remainder together as
`{Quotient[m, n], Mod[m, n]}`. Because the quotient is floored and the remainder
takes the sign of the divisor `n`, the two always reconstruct the dividend:
`n q + r == m`. So `QuotientRemainder[-17, 5] = {-4, 3}` (a non-negative
remainder) while `QuotientRemainder[17, -5] = {-4, -3}` (a non-positive one).
