### Worked examples

```mathematica
In[1]:= FromContinuedFraction[{3, 7, 15, 1, 292}]
Out[1]= 103993/33102
```

```mathematica
In[1]:= FromContinuedFraction[{1, {1}}]
Out[1]= 1/2 (1 + Sqrt[5])
```

```mathematica
In[1]:= FromContinuedFraction[{1, 2, {2}}]
Out[1]= Sqrt[2]
```

```mathematica
In[1]:= FromContinuedFraction[{a, b, c}]
Out[1]= (a + (1 + a b) c)/(1 + b c)
```

### Notes

For a finite list, `FromContinuedFraction` reconstructs the rational (or
symbolic) convergent `a1 + 1/(a2 + 1/(a3 + ...))`. The classic terms
`{3, 7, 15, 1, 292}` give `103993/33102`, the celebrated convergent of Pi.

A trailing sublist marks the periodic part of a *quadratic irrational*: the
purely periodic `{1, {1}}` returns the golden ratio `(1 + Sqrt[5])/2`, and
`{1, 2, {2}}` recovers `Sqrt[2]` from its eventually-periodic expansion
`[1; 2, 2, 2, ...]`. With symbolic terms the result is the exact nested form,
left un-expanded. `FromContinuedFraction` is the inverse of `ContinuedFraction`.
