---
status: Stable
references:
  - "Richard Bird, *Introduction to Functional Programming using Haskell*, 2nd ed. (the foldl/reduce accumulator pattern)."
---
### Worked examples

```mathematica
In[1]:= Fold[f, x, {a, b, c}]
Out[1]= f[f[f[x, a], b], c]
```

```mathematica
In[1]:= Fold[Plus, 0, {1, 2, 3, 4}]
Out[1]= 10
```

```mathematica
In[1]:= Fold[#1*10 + #2 &, 0, {1, 2, 3}]
Out[1]= 123
```

### Notes

`Fold[f, x, list]` is a left fold: it threads an accumulator (initial value `x`)
through `list`, applying `f[acc, elem]` at each step. The classic `Fold[Plus, 0,
list]` sums a list, while `#1*10 + #2 &` shows the accumulator (`#1`) and current
element (`#2`) being combined to build a number digit by digit. The two-argument
form `Fold[f, list]` uses the list's first element as the seed. `Fold[f, x, {}]`
returns `x` unchanged.
