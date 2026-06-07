---
status: Stable
references:
  - "Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.3 (sequences as conventional interfaces; filtering)."
---
### Worked examples

```mathematica
In[1]:= Select[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {2, 4, 6}
```

```mathematica
In[1]:= Select[Range[10], # > 5 &]
Out[1]= {6, 7, 8, 9, 10}
```

```mathematica
In[1]:= Select[{1, 2, 3, 4, 5}, PrimeQ, 2]
Out[1]= {2, 3}
```

### Notes

`Select[list, crit]` keeps the elements for which `crit[elem]` returns `True`;
any other result (including `False` or an unevaluated predicate) drops the
element. The criterion is usually a predicate symbol such as `EvenQ` or
`PrimeQ`, or a pure function like `# > 5 &`. The optional third argument caps the
number of elements returned, which lets `Select` stop early once enough matches
are found.
