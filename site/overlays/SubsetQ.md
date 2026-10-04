### Worked examples

```mathematica
In[1]:= SubsetQ[{1, 2, 3, 4}, {2, 4}]  (* is the second a subset of the first? *)
In[2]:= SubsetQ[{1, 2, 3}, {4}]
```

### Notes

`SubsetQ[a, b]` is `True` when every element of `b` occurs in `a`; the test is built
on a hash set, so it is O(|a| + |b|). Multiplicity is ignored. Lists and associations
(compared by value) mix freely, but two other expressions must share a head.
