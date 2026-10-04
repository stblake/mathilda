### Worked examples

```mathematica
In[1]:= CountDistinct[{a, b, a, c, b, a}]
```

```mathematica
In[1]:= CountDistinct[{1, 1, 2, 3, 3, 3, 4}]
```

```mathematica
In[1]:= CountDistinct[<|x -> 1, y -> 1, z -> 2|>]  (* over an association, the distinct values *)
```

### Notes

`CountDistinct[list]` is the number of distinct elements — `Length[Union[list]]`
computed in a single hash pass, without building the sorted set. Over an
association it counts the distinct *values*. It answers the "how many different
things are here?" question directly; use `CountDistinctBy` to count distinct
values of a function of each element.
