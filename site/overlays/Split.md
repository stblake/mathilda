### Worked examples

```mathematica
In[1]:= Split[{1, 1, 2, 3, 3, 3, 1}]
Out[1]= {{1, 1}, {2}, {3, 3, 3}, {1}}
```

```mathematica
In[1]:= Split[Sort[{3, 1, 1, 2, 3, 2, 1}]]
Out[1]= {{1, 1, 1}, {2, 2}, {3, 3}}
```

Sorting first turns `Split` into a run-length / tally tool, grouping all equal elements together.

```mathematica
In[1]:= Split[{1, 2, 4, 7, 8, 10, 11}, (#2 - #1 == 1) &]
Out[1]= {{1, 2}, {4}, {7, 8}, {10, 11}}
```

With a custom test, `Split` extracts maximal runs of *consecutive* integers.

### Notes

`Split[list]` breaks `list` into runs of consecutive identical elements. `Split[list, test]` instead starts a new run whenever `test[ei, ej]` is not `True` for adjacent elements `ei, ej`, so a relational test like `#2 - #1 == 1` groups arithmetic runs and `Greater` would group strictly descending runs. The concatenation of the result always reproduces the original list.
