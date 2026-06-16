### Worked examples

```mathematica
In[1]:= Union[{3, 1, 2, 1, 3}]
Out[1]= {1, 2, 3}
```

With several arguments `Union` performs a true set union, merging and deduplicating across all lists:

```mathematica
In[1]:= Union[{1, 2, 3}, {2, 3, 4}, {5}]
Out[1]= {1, 2, 3, 4, 5}
```

Deduplication uses canonical structural equality, so symbolic and compound elements are handled too, and the result is returned in canonical sorted order:

```mathematica
In[1]:= Union[{x, Sin[y], x, 1, Sin[y]}]
Out[1]= {1, x, Sin[y]}
```

### Notes

`Union[list]` gives the sorted list of distinct elements; `Union[l1, l2, ...]` gives the set union. Comparison is by canonical structural equality, and the output is always sorted with duplicates removed.
